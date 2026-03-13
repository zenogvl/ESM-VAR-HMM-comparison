# ============================================================
# Multilevel VAR Pipeline using mlVAR
# ============================================================
# Structure:
#   1. kalman_filter()         - Impute missing values per variable via Kalman smoother
#   2. prepare_mlvar_data()    - Train/test split with optional Kalman imputation
#   3. fit_mlvar()             - Fit mlVAR model on training data
#   4. predict_mlvar()         - Predict on test data, compute per-variable RMSE
#   5. run_mlvar_pipeline()    - Outer loop over multiple datasets
#
# Design principles (to match mHMM pipeline and support 40 datasets):
#   - All column names (ID, day, beep) are passed as arguments, not assumed
#   - No night rows are inserted: mlVAR handles overnight breaks natively
#     via the dayvar/beepvar arguments, so no manual row insertion is needed
#   - Kalman imputation is applied per individual on the raw data before
#     the train/test split
#   - Data is NOT re-scaled inside this pipeline: your data is already on [0,1]
#     and mlVAR's scale=FALSE preserves that. See note on centering below.
#   - The model is fitted on train data; predictions use the subject-level
#     temporal coefficients extracted with mlVAR::getNet()
#
# ── Centering note ──────────────────────────────────────────
# mlVAR internally person-mean-centres predictors before estimating
# temporal effects when scale=TRUE (which also standardises to SD=1).
# Since your data is already rescaled to [0,1] you do NOT want SD
# rescaling, but person-mean centering is still advisable for VAR
# models because it removes stable between-person differences and
# ensures the lag-1 coefficients reflect within-person dynamics.
# We therefore pass scale=FALSE to mlVAR and centre predictors
# manually inside predict_mlvar() — exactly mirroring what the
# ResAnalysis function in your estimate_ml_var.R file does.
# ────────────────────────────────────────────────────────────
#
# ── Computation time note ───────────────────────────────────
# mlVAR uses lme4 under the hood (one mixed model per variable).
# For datasets with many variables or participants this is fast
# compared to MCMC. The bottleneck is usually the number of
# variables × the complexity of the random-effects structure.
# Setting contemporaneous = "correlated" and temporal = "correlated"
# is most thorough but slowest; use temporal = "orthogonal" if
# computation time becomes a concern.
# ────────────────────────────────────────────────────────────
# ============================================================

library(mlVAR)
library(dplyr)

# ============================================================
# 1. kalman_filter()
# ============================================================
# Applies a local-level (random walk + noise) Kalman smoother to
# a data frame of time series to impute missing values.
# Each column is treated as an independent univariate series.
#
# This is the same approach used in estimate_vars.R (fit_var_models).
# The function is kept simple and self-contained so it can be swapped
# for a multivariate implementation later if desired.
#
# Arguments:
#   x   - data frame or matrix with one column per variable,
#         rows ordered by time within a single individual
#
# Returns: data frame with NAs replaced by Kalman-smoothed estimates

kalman_filter <- function(x) {
  x <- as.data.frame(x)
  for (col in colnames(x)) {
    ts_col <- x[[col]]
    if (any(is.na(ts_col)) && !all(is.na(ts_col))) {
      # StructTS requires the first value to be non-NA.
      # Forward-fill leading NAs with the first observed value,
      # and backward-fill trailing NAs with the last observed value,
      # so the smoother has valid anchors at both ends.
      first_obs <- which(!is.na(ts_col))[1]
      last_obs  <- tail(which(!is.na(ts_col)), 1)
      if (first_obs > 1)
        ts_col[seq_len(first_obs - 1)] <- ts_col[first_obs]
      if (last_obs < length(ts_col))
        ts_col[seq(last_obs + 1, length(ts_col))] <- ts_col[last_obs]

      fit      <- stats::StructTS(ts_col, type = "level")
      smoothed <- stats::tsSmooth(fit)
      x[[col]] <- as.numeric(smoothed[, 1])
    }
  }
  return(x)
}


# ============================================================
# 2. prepare_mlvar_data()
# ============================================================
# For each individual:
#   1. Optionally imputes missing values via Kalman smoother
#   2. Splits into train (first train_proportion rows) and test (remainder)
#
# No night rows are inserted: mlVAR handles overnight breaks natively
# via dayvar/beepvar, and predict_mlvar() uses the day column directly
# to avoid predicting across night boundaries.
#
# Individuals with fewer than min_obs complete observations (after
# Kalman or listwise deletion) are dropped, as in estimate_vars.R.
#
# Arguments:
#   data             - Full data frame: must contain id_col, 'day', variables
#   variables        - Character vector of variable names
#   id_col           - Name of the subject ID column (default "ID")
#   na_method        - "kalman_filter" or "listwise_deletion"
#   train_proportion - Proportion of rows used for training (default 0.80)
#   min_obs          - Minimum complete observations required per person
#
# Returns: named list with
#   $train      - train data (id_col + day + beep [if present] + variables)
#   $test       - test data, same columns
#   $ids_kept   - IDs that passed the min_obs threshold
#   $n_t_train  - named vector of train row counts per person
#   $n_t_test   - named vector of test row counts per person

prepare_mlvar_data <- function(data,
                               variables,
                               id_col           = "ID",
                               na_method        = c("kalman_filter",
                                                    "listwise_deletion"),
                               train_proportion = 0.80,
                               min_obs          = 20) {

  na_method <- match.arg(na_method)

  # --- Identify individuals with enough complete observations ---
  ids_enough <- data %>%
    select(all_of(c(id_col, variables))) %>%
    na.omit() %>%
    count(.data[[id_col]]) %>%
    filter(n > min_obs) %>%
    pull(.data[[id_col]])

  data <- data[data[[id_col]] %in% ids_enough, ]

  # --- Optional Kalman imputation (on real beep sequences, before night rows) ---
  if (na_method == "kalman_filter") {
    message("Applying Kalman filter for missing data imputation...")
    for (id in unique(data[[id_col]])) {
      mask <- data[[id_col]] == id
      data[mask, variables] <- kalman_filter(data[mask, variables])
    }
  }

  ids       <- sort(unique(data[[id_col]]))
  train_list <- list()
  test_list  <- list()
  n_t_train  <- integer(length(ids))
  n_t_test   <- integer(length(ids))
  names(n_t_train) <- as.character(ids)
  names(n_t_test)  <- as.character(ids)

  # Columns to keep: ID, day, beep (if present), variables
  keep_cols <- unique(c(id_col, "day",
                        if ("beep" %in% colnames(data)) "beep",
                        variables))
  keep_cols <- keep_cols[keep_cols %in% colnames(data)]

  for (id in ids) {
    data_i <- data[data[[id_col]] == id, keep_cols, drop = FALSE]

    # Train / test split
    n_rows      <- nrow(data_i)
    split_index <- floor(n_rows * train_proportion)

    # Ensure we have at least 1 test row
    if (split_index >= n_rows) split_index <- n_rows - 1

    train_i <- data_i[1:split_index, ]
    test_i  <- data_i[(split_index + 1):n_rows, ]

    train_list[[as.character(id)]] <- train_i
    test_list[[as.character(id)]]  <- test_i
    n_t_train[as.character(id)]    <- nrow(train_i)
    n_t_test[as.character(id)]     <- nrow(test_i)
  }

  train_data <- do.call(rbind, train_list)
  test_data  <- do.call(rbind, test_list)
  rownames(train_data) <- NULL
  rownames(test_data)  <- NULL

  return(list(
    train      = train_data,
    test       = test_data,
    ids_kept   = ids,
    n_t_train  = n_t_train,
    n_t_test   = n_t_test
  ))
}


# ============================================================
# 3. fit_mlvar()
# ============================================================
# Fits a multilevel VAR(1) model using mlVAR::mlVAR() on the
# training data produced by prepare_mlvar_data().
#
# scale = FALSE because your data is already on [0,1]. Person-mean
# centering is handled internally by mlVAR when scale = FALSE as
# long as the dayvar / beepvar arguments are supplied.
#
# Arguments:
#   prepared_data    - Output of prepare_mlvar_data()
#   variables        - Character vector of variable names
#   id_col           - Name of the ID column (default "ID")
#   contemporaneous  - Random-effects structure for contemporaneous
#                      network: "correlated" (default) or "orthogonal"
#   temporal         - Random-effects structure for temporal network:
#                      "correlated" (default) or "orthogonal"
#                      Use "orthogonal" to reduce computation time for
#                      datasets with many variables
#   lags             - Number of lags (default 1; VAR(1))
#   verbose          - Print mlVAR progress (default FALSE)
#
# Returns: list with
#   $model         - fitted mlVAR object
#   $prepared_data - passed through for use in predict_mlvar()

fit_mlvar <- function(prepared_data,
                      variables,
                      id_col          = "ID",
                      contemporaneous = "correlated",
                      temporal        = "correlated",
                      lags            = 1,
                      verbose         = FALSE) {

  train    <- prepared_data$train
  has_beep <- "beep" %in% colnames(train)
  has_day  <- "day"  %in% colnames(train)

  # --- Check for zero-variance variables in training data ---
  # A variable with no variance causes L-BFGS-B to receive non-finite
  # likelihood values and crash. Drop such variables with a warning.
  vars_sd   <- sapply(variables, function(v) sd(train[[v]], na.rm = TRUE))
  zero_vars <- names(vars_sd[is.na(vars_sd) | vars_sd == 0])
  if (length(zero_vars) > 0) {
    warning(sprintf(
      "Dropping %d zero-variance variable(s) from model: %s",
      length(zero_vars), paste(zero_vars, collapse = ", ")
    ))
    variables <- setdiff(variables, zero_vars)
  }
  if (length(variables) < 2) {
    stop("Fewer than 2 variables remain after dropping zero-variance columns.")
  }

  # --- Fit with L-BFGS-B, fall back to bobyqa on numerical failure ---
  # L-BFGS-B is mlVAR's default optimiser but can fail with non-finite
  # gradients on difficult likelihood surfaces. bobyqa is more robust.
  message("Fitting mlVAR model...")
  t_start <- proc.time()

  fit_args <- list(
    data            = train,
    vars            = variables,
    idvar           = id_col,
    lags            = lags,
    dayvar          = if (has_day)  "day"  else NULL,
    beepvar         = if (has_beep) "beep" else NULL,
    contemporaneous = contemporaneous,
    temporal        = temporal,
    scale           = FALSE,
    verbose         = verbose
  )

  model <- tryCatch(
    do.call(mlVAR::mlVAR, fit_args),
    error = function(e) {
      if (grepl("L-BFGS-B|finite|fn", conditionMessage(e))) {
        message("  L-BFGS-B failed, retrying with bobyqa optimiser...")
        fit_args$optimizer <<- "bobyqa"
        do.call(mlVAR::mlVAR, fit_args)
      } else {
        stop(e)
      }
    }
  )

  elapsed <- round((proc.time() - t_start)["elapsed"], 1)
  message(sprintf("mlVAR fitted in %.1f seconds.", elapsed))

  return(list(
    model         = model,
    prepared_data = prepared_data,
    variables     = variables,
    id_col        = id_col
  ))
}


# ============================================================
# 4. predict_mlvar()
# ============================================================
# Predicts on the test set using subject-level temporal coefficients
# from the fitted mlVAR model, then computes per-variable RMSE.
#
# Prediction approach (mirrors ResAnalysis in estimate_ml_var.R):
#   For each individual j and time t:
#     y_hat[t, i] = intercept[j, i] + phi[j][i, ] %*% y_centred[t-1, ]
#   where y_centred uses person means estimated from the TRAIN set.
#   Predictions are only made when t-1 is in the same day (overnight
#   transitions are identified via the day column, not NA rows).
#
# Note on centering:
#   We person-mean-centre using means from the TRAIN set so that the
#   test predictions are properly out-of-sample. This matches the
#   centering mlVAR applied during estimation.
#
# Arguments:
#   fit_result  - Output of fit_mlvar()
#
# Returns: list with
#   $rmse        - named numeric vector: RMSE per variable (pooled across people)
#   $rmse_person - data frame: per-person per-variable RMSE
#   $predictions - named list of prediction matrices (one per person)

predict_mlvar <- function(fit_result) {

  model      <- fit_result$model
  test_data  <- fit_result$prepared_data$test
  train_data <- fit_result$prepared_data$train
  variables  <- fit_result$variables
  id_col     <- fit_result$id_col
  ids        <- fit_result$prepared_data$ids_kept
  p          <- length(variables)

  # Containers for pooled RMSE computation
  ss_total   <- setNames(numeric(p), variables)  # sum of squared errors
  n_total    <- 0L

  rmse_person_rows <- list()
  all_predictions  <- list()

  for (j in seq_along(ids)) {

    id <- ids[j]
    id_chr <- as.character(id)

    # --- Subject-level temporal coefficients and intercepts ---
    # getNet() returns the p x p temporal coefficient matrix for subject j.
    # Rows = outcome variable, columns = predictor (lag-1) variable.
    phi_j  <- tryCatch(
      mlVAR::getNet(model, type = "temporal", subject = j, verbose = FALSE),
      error = function(e) {
        warning(sprintf("Could not get temporal network for subject %s: %s",
                        id_chr, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(phi_j)) next

    # Subject-level intercepts (mu = person means in mlVAR's parameterisation)
    mu_j <- tryCatch(
      model$results$mu$subject[[j]],
      error = function(e) rep(0, p)
    )
    if (is.null(mu_j) || length(mu_j) != p) mu_j <- rep(0, p)

    # --- Person means from TRAIN set for centering ---
    train_j   <- train_data[train_data[[id_col]] == id, variables, drop = FALSE]
    col_means <- colMeans(train_j, na.rm = TRUE)

    # --- Test data for this person ---
    test_j  <- test_data[test_data[[id_col]] == id, , drop = FALSE]
    N       <- nrow(test_j)
    obs_mat <- as.matrix(test_j[, variables, drop = FALSE])

    # Centre test observations using train means
    obs_centred <- sweep(obs_mat, 2, col_means, "-")

    # --- Identify predictable rows ---
    # A row t is predictable if t-1 is in the same day.
    # We use the day column directly — no night NA rows exist in the data.
    day_vec     <- test_j[["day"]]
    predictable <- rep(FALSE, N)
    for (t in 2:N) {
      if (!is.na(day_vec[t]) && !is.na(day_vec[t - 1]) &&
          day_vec[t] == day_vec[t - 1]) {
        predictable[t] <- TRUE
      }
    }

    # --- Predict ---
    pred_mat <- matrix(NA_real_, nrow = N, ncol = p)
    colnames(pred_mat) <- variables

    for (t in which(predictable)) {
      # VAR(1): y_hat[t] = intercept + phi %*% y_centred[t-1]
      pred_mat[t, ] <- mu_j + phi_j %*% obs_centred[t - 1, ]
    }

    all_predictions[[id_chr]] <- pred_mat

    # --- Per-person RMSE ---
    person_rmse <- setNames(rep(NA_real_, p), variables)
    for (v in variables) {
      obs_v  <- obs_mat[, v]
      pred_v <- pred_mat[, v]
      mask   <- predictable & !is.na(obs_v) & !is.na(pred_v)
      if (sum(mask) > 0) {
        person_rmse[v] <- sqrt(mean((obs_v[mask] - pred_v[mask])^2))
        ss_total[v]    <- ss_total[v] + sum((obs_v[mask] - pred_v[mask])^2)
        n_total        <- n_total + sum(mask)
      }
    }

    rmse_person_rows[[id_chr]] <- c(setNames(id, id_col), person_rmse)
  }

  # --- Pooled RMSE: sqrt( total SS / total N ) per variable ---
  n_vars     <- length(variables)
  rmse_pooled <- sqrt(ss_total / (n_total / n_vars))

  message("RMSE per variable (pooled across individuals):")
  print(round(rmse_pooled, 4))

  # Assemble per-person RMSE data frame
  rmse_person_df <- do.call(rbind, lapply(rmse_person_rows, function(r) {
    as.data.frame(as.list(r), stringsAsFactors = FALSE)
  }))
  rownames(rmse_person_df) <- NULL
  for (v in variables) {
    rmse_person_df[[v]] <- as.numeric(rmse_person_df[[v]])
  }

  return(list(
    rmse        = rmse_pooled,
    rmse_person = rmse_person_df,
    predictions = all_predictions
  ))
}


# ============================================================
# 5. run_mlvar_pipeline()
# ============================================================
# Outer loop over multiple datasets. For each dataset:
#   1. Loads and validates data
#   2. Prepares data (Kalman / listwise + train/test split)
#   3. Fits mlVAR on train data
#   4. Predicts on test data and computes RMSE
#   5. Saves outputs to output_base_dir / study_name /
#
# Arguments:
#   dataset_paths    - Named character vector: names = study names,
#                      values = file paths to CSV files
#   variables_list   - Named list: study name -> character vector of variables
#   output_base_dir  - Base directory for outputs
#   id_col           - Name of the ID column (default "ID"); can be a single
#                      string (applied to all datasets) or a named vector
#   na_method        - "kalman_filter" or "listwise_deletion"
#   train_proportion - Proportion of rows for training (default 0.80)
#   min_obs          - Minimum complete observations per person (default 20)
#   contemporaneous  - mlVAR contemporaneous structure (default "correlated")
#   temporal         - mlVAR temporal structure (default "correlated");
#                      use "orthogonal" to speed up large datasets
#   lags             - Number of lags (default 1)
#
# Returns: named list of results per study, each with:
#   $fit_result    - output of fit_mlvar()
#   $pred_result   - output of predict_mlvar()

run_mlvar_pipeline <- function(dataset_paths,
                               variables_list,
                               output_base_dir,
                               id_col           = "ID",
                               na_method        = c("kalman_filter",
                                                    "listwise_deletion"),
                               train_proportion = 0.80,
                               min_obs          = 20,
                               contemporaneous  = "correlated",
                               temporal         = "correlated",
                               lags             = 1) {

  na_method   <- match.arg(na_method)
  study_names <- names(dataset_paths)

  # Guard against unnamed dataset_paths (would silently skip all studies)
  if (is.null(study_names) || any(study_names == "")) {
    stop("dataset_paths must be a fully named character vector (names = study labels).")
  }

  results <- list()

  for (study in study_names) {

    message("\n", strrep("=", 60))
    message("Processing study: ", study)
    message(strrep("=", 60))

    study_start <- proc.time()

    # Wrap each study in tryCatch so one failure doesn't abort the whole loop
    study_result <- tryCatch({

      # --- Load data ---
      data_study <- read.csv(dataset_paths[[study]])

      # --- Resolve ID column name for this study ---
      # id_col can be a single string or a named vector (one entry per study)
      id_col_study <- if (length(id_col) > 1 && !is.null(names(id_col))) {
        id_col[[study]]
      } else {
        id_col
      }

      # --- Validate required columns ---
      required     <- c(id_col_study, "day", variables_list[[study]])
      missing_cols <- setdiff(required, colnames(data_study))
      if (length(missing_cols) > 0) {
        stop(sprintf("missing columns: %s", paste(missing_cols, collapse = ", ")))
      }

      variables <- variables_list[[study]]

      # --- Prepare data ---
      message("Preparing data (na_method = '", na_method, "')...")
      prepared <- prepare_mlvar_data(
        data             = data_study,
        variables        = variables,
        id_col           = id_col_study,
        na_method        = na_method,
        train_proportion = train_proportion,
        min_obs          = min_obs
      )

      message(sprintf("  %d individuals retained after min_obs filter.",
                      length(prepared$ids_kept)))

      # --- Output directory ---
      study_dir <- file.path(output_base_dir, study)
      dir.create(study_dir, recursive = TRUE, showWarnings = FALSE)

      # --- Fit model ---
      fit_result <- fit_mlvar(
        prepared_data   = prepared,
        variables       = variables,
        id_col          = id_col_study,
        contemporaneous = contemporaneous,
        temporal        = temporal,
        lags            = lags
      )

      saveRDS(fit_result, file.path(study_dir, "mlvar_fit.rds"))

      # --- Predict + RMSE ---
      message("Computing predictions and RMSE on test data...")
      pred_result <- predict_mlvar(fit_result)

      saveRDS(pred_result, file.path(study_dir, "mlvar_pred.rds"))
      write.csv(
        pred_result$rmse_person,
        file.path(study_dir, paste0("RMSE_person_", study, ".csv")),
        row.names = FALSE
      )

      rmse_summary <- data.frame(
        study    = study,
        variable = names(pred_result$rmse),
        rmse     = as.numeric(pred_result$rmse)
      )
      write.csv(
        rmse_summary,
        file.path(study_dir, paste0("RMSE_pooled_", study, ".csv")),
        row.names = FALSE
      )

      # --- Timing ---
      elapsed      <- proc.time() - study_start
      elapsed_mins <- round(elapsed["elapsed"] / 60, 1)
      message(sprintf("Study '%s' complete in %.1f min.", study, elapsed_mins))

      list(fit_result = fit_result, pred_result = pred_result)

    }, error = function(e) {
      message(sprintf("ERROR in study '%s': %s  -- skipping.", study, conditionMessage(e)))
      NULL
    })

    results[[study]] <- study_result
  }

  failed <- sum(sapply(results, is.null))
  if (failed > 0) {
    message(sprintf("\n%d of %d studies failed (stored as NULL in results).",
                    failed, length(study_names)))
  }

  return(results)
}


# ============================================================
# 6. collect_mlvar_rmse()
# ============================================================
# Extracts the pooled per-variable RMSE from run_mlvar_pipeline()
# output into a simple named list that mirrors the structure used
# for mHMM results, so both can be passed directly to compare_rmse.R.
#
# Arguments:
#   mlvar_results - output of run_mlvar_pipeline()
#
# Returns: named list; one named numeric vector per study, e.g.:
#   list(
#     Rowland_2020 = c(Happy = 0.12, Sad = 0.15, ...),
#     Study_B      = c(Happy = 0.11, ...)
#   )

collect_mlvar_rmse <- function(mlvar_results) {
  out <- lapply(mlvar_results, function(study_res) {
    if (is.null(study_res)) return(NULL)
    study_res$pred_result$rmse
  })
  # Drop failed studies (NULL entries)
  out <- Filter(Negate(is.null), out)
  return(out)
}




selected_datasets <- dataset_names[c(3:8, 28)]
dataset_paths_selected <- setNames(
  paste0("Data/CleanRescaled/Rescaled_", selected_datasets, ".csv"),
  selected_datasets
)


dataset_paths <- c(
  Rowland_2020 = "Data/CleanRescaled/Rescaled_Rowland_2020.csv"
)

variables_list <- list(
  Rowland_2020 = c("Happy", "Excited", "Relaxed", "Satisfied",
                   "Angry", "Anxious", "Depressed", "Sad")
)

# --- Run full pipeline ---
mlvar_results <- run_mlvar_pipeline(
  dataset_paths    = dataset_paths_selected,
  variables_list   = variables_list[dataset_names[c(3:8,28)]],
  output_base_dir  = "Output/mlVAR_results/",
  id_col           = "ID",            # change if studies use different column names
  na_method        = "listwise_deletion", # or "listwise_deletion"
  train_proportion = 0.80,
  min_obs          = 20,
  contemporaneous  = "correlated",
  temporal         = "correlated"     # use "orthogonal" if too slow
)




dataset_paths <- c(
  Bernstein_2019  = "Data/CleanRescaled/Rescaled_Bernstein_2019.csv",
  Blanke_2017     = "Data/CleanRescaled/Rescaled_Blanke_2017.csv",
  Brans_2013_ds1  = "Data/CleanRescaled/Rescaled_Brans_2013_ds1.csv",
  Brans_2013_ds2  = "Data/CleanRescaled/Rescaled_Brans_2013_ds2.csv",
  Bringmann_2013  = "Data/CleanRescaled/Rescaled_Bringmann_2013.csv",
  Bringmann_2016  = "Data/CleanRescaled/Rescaled_Bringmann_2016.csv",
  Rowland_2020    = "Data/CleanRescaled/Rescaled_Rowland_2020.csv"
)

variables_list <- list(
  Bernstein_2019 = c("Cheerful", "Anxiety", "Content", "Sad"),
  Blanke_2017    = c("Nervous", "Downhearted", "Distressed", "Happy", "Relaxed", "Content"),
  Brans_2013_ds1 = c("Angry", "Depressed", "Anxiety", "Stressed", "Relaxed", "Happy"),
  Brans_2013_ds2 = c("Angry", "Depressed", "Anxiety", "Sad", "Relaxed", "Happy"),
  Bringmann_2013 = c("Excited", "Worried", "Anxiety", "Sad", "Relaxed"),
  Bringmann_2016 = c("Angry", "Depressed", "Sad", "Anxious", "Relaxed", "Happy"),
  Rowland_2020   = c("Happy", "Excited", "Relaxed", "Satisfied", "Angry", "Anxious", "Depressed", "Sad")
)


mlvar_rmse_list <- collect_mlvar_rmse(mlvar_results)


# ============================================================
# Example usage
# ============================================================

# dataset_paths <- c(
#   Rowland_2020 = "Data/CleanRescaled/Rescaled_Rowland_2020.csv",
#   Study_B      = "Data/CleanRescaled/Rescaled_Study_B.csv"
# )
#
# variables_list <- list(
#   Rowland_2020 = c("Happy", "Excited", "Relaxed", "Satisfied",
#                    "Angry", "Anxious", "Depressed", "Sad"),
#   Study_B      = c("Happy", "Sad", "Anxious")
# )
#
# # --- Run full pipeline ---
# all_results <- run_mlvar_pipeline(
#   dataset_paths    = dataset_paths,
#   variables_list   = variables_list,
#   output_base_dir  = "Output/mlVAR_results/",
#   id_col           = "ID",            # change if studies use different column names
#   na_method        = "kalman_filter", # or "listwise_deletion"
#   train_proportion = 0.80,
#   min_obs          = 20,
#   contemporaneous  = "correlated",
#   temporal         = "correlated"     # use "orthogonal" if too slow
# )
#
# # --- Inspect a specific study ---
# rmse_rowland <- all_results$Rowland_2020$pred_result$rmse
# print(rmse_rowland)
#
# # --- Run individual steps manually (e.g. after checking model) ---
# prepared <- prepare_mlvar_data(
#   data      = read.csv("Data/CleanRescaled/Rescaled_Rowland_2020.csv"),
#   variables = variables_list$Rowland_2020,
#   id_col    = "ID",
#   na_method = "kalman_filter"
# )
#
# fit <- fit_mlvar(prepared, variables_list$Rowland_2020)
# pred <- predict_mlvar(fit)
