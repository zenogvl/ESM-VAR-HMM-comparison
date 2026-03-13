# ============================================================
# Multilevel HMM Pipeline using mHMMbayes
# ============================================================
# Structure:
#   1. add_night_blocks()       - Insert NA night gaps between days
#   2. prepare_mhmm_data()      - Format data, train/test split
#   3. make_start_values()      - Auto-generate starting values
#   4. make_priors()            - Auto-generate weakly informative priors
#   5. fit_mhmm()               - Fit models across m states and chains
#   6. check_mhmm()             - Trace plots, Gelman-Rubin, label switching
#   7. predict_mhmm()           - Viterbi decoding + RMSE on test data
#   8. run_mhmm_pipeline()      - Outer wrapper to loop over datasets
# ============================================================

library(mHMMbayes)
library(parallel)
library(RColorBrewer)
library(coda)
library(scales)


# ============================================================
# 1. add_night_blocks()
# ============================================================
# Inserts rows of NAs between days to represent the night gap.
# This is required by mHMMbayes to break the temporal continuity
# between days, so the model does not assume continuity overnight.
#
# Arguments:
#   data                 - Data frame for a single individual,
#                          must contain a column named 'day'
#   n_night_obs          - Number of NA rows to insert per night gap
#
# Returns: data frame with night block NAs inserted between days

add_night_blocks <- function(data, n_night_obs) {
  
  # Create a blank night block with the same columns as the data
  night_block <- as.data.frame(
    matrix(NA, nrow = n_night_obs, ncol = ncol(data))
  )
  colnames(night_block) <- colnames(data)
  
  # Preserve ID so mHMMbayes can identify the subject
  if ("ID" %in% colnames(data)) {
    night_block$ID <- data$ID[1]
  }
  
  # Walk through rows and insert a night block whenever the day changes
  result <- data
  offset  <- 0  # tracks how many rows we have already inserted
  
  for (i in 1:(nrow(data) - 1)) {
    
    current_row <- i + offset
    
    if (data$day[i] != data$day[i + 1]) {
      
      # Split result at the current position and insert the night block
      result <- rbind(
        result[1:current_row, ],
        night_block,
        result[(current_row + 1):nrow(result), ]
      )
      
      offset <- offset + n_night_obs
    }
  }
  
  return(result)
}


# ============================================================
# 2. prepare_mhmm_data()
# ============================================================
# Applies add_night_blocks() per individual, then splits each
# person's data into a train set (first 80%) and test set (last 20%).
# Night block rows are added before the split so the temporal
# structure is preserved in both sets.
#
# Arguments:
#   data             - Full data frame with columns: ID, day, and
#                      all variable columns
#   variables        - Character vector of dependent variable names
#   n_night_obs      - Number of NA rows per night gap
#   train_proportion - Proportion of rows per person used for
#                      training (default 0.80)
#
# Returns: a named list with:
#   $train           - data frame ready for mHMM(), subj_id first
#   $test            - data frame for prediction
#   $n_t_train       - named integer vector: obs per person in train
#   $n_t_test        - named integer vector: obs per person in test

prepare_mhmm_data <- function(data,
                              variables,
                              n_night_obs,
                              train_proportion = 0.80) {
  
  ids <- sort(unique(data$ID))
  
  train_list <- list()
  test_list  <- list()
  n_t_train  <- integer(length(ids))
  n_t_test   <- integer(length(ids))
  names(n_t_train) <- ids
  names(n_t_test)  <- ids
  
  for (id in ids) {
    
    # Subset one individual
    data_i <- data[data$ID == id, ]
    
    # Insert night blocks between days
    data_i_nights <- add_night_blocks(data_i, n_night_obs)
    
    # Train / test split: last 20% of rows are held out
    n_rows      <- nrow(data_i_nights)
    split_index <- floor(n_rows * train_proportion)
    
    train_i <- data_i_nights[1:split_index, ]
    test_i  <- data_i_nights[(split_index + 1):n_rows, ]
    
    # Keep only the columns mHMMbayes needs: ID + dependent variables
    train_i <- train_i[, c("ID", variables)]
    test_i  <- test_i[,  c("ID", variables)]
    
    train_list[[as.character(id)]] <- train_i
    test_list[[as.character(id)]]  <- test_i
    
    n_t_train[as.character(id)] <- nrow(train_i)
    n_t_test[as.character(id)]  <- nrow(test_i)
  }
  
  train_data <- do.call(rbind, train_list)
  test_data  <- do.call(rbind, test_list)
  
  # Rename ID column to subj_id as expected by mHMMbayes
  colnames(train_data)[colnames(train_data) == "ID"] <- "subj_id"
  colnames(test_data)[colnames(test_data)   == "ID"] <- "subj_id"
  
  rownames(train_data) <- NULL
  rownames(test_data)  <- NULL
  
  return(list(
    train    = train_data,
    test     = test_data,
    n_t_train = n_t_train,
    n_t_test  = n_t_test
  ))
}


# ============================================================
# 3. make_start_values()
# ============================================================
# Automatically generates starting values for:
#   - The transition probability matrix (gamma): uniform
#   - The emission distributions: means spread evenly across the
#     observed range of each variable, SDs set to 1/6 of the range
#
# This works for any set of variables and any m, which is important
# when running across 40 different datasets.
#
# Arguments:
#   train_data   - Training data frame (subj_id + variables)
#   variables    - Character vector of dependent variable names
#   m            - Number of hidden states
#
# Returns: list(gamma = matrix, emiss = list of m x 2 matrices)

make_start_values <- function(train_data, variables, m) {
  
  # Transition matrix: equal probability of staying/switching
  start_gamma <- matrix(1 / m, nrow = m, ncol = m)
  
  # Emission starting values: spread means across observed range
  start_emiss <- lapply(variables, function(v) {
    
    obs_vals <- na.omit(train_data[[v]])
    v_min    <- quantile(obs_vals, 0.10)
    v_max    <- quantile(obs_vals, 0.90)
    
    # Evenly spaced means between 10th and 90th percentile
    means <- seq(v_min, v_max, length.out = m)
    
    # SD as a fraction of the observed range
    sds   <- rep((v_max - v_min) / 6, m)
    
    cbind(means, sds)
  })
  
  names(start_emiss) <- variables
  
  return(list(gamma = start_gamma, emiss = start_emiss))
}


# ============================================================
# 4. make_priors()
# ============================================================
# Generates weakly informative priors for continuous emission
# distributions. All prior parameters are exposed as arguments
# so you can override them per dataset if needed.
#
# Arguments:
#   train_data   - Training data frame (subj_id + variables)
#   variables    - Character vector of dependent variable names
#   m            - Number of hidden states
#   emiss_K0     - Prior sample size for emission mean (default 1)
#   emiss_nu     - Degrees of freedom for between-subject variance (default 1)
#   emiss_V_sd   - SD used to set prior between-subject variance (default NULL:
#                  auto = 1/4 of observed range per variable)
#   emiss_a0     - Shape of emission variance prior (default 1.5)
#   emiss_b0     - Scale of emission variance prior (default NULL: auto)
#
# Returns: output of prior_emiss_cont()

make_priors <- function(train_data,
                        variables,
                        m,
                        emiss_K0    = 1,
                        emiss_nu    = 1,
                        emiss_V_sd  = NULL,
                        emiss_a0    = 1.5,
                        emiss_b0    = NULL) {
  
  n_dep <- length(variables)
  
  
  v <-  "Happy" 
  # Compute per-variable prior means (evenly spread over observed range)
  emiss_mu0 <- lapply(variables, function(v) {
    obs_vals <- na.omit(train_data[[v]])
    v_min    <- quantile(obs_vals, 0.10)
    v_max    <- quantile(obs_vals, 0.90)
    matrix(seq(v_min, v_max, length.out = m), nrow = 1)
  })
  names(emiss_mu0) <- variables
  
  # Between-subject variance prior: auto if not supplied
  emiss_V <- lapply(variables, function(v) {
    if (!is.null(emiss_V_sd)) {
      rep(emiss_V_sd^2, m)
    } else {
      obs_vals <- na.omit(train_data[[v]])
      auto_sd  <- (quantile(obs_vals, 0.90) - quantile(obs_vals, 0.10)) / 4
      rep(auto_sd^2, m)
    }
  })
  
  # Emission variance prior scale: auto if not supplied
  b0 <- lapply(variables, function(v) {
    if (!is.null(emiss_b0)) {
      rep(emiss_b0, m)
    } else {
      obs_vals <- na.omit(train_data[[v]])
      auto_b   <- ((quantile(obs_vals, 0.90) - quantile(obs_vals, 0.10)) / 4)
      rep(auto_b, m)
    }
  })
  
  prior <- mHMMbayes::prior_emiss_cont(
    gen        = list(m = m, n_dep = n_dep),
    emiss_mu0  = emiss_mu0,
    emiss_K0   = rep(list(emiss_K0),  n_dep),
    emiss_V    = emiss_V,
    emiss_nu   = rep(list(emiss_nu),  n_dep),
    emiss_a0   =  ,
    emiss_b0   = b0
  )
  
  return(prior)
}


# ============================================================
# 5. fit_mhmm()
# ============================================================
# Fits multilevel HMMs for a range of state numbers (m_range)
# and across multiple chains. Chains are run in parallel using
# parallel::mclapply (unix) or parallel::parLapply (windows).
# Parallelising by chain is preferred over parallelising by dataset
# because it keeps convergence assessment (multi-chain Gelman-Rubin)
# naturally intact.
#
# Arguments:
#   data             - Full data frame (ID, day, variables, ...)
#   variables        - Character vector of dependent variable names
#   m_range          - Integer vector of state numbers to fit (e.g. 2:5)
#   n_chains         - Number of MCMC chains (default 4)
#   n_iter           - Number of Gibbs iterations (default 2000)
#   burn_in          - Burn-in period (default 500)
#   n_night_obs      - Number of NA rows per night gap
#   train_proportion - Proportion of data used for training (default 0.80)
#   output_dir       - Directory to save fitted model .rds files
#   prior_args       - Named list of arguments passed to make_priors()
#                      (use to override defaults per dataset)
#   n_cores          - Number of cores for parallelisation (default: detected - 1)
#
# Returns: list with:
#   $models          - nested list: models[[chain]][[m]]
#   $prepared_data   - output of prepare_mhmm_data()
#   $aic_table       - data frame with AIC per m

fit_mhmm <- function(data,
                     variables,
                     m_range          = 2:4,
                     n_chains         = 4,
                     n_iter           = 2000,
                     burn_in          = 500,
                     n_night_obs      = 8,
                     train_proportion = 0.80,
                     output_dir       = NULL,
                     prior_args       = list(),
                     n_cores          = NULL) {
  
  n_dep <- length(variables)
  
  # --- Prepare data (night blocks + train/test split) ---
  message("Preparing data...")
  prepared   <- prepare_mhmm_data(data, variables, n_night_obs, train_proportion)
  train_data <- prepared$train
  
  # Note: m = 1 is supported by mHMMbayes via a 1x1 transition matrix
  # (the single entry is always 1 — the model stays in state 1 with
  # certainty). No special handling needed; the loop below covers all m.
  fit_result_m1 <- NULL  # kept for AIC table compatibility; unused for m >= 2
  
  # --- Set up parallel backend ---
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  # Cap cores at number of chains — no benefit running more
  n_cores <- min(n_cores, n_chains)
  
  # --- Define the work for one chain ---
  # Each chain fits all models in m_range independently.
  # We pass everything needed as a self-contained list so
  # parallel workers do not need access to the parent environment.
  
  fit_one_chain <- function(chain_id) {
    
    set.seed(chain_id)
    chain_models <- list()
    
    for (m in m_range) {
      
      message(sprintf("  Chain %d | m = %d | starting...", chain_id, m))
      
      # Build starting values and priors fresh for this m
      sv <- make_start_values(train_data, variables, m)
      pr <- do.call(make_priors,
                    c(list(train_data = train_data,
                           variables  = variables,
                           m          = m),
                      prior_args))
      
      # Fit the model
      model <- mHMMbayes::mHMM(
        s_data          = train_data,
        data_distr      = "continuous",
        gen             = list(m = m, n_dep = n_dep),
        start_val       = c(list(sv$gamma), sv$emiss),
        emiss_hyp_prior = pr,
        mcmc            = list(J = n_iter, burn_in = burn_in)
      )
      
      chain_models[[as.character(m)]] <- model
      message(sprintf("  Chain %d | m = %d | done", chain_id, m))
    }
    
    return(chain_models)
  }
  
  # --- Run chains: parallel or sequential ---
  #
  # Parallelisation strategy by platform:
  #
  #   Windows (any):  socket cluster via parLapply works in both RStudio
  #                   and terminal. Each worker is a fresh R session, so
  #                   all functions and packages must be explicitly exported.
  #
  #   Unix/Mac in RStudio: mclapply forks the process, which hangs silently
  #                   inside RStudio. Fall back to sequential in that case.
  #
  #   Unix/Mac in terminal: mclapply works correctly and is most efficient.
  
  os_type    <- .Platform$OS.type
  in_rstudio <- identical(Sys.getenv("RSTUDIO"), "1")
  
  use_parallel <- n_cores > 1
  
  if (use_parallel && os_type == "windows") {
    
    # Windows socket cluster — each worker is a blank R session.
    # We test the cluster before committing to it: if the library load
    # or a smoke-test fails, we fall back to sequential automatically
    # rather than hanging silently.
    
    message(sprintf(
      "Fitting %d chains in parallel on %d cores (Windows socket cluster)...",
      n_chains, n_cores
    ))
    
    cl <- parallel::makeCluster(n_cores)
    
    # --- Smoke test: verify workers can load mHMMbayes ---
    worker_ok <- tryCatch({
      parallel::clusterEvalQ(cl, {
        library(mHMMbayes)
        TRUE
      })
      TRUE
    }, error = function(e) {
      message("  Worker smoke test failed: ", conditionMessage(e))
      FALSE
    })
    
    if (!worker_ok) {
      
      parallel::stopCluster(cl)
      message("  Falling back to sequential — worker setup failed.")
      message("  To debug: run parallel::makeCluster(1) and clusterEvalQ(cl, library(mHMMbayes)) manually.")
      all_chains <- lapply(seq_len(n_chains), fit_one_chain)
      
    } else {
      
      on.exit(parallel::stopCluster(cl), add = TRUE)
      
      # Export all objects and functions that fit_one_chain needs.
      # Socket workers start as blank sessions — nothing from the parent
      # environment is visible unless explicitly exported here.
      parallel::clusterExport(
        cl,
        varlist = c(
          # data and settings
          "train_data", "variables", "m_range", "n_dep",
          "n_iter", "burn_in", "prior_args",
          # helper functions called inside fit_one_chain
          "make_start_values", "make_priors",
          # fit_one_chain itself
          "fit_one_chain"
        ),
        envir = environment()
      )
      
      # Run chains in parallel, catching any per-chain errors
      all_chains <- parallel::parLapply(cl, seq_len(n_chains), function(chain_id) {
        tryCatch(
          fit_one_chain(chain_id),
          error = function(e) {
            list(error = conditionMessage(e), chain_id = chain_id)
          }
        )
      })
      
      # Check for per-chain errors and report them clearly
      for (i in seq_along(all_chains)) {
        if (!is.null(all_chains[[i]]$error)) {
          stop(sprintf(
            "Chain %d failed with error: %s

This usually means a worker could not find a function or package. Check that all helper functions are exported.",
            all_chains[[i]]$chain_id,
            all_chains[[i]]$error
          ))
        }
      }
    }
    
  } else if (use_parallel && !in_rstudio) {
    
    # Unix/Mac running from terminal — mclapply is safe here
    message(sprintf(
      "Fitting %d chains in parallel on %d cores (Unix fork)...",
      n_chains, n_cores
    ))
    
    all_chains <- parallel::mclapply(
      X           = seq_len(n_chains),
      FUN         = fit_one_chain,
      mc.cores    = n_cores,
      mc.set.seed = FALSE
    )
    
  } else {
    
    # Sequential fallback:
    #   - Unix/Mac inside RStudio (mclapply would hang)
    #   - n_cores = 1
    if (in_rstudio && os_type != "windows") {
      message(sprintf(
        "Fitting %d chains sequentially (mclapply disabled inside RStudio on Unix)...",
        n_chains
      ))
      message("Tip: run via Rscript in a terminal for parallel speedup.")
    } else {
      message(sprintf("Fitting %d chains sequentially (n_cores = 1)...", n_chains))
    }
    
    all_chains <- lapply(seq_len(n_chains), fit_one_chain)
  }
  
  names(all_chains) <- paste0("chain_", seq_len(n_chains))
  
  # --- Save per-chain results ---
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (chain_id in seq_len(n_chains)) {
      saveRDS(
        all_chains[[chain_id]],
        file = file.path(output_dir, sprintf("chain_%d.rds", chain_id))
      )
    }
    message("Chain results saved to: ", output_dir)
  }
  
  # --- Compute AIC table ---
  # AIC = -2 * log-likelihood + 2 * n_parameters
  #
  # mHMMbayes does not store the log-likelihood in the model object.
  # We therefore compute it from the posterior mean parameters:
  #   - Emission means and SDs are taken from emiss_mu_bar and emiss_sd_bar
  #     (posterior means across all sampled iterations)
  #   - Transition probabilities are taken from gamma_prob_bar
  #   - The most likely state sequence is decoded via the Viterbi algorithm
  #   - The LL is then sum of log p(obs | state, emission params) over all
  #     non-missing observations, weighted by the decoded state sequence
  #
  # Free parameters per model:
  #   gamma:  m * (m-1)    [0 for m=1, since the only value is fixed at 1]
  #   emiss:  m * n_dep * 2  [mean + SD per state per variable]
  
  # --- Helper: compute LL and AIC directly from model object ---
  # Based on the mHMMbayes print.mHMM source code, the per-subject LL is
  # stored in x$PD_subj[[s]]$log_likl as a [J x 1] matrix. The print method
  # takes the median over post-burn-in iterations for each subject.
  # We replicate this exactly, then sum over subjects for the total LL.
  #
  # AIC formula from mHMMbayes source (continuous case):
  #   n_par = m * n_dep * 2 + (m-1) * m
  #   AIC   = 2 * n_par - 2 * LL   (computed per subject, then averaged)
  #
  # We follow the same per-subject approach to stay consistent with what
  # mHMMbayes itself reports, storing both the per-subject average and the
  # total (summed) LL so both are available for downstream use.
  
  compute_aic_from_model <- function(model) {
    
    input   <- model$input
    n_subj  <- input$n_subj
    burn_in <- input$burn_in
    J       <- input$J
    m_val   <- input$m
    n_dep_m <- input$n_dep
    n_vary  <- input$n_vary
    
    # Per-subject LL: median over post-burn-in Gibbs samples
    ll_per_subj <- numeric(n_subj)
    for (s in seq_len(n_subj)) {
      ll_per_subj[s] <- median(
        model$PD_subj[[s]]$log_likl[(burn_in + 1):J, 1]
      )
    }
    
    # AIC per subject (continuous emission formula from mHMMbayes source)
    n_par     <- m_val * n_dep_m * 2 + (m_val - 1) * m_val
    aic_per_subj  <- 2 * n_par - 2 * ll_per_subj
    aicc_per_subj <- ((2 * n_vary * n_par) / (n_vary - n_par - 1)) - 2 * ll_per_subj
    
    list(
      ll_mean   = mean(ll_per_subj),   # average LL per subject (matches print output)
      ll_total  = sum(ll_per_subj),    # total LL across all subjects
      aic_mean  = mean(aic_per_subj),  # average AIC per subject (matches print output)
      aicc_mean = mean(aicc_per_subj), # average AICc per subject
      n_par     = n_par
    )
  }
  
  aic_rows <- list()
  
  for (m_val in m_range) {
    
    model <- all_chains[[1]][[as.character(m_val)]]
    
    fit_stats <- tryCatch(
      compute_aic_from_model(model),
      error = function(e) {
        warning(sprintf("AIC computation failed for m = %d: %s", m_val,
                        conditionMessage(e)))
        NULL
      }
    )
    
    if (is.null(fit_stats)) {
      aic_rows[[as.character(m_val)]] <- data.frame(
        m         = m_val,
        ll_mean   = NA_real_,
        ll_total  = NA_real_,
        n_par     = NA_integer_,
        aic_mean  = NA_real_,
        aicc_mean = NA_real_
      )
    } else {
      aic_rows[[as.character(m_val)]] <- data.frame(
        m         = m_val,
        ll_mean   = fit_stats$ll_mean,
        ll_total  = fit_stats$ll_total,
        n_par     = fit_stats$n_par,
        aic_mean  = fit_stats$aic_mean,
        aicc_mean = fit_stats$aicc_mean
      )
    }
  }
  
  aic_table           <- do.call(rbind, aic_rows)
  rownames(aic_table) <- NULL
  
  # Select best m using average AIC per subject (consistent with mHMMbayes output)
  best_m_aic <- aic_table$m[which.min(aic_table$aic_mean)]
  message("\nAIC table:")
  print(aic_table)
  message("Best m by AIC: ", best_m_aic)
  
  return(list(
    models        = all_chains,
    prepared_data = prepared,
    aic_table     = aic_table
  ))
}


# ============================================================
# 6. check_mhmm()
# ============================================================
# Produces diagnostic plots and statistics to assess:
#   (a) Convergence   — trace plots + Gelman-Rubin Rhat per m
#   (b) Label switching — per-subject emission trace plots
#
# This function is intentionally separate from fit_mhmm() so
# you can inspect diagnostics before proceeding to prediction.
#
# Arguments:
#   fit_result    - Output of fit_mhmm()
#   variables     - Character vector of dependent variable names
#   m_range       - Integer vector of state numbers that were fitted
#   output_dir    - Directory to save PDFs (NULL = no saving)
#   n_chains      - Number of chains in fit_result
#
# Returns: list with $gr_stats (Gelman-Rubin per m) printed to console

check_mhmm <- function(fit_result,
                       variables,
                       m_range   = NULL,
                       output_dir = NULL,
                       n_chains  = NULL) {
  
  all_chains <- fit_result$models
  
  # Infer m_range and n_chains from the data if not supplied
  if (is.null(m_range)) {
    m_range  <- as.integer(names(all_chains[[1]]))
  }
  if (is.null(n_chains)) {
    n_chains <- length(all_chains)
  }
  
  n_dep <- length(variables)
  cols  <- RColorBrewer::brewer.pal(max(4, n_chains), "Set1")[1:n_chains]
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ---- (a) Trace Plots & Gelman-Rubin ----
  gr_stats <- list()
  
  for (m in m_range) {
    
    chain_models <- lapply(all_chains, function(ch) ch[[as.character(m)]])
    
    # --- Trace plot ---
    # For m = 1: no transition matrix exists, so only emission traces are plotted.
    # For m >= 2: one extra panel at the top shows the logit TPM trace.
    plot_file <- if (!is.null(output_dir)) {
      file.path(output_dir, sprintf("traceplots_m%d.pdf", m))
    } else NULL
    
    n_panels <- if (m == 1) n_dep else n_dep + 1
    if (!is.null(plot_file)) pdf(plot_file, width = 10, height = 3 * n_panels)
    par(mfrow = c(n_panels, 1), mar = c(3, 4, 2, 1))
    
    # Transition probability trace — only for m >= 2
    if (m >= 2) {
      plot(
        chain_models[[1]]$gamma_int_bar[, 1],
        type  = "l", col = cols[1],
        main  = sprintf("m = %d | Logit TPM [1,2]", m),
        xlab  = "Iteration", ylab  = "Logit gamma[1,2]",
        ylim  = range(sapply(chain_models,
                             function(ch) range(ch$gamma_int_bar[, 1])))
      )
      for (chain_id in 2:n_chains) {
        lines(chain_models[[chain_id]]$gamma_int_bar[, 1], col = cols[chain_id])
      }
      legend("topright", legend = paste0("Chain ", 1:n_chains),
             col = cols, lty = 1, bty = "n", cex = 0.8)
    }
    
    # Emission mean traces for all m.
    # For m >= 2 we show the highest state (state m) as it is most informative
    # for label switching. For m = 1 there is only state 1.
    state_to_plot <- m
    for (v in variables) {
      all_vals <- sapply(chain_models,
                         function(ch) ch$emiss_mu_bar[[v]][, state_to_plot])
      plot(
        chain_models[[1]]$emiss_mu_bar[[v]][, state_to_plot],
        type  = "l", col = cols[1],
        main  = sprintf("Emission mean: %s | State %d", v, state_to_plot),
        xlab  = "Iteration", ylab  = "Mean",
        ylim  = range(all_vals)
      )
      for (chain_id in 2:n_chains) {
        lines(chain_models[[chain_id]]$emiss_mu_bar[[v]][, state_to_plot],
              col = cols[chain_id])
      }
      # Add legend on first emission panel for m = 1 (no gamma panel above)
      if (m == 1 && v == variables[1]) {
        legend("topright", legend = paste0("Chain ", 1:n_chains),
               col = cols, lty = 1, bty = "n", cex = 0.8)
      }
    }
    
    if (!is.null(plot_file)) {
      dev.off()
      message("Saved: ", plot_file)
    }
    
    # --- Gelman-Rubin Rhat ---
    # For m = 1 there are no transition parameters, so we only compute Rhat
    # for emission means. For m >= 2 we also compute Rhat for the TPM.
    
    # Transition matrix Rhat — only for m >= 2
    gr_gamma <- if (m >= 2) {
      tryCatch({
        gamma_mcmc <- lapply(chain_models, function(ch) {
          coda::mcmc(ch$gamma_int_bar)
        })
        coda::gelman.diag(coda::mcmc.list(gamma_mcmc))$psrf
      }, error = function(e) NULL)
    } else {
      NULL
    }
    
    gr_emiss <- tryCatch({
      emiss_mcmc <- lapply(chain_models, function(ch) {
        # Bind all emission mean traces column-wise
        do.call(cbind, lapply(variables, function(v) ch$emiss_mu_bar[[v]]))
      })
      emiss_mcmc_list <- coda::mcmc.list(lapply(emiss_mcmc, coda::mcmc))
      coda::gelman.diag(emiss_mcmc_list)$psrf
    }, error = function(e) NULL)
    
    gr_stats[[as.character(m)]] <- list(
      gamma  = gr_gamma,
      emiss  = gr_emiss
    )
    
    # Print summary to console
    message(sprintf("\n--- Gelman-Rubin Rhat | m = %d ---", m))
    if (!is.null(gr_gamma)) {
      message("  Transition matrix (mean Rhat): ",
              round(mean(gr_gamma[, 1], na.rm = TRUE), 3))
    }
    if (!is.null(gr_emiss)) {
      message("  Emission means    (mean Rhat): ",
              round(mean(gr_emiss[, 1], na.rm = TRUE), 3))
    }
    message("  (Rhat < 1.1 indicates good convergence)")
  }
  
  # ---- (b) Label Switching Plots ----
  # For each m, plot the per-subject emission traces for each state.
  # If chains agree on which state is which, there is no label switching.
  
  for (m in m_range) {
    
    # Label switching requires >= 2 states — with only one state there are
    # no labels that can switch, so we skip and note it in the console.
    if (m < 2) {
      message("  m = 1: label switching not applicable (only one state).")
      next
    }
    
    chain_models <- lapply(all_chains, function(ch) ch[[as.character(m)]])
    state_cols   <- RColorBrewer::brewer.pal(max(3, m), "Dark2")[1:m]
    
    # Layout: n_dep rows (one per variable) x subjects_per_page columns.
    # Each page shows subjects_per_page subjects side by side, so the total
    # number of pages = ceiling(n_subjects / subjects_per_page).
    subjects_per_page <- 5
    
    plot_file <- if (!is.null(output_dir)) {
      file.path(output_dir, sprintf("label_switching_m%d.pdf", m))
    } else NULL
    
    if (!is.null(plot_file)) {
      # Width: one column per subject on the page; height: one row per variable
      pdf(plot_file,
          width  = subjects_per_page * 3,
          height = n_dep * 2.5)
      message("Saving label switching plot: ", plot_file)
    }
    
    # Use chain 1 only — label switching shows up as state traces crossing
    model_m    <- chain_models[[1]]
    n_subjects <- length(model_m$PD_subj)
    n_pages    <- ceiling(n_subjects / subjects_per_page)
    
    for (page in seq_len(n_pages)) {
      
      # Indices of subjects on this page
      subj_start <- (page - 1) * subjects_per_page + 1
      subj_end   <- min(page * subjects_per_page, n_subjects)
      subj_page  <- subj_start:subj_end
      n_on_page  <- length(subj_page)
      
      # Grid: n_dep rows x n_on_page columns
      par(mfrow = c(n_dep, n_on_page), mar = c(2, 3, 2, 1), oma = c(0, 0, 2, 0))
      
      for (v_idx in seq_along(variables)) {
        v <- variables[v_idx]
        
        # Columns in PD_subj for this variable: v_idx*m - m + 1 : v_idx*m
        col_start <- v_idx * m - m + 1
        col_end   <- v_idx * m
        
        for (subj_idx in subj_page) {
          
          subj_traces <- model_m$PD_subj[[subj_idx]]$cont_emiss[, col_start:col_end, drop = FALSE]
          
          matplot(
            subj_traces,
            type  = "l", lty = 1,
            col   = state_cols,
            main  = sprintf("S%d | %s", subj_idx, v),
            xlab  = "", ylab = if (subj_idx == subj_start) v else "",
            ylim  = c(0, 1),
            cex.main = 0.8, cex.axis = 0.7
          )
        }
      }
      
      # Page title
      title(
        main = sprintf("Label switching | m = %d | Page %d of %d", m, page, n_pages),
        outer = TRUE, cex.main = 1
      )
    }
    
    if (!is.null(plot_file)) dev.off()
  }
  
  # --- Save Gelman-Rubin results ---
  # Saves two formats:
  #   gr_stats.rds     — full psrf objects, useful for detailed inspection
  #   gr_summary.csv   — tidy table of mean Rhat per m, easy to report
  
  if (!is.null(output_dir)) {
    
    # Save full GR objects as RDS
    saveRDS(gr_stats, file.path(output_dir, "gr_stats.rds"))
    
    # Build a tidy summary table: one row per m, mean Rhat for gamma and emiss
    gr_rows <- lapply(names(gr_stats), function(m_chr) {
      gs <- gr_stats[[m_chr]]
      data.frame(
        m             = as.integer(m_chr),
        gamma_rhat    = if (!is.null(gs$gamma))
          round(mean(gs$gamma[, 1],  na.rm = TRUE), 3)
        else NA_real_,
        emiss_rhat    = if (!is.null(gs$emiss))
          round(mean(gs$emiss[, 1],  na.rm = TRUE), 3)
        else NA_real_,
        gamma_rhat_max = if (!is.null(gs$gamma))
          round(max(gs$gamma[, 1],   na.rm = TRUE), 3)
        else NA_real_,
        emiss_rhat_max = if (!is.null(gs$emiss))
          round(max(gs$emiss[, 1],   na.rm = TRUE), 3)
        else NA_real_
      )
    })
    
    gr_summary <- do.call(rbind, gr_rows)
    write.csv(gr_summary,
              file.path(output_dir, "gr_summary.csv"),
              row.names = FALSE)
    
    message("Gelman-Rubin results saved to: ", output_dir)
    message("\nGelman-Rubin summary (Rhat < 1.1 = good convergence):")
    print(gr_summary)
  }
  
  return(invisible(gr_stats))
}


# ============================================================
# 7. predict_mhmm()
# ============================================================
# Uses the Viterbi algorithm to decode the most likely state
# sequence on the test data, then predicts emissions using the
# group-level emission means, and computes RMSE per variable.
#
# Note: mHMMbayes does not have a built-in predict() function,
# so we use the posterior group-level emission means (after
# burn-in) as fixed parameters for prediction.
#
# Arguments:
#   fit_result    - Output of fit_mhmm()
#   variables     - Character vector of dependent variable names
#   best_m        - Number of states to use for prediction.
#                   If NULL, selected automatically by AIC.
#   chain_id      - Which chain to use for prediction (default 1)
#
# Returns: list with $rmse, $states, $predictions, $best_m

predict_mhmm <- function(fit_result,
                         variables,
                         best_m   = NULL,
                         chain_id = 1) {
  
  aic_table  <- fit_result$aic_table
  test_data  <- fit_result$prepared_data$test
  all_chains <- fit_result$models
  
  # --- Select best m ---
  if (is.null(best_m)) {
    best_m <- aic_table$m[which.min(aic_table$AIC)]
    message("Auto-selected m = ", best_m, " by AIC")
  } else {
    message("Using user-specified m = ", best_m)
  }
  
  model <- all_chains[[chain_id]][[as.character(best_m)]]
  
  if (is.null(model)) {
    stop(sprintf("No fitted model found for m = %d in chain %d", best_m, chain_id))
  }
  
  # --- Extract group-level posterior means (after burn-in) ---
  # emiss_mu_bar contains the full chain of group-level means [J x m]
  # We take the posterior mean (columns = states, rows = iterations)
  
  emiss_means <- sapply(variables, function(v) {
    colMeans(model$emiss_mu_bar[[v]])  # length m vector: mean per state
  })
  # emiss_means is now [m x n_dep]
  
  # --- Extract group-level transition matrix ---
  # gamma_prob_bar contains [J x m*m] posterior samples
  # We take the posterior mean and reshape
  gamma_mean <- matrix(
    colMeans(model$gamma_prob_bar),
    nrow = best_m, ncol = best_m, byrow = TRUE
  )
  
  # --- Viterbi decoding per subject on test data ---
  ids <- unique(test_data$subj_id)
  
  all_states      <- list()
  all_predictions <- list()
  rmse_per_var    <- setNames(numeric(length(variables)), variables)
  n_pred_total    <- 0
  
  for (id in ids) {
    
    test_i <- test_data[test_data$subj_id == id, variables, drop = FALSE]
    
    # Skip subjects with too many NAs to decode
    if (mean(is.na(test_i)) > 0.9) next
    
    # Run Viterbi: find most likely state sequence
    # We implement a simple forward-Viterbi here since mHMMbayes
    # does not expose a standalone predict function
    n_obs  <- nrow(test_i)
    states <- rep(NA_integer_, n_obs)
    
    # Initialise: equal initial state probabilities
    log_delta <- log(rep(1 / best_m, best_m))
    psi       <- matrix(0L, n_obs, best_m)
    
    for (t in seq_len(n_obs)) {
      
      obs_t <- as.numeric(test_i[t, ])
      
      if (all(is.na(obs_t))) {
        # NA row (night gap): carry forward without update
        psi[t, ] <- seq_len(best_m)
        next
      }
      
      # Log emission probability under Gaussian emission for each state
      log_emit <- sapply(seq_len(best_m), function(s) {
        sum(sapply(seq_along(variables), function(v_idx) {
          x_v <- obs_t[v_idx]
          if (is.na(x_v)) return(0)  # ignore missing variables
          dnorm(x_v, mean = emiss_means[s, v_idx],
                sd   = sqrt(1 / best_m),  # placeholder SD
                log  = TRUE)
        }))
      })
      
      if (t == 1) {
        log_delta <- log_delta + log_emit
        psi[t, ]  <- seq_len(best_m)
      } else {
        # Transition: log(gamma) + previous delta
        trans_mat <- log(gamma_mean) + matrix(log_delta, best_m, best_m, byrow = FALSE)
        log_delta <- apply(trans_mat, 2, max) + log_emit
        psi[t, ]  <- apply(trans_mat, 2, which.max)
      }
    }
    
    # Backtrack
    decoded        <- integer(n_obs)
    decoded[n_obs] <- which.max(log_delta)
    for (t in (n_obs - 1):1) {
      decoded[t] <- psi[t + 1, decoded[t + 1]]
    }
    
    all_states[[as.character(id)]] <- decoded
    
    # --- Predict emission values using state means ---
    pred_i <- matrix(NA_real_, n_obs, length(variables))
    colnames(pred_i) <- variables
    for (t in seq_len(n_obs)) {
      if (!is.na(decoded[t])) {
        pred_i[t, ] <- emiss_means[decoded[t], ]
      }
    }
    all_predictions[[as.character(id)]] <- pred_i
    
    # --- Accumulate squared errors ---
    for (v in variables) {
      obs_v    <- test_i[[v]]
      pred_v   <- pred_i[, v]
      mask     <- !is.na(obs_v) & !is.na(pred_v)
      rmse_per_var[v] <- rmse_per_var[v] + sum((obs_v[mask] - pred_v[mask])^2)
      n_pred_total    <- n_pred_total + sum(mask)
    }
  }
  
  # Finalise RMSE
  rmse_per_var <- sqrt(rmse_per_var / (n_pred_total / length(variables)))
  
  message("RMSE per variable:")
  print(round(rmse_per_var, 4))
  
  return(list(
    best_m      = best_m,
    rmse        = rmse_per_var,
    states      = all_states,
    predictions = all_predictions,
    aic_table   = aic_table
  ))
}


# ============================================================
# 8. run_mhmm_pipeline()
# ============================================================
# Outer loop over multiple datasets. For each dataset:
#   1. Fits models via fit_mhmm()
#   2. Saves diagnostics via check_mhmm()
#   3. Stores fitted objects and AIC tables
#
# Prediction (predict_mhmm) is intentionally left out of this
# loop so you can inspect convergence and label switching first,
# then call predict_mhmm() manually per dataset.
#
# Arguments:
#   dataset_paths    - Named character vector: names = study names,
#                      values = file paths to CSV files
#   variables_list   - Named list of character vectors: study name ->
#                      variable names for that study
#   output_base_dir  - Base directory for all outputs
#   m_range          - Integer vector of states to fit (e.g. 2:4)
#   n_chains         - Number of MCMC chains
#   n_iter           - Number of Gibbs iterations
#   burn_in          - Burn-in
#   n_night_obs      - Number of NA rows per night gap
#   train_proportion - Training proportion (default 0.80)
#   prior_args_list  - Optional named list of prior_args per study;
#                      if NULL, auto priors are used for all datasets
#   n_cores          - Cores for parallelisation
#
# Returns: named list of fit_mhmm() outputs, one per study

run_mhmm_pipeline <- function(dataset_paths,
                              variables_list,
                              output_base_dir,
                              m_range          = 2:4,
                              n_chains         = 4,
                              n_iter           = 2000,
                              burn_in          = 500,
                              n_night_obs      = 8,
                              train_proportion = 0.80,
                              prior_args_list  = NULL,
                              n_cores          = NULL) {
  
  study_names <- names(dataset_paths)
  results     <- list()
  
  for (study in study_names) {
    
    message("\n", strrep("=", 60))
    message("Processing study: ", study)
    message(strrep("=", 60))
    
    # Record start time so we can report how long estimation took
    study_start_time <- proc.time()
    
    # Load data
    data_study <- read.csv(dataset_paths[[study]])
    
    # Ensure ID and day columns exist (adapt column names as needed)
    if (!"ID" %in% colnames(data_study)) {
      stop(sprintf("Study '%s': no 'ID' column found.", study))
    }
    if (!"day" %in% colnames(data_study)) {
      stop(sprintf("Study '%s': no 'day' column found.", study))
    }
    
    variables  <- variables_list[[study]]
    prior_args <- if (!is.null(prior_args_list)) prior_args_list[[study]] else list()
    
    # Output directories
    fit_dir   <- file.path(output_base_dir, study, "chains")
    check_dir <- file.path(output_base_dir, study, "diagnostics")
    
    # Fit models
    fit_result <- fit_mhmm(
      data             = data_study,
      variables        = variables,
      m_range          = m_range,
      n_chains         = n_chains,
      n_iter           = n_iter,
      burn_in          = burn_in,
      n_night_obs      = n_night_obs,
      train_proportion = train_proportion,
      output_dir       = fit_dir,
      prior_args       = prior_args,
      n_cores          = n_cores
    )
    
    # Save diagnostics
    check_mhmm(
      fit_result = fit_result,
      variables  = variables,
      m_range    = m_range,
      output_dir = check_dir,
      n_chains   = n_chains
    )
    
    # Store and save
    results[[study]] <- fit_result
    saveRDS(fit_result, file.path(output_base_dir, study, "fit_result.rds"))
    
    # Report elapsed time for this study
    elapsed       <- proc.time() - study_start_time
    elapsed_mins  <- round(elapsed["elapsed"] / 60, 1)
    elapsed_secs  <- round(elapsed["elapsed"], 0)
    message(sprintf(
      "Study '%s' complete. Time elapsed: %s min (%s sec).",
      study, elapsed_mins, elapsed_secs
    ))
    message("AIC table:")
    print(fit_result$aic_table)
  }
  
  return(invisible(results))
}

selected_datasets <- dataset_names[c(3:8, 28)]
dataset_paths_selected <- setNames(
  paste0("Data/CleanRescaled/Rescaled_", selected_datasets, ".csv"),
  selected_datasets
)

# Run the full pipeline (fit + diagnostics)
all_results <- run_mhmm_pipeline(
  dataset_paths   = dataset_paths_selected,
  variables_list  = variables_list[dataset_names[c(3:8,28)]],
  output_base_dir = "Output/mHMM_results/",
  m_range         = 1:4,
  n_chains        = 4,
  n_iter          = 2000,
  burn_in         = 500,
  n_night_obs     = 8
)

all_results <- run_mhmm_pipeline(
  dataset_paths   = dataset_paths_selected[7],
  variables_list  = variables_list[dataset_names[c(3:8,28)]][7],
  output_base_dir = "Output/mHMM_results/",
  m_range         = 1:3,
  n_chains        = 4,
  n_iter          = 2000,
  burn_in         = 500,
  n_night_obs     = 8
)





str(all_results$Rowland_test$models$chain_1[["2"]], max.level = 1)


all_results <- run_mhmm_pipeline(
  dataset_paths   = c(Rowland_test = "Data/CleanRescaled/Rescaled_Rowland_2020.csv"),
  variables_list  = list(Rowland_test = variables_list$Rowland_2020) ,
  output_base_dir = "Output/mHMM_results/",
  m_range         = 1:2,
  n_chains        = 4,
  n_iter          = 40,
  burn_in         = 10,
  n_night_obs     = 8
)


# Load all fit results from disk
selected_datasets <- dataset_names[c(3:8, 28)]  # adjust indices to match your run

mlhmm_results <- setNames(
  lapply(selected_datasets, function(study) {
    readRDS(file.path("Output/mHMM_results", study, "fit_result.rds"))
  }),
  selected_datasets
)

# Quick check — should print the AIC table for each study
for (study in names(mlhmm_results)) {
  message("Loaded: ", study)
  print(mlhmm_results[[study]]$aic_table)
}

best_m_list <- list(
  Bernstein_2019  = 2,
  Blanke_2017     = 2,
  Brans_2013_ds1  = 3,
  Brans_2013_ds2  = 3,
  Bringmann_2013  = 4,
  Bringmann_2016  = 3,
  Rowland_2020    = 3
)

# Loop over datasets and predict using the manually selected m
for (study in names(best_m_list)) {
  
  message("Predicting study: ", study)
  
  pred_result <- predict_mhmm(
    fit_result = mlhmm_results[[study]],
    variables  = variables_list[[study]],
    best_m     = best_m_list[[study]]
  )
  
  # Save prediction results per study
  saveRDS(pred_result, file.path("Output/mHMM_results", study, "pred_result.rds"))
  
  message("  RMSE for ", study, ":")
  print(pred_result$rmse)
}


mlhmm_rmse_list <- setNames(
  lapply(selected_datasets, function(study) {
    readRDS(file.path("Output/mHMM_results", study, "pred_result.rds"))$rmse
  }),
  selected_datasets
)

# Print each study's RMSE for easy inspection
for (study in names(mlhmm_rmse_list)) {
  message(study, ":")
  print(round(mlhmm_rmse_list[[study]], 4))
}





# ============================================================
# Example usage
# ============================================================

# # Define datasets and their variables
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
# # Run the full pipeline (fit + diagnostics)
# all_results <- run_mhmm_pipeline(
#   dataset_paths   = dataset_paths,
#   variables_list  = variables_list,
#   output_base_dir = "Output/mHMM_results/",
#   m_range         = 2:4,
#   n_chains        = 4,
#   n_iter          = 2000,
#   burn_in         = 500,
#   n_night_obs     = 8
# )
#
# # --- After inspecting diagnostics, run prediction ---
# # Automatic m selection by AIC:
# pred_rowland <- predict_mhmm(all_results$Rowland_2020,
#                              variables = variables_list$Rowland_2020)
#
# # Manual m override after inspection:
# pred_rowland <- predict_mhmm(all_results$Rowland_2020,
#                              variables = variables_list$Rowland_2020,
#                              best_m    = 3)
