# ============================================================
# RMSE Comparison: multilevel HMM vs multilevel VAR
# ============================================================
# Works with two simple named lists of RMSE vectors:
#
#   mhmm_rmse  - built manually from predict_mhmm() outputs:
#     list(
#       Rowland_2020 = predict_mhmm(fit, variables)$rmse,
#       Study_B      = predict_mhmm(fit, variables)$rmse
#     )
#
#   mlvar_rmse - built with collect_mlvar_rmse() from mlvar_pipeline.R:
#     mlvar_rmse <- collect_mlvar_rmse(run_mlvar_pipeline(...))
#
# Both lists have the structure:
#   named list -> named numeric vector (one value per variable)
#
# Functions:
#   1. rmse_to_long()       - Stack both lists into a long data frame
#   2. summarise_rmse()     - Wide table with diff and winner per variable
#   3. plot_rmse_compare()  - Faceted dot plot, one panel per study
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)


# ============================================================
# 1. rmse_to_long()
# ============================================================
# Converts two named RMSE lists into a single long-format data
# frame ready for summarising and plotting.
#
# Arguments:
#   mhmm_rmse  - named list of named numeric vectors (one per study)
#   mlvar_rmse - named list of named numeric vectors (one per study)
#   studies    - character vector of studies to include;
#                defaults to the intersection of both lists
#
# Returns: data frame with columns: study, variable, model, rmse

rmse_to_long <- function(mhmm_rmse, mlvar_rmse, studies = NULL) {

  if (is.null(studies)) {
    studies <- intersect(names(mhmm_rmse), names(mlvar_rmse))
  }
  if (length(studies) == 0) {
    stop("No matching study names found between mhmm_rmse and mlvar_rmse.")
  }

  rows <- list()
  for (study in studies) {
    for (model_name in c("mHMM", "mlVAR")) {
      vec <- if (model_name == "mHMM") mhmm_rmse[[study]] else mlvar_rmse[[study]]
      if (is.null(vec)) {
        warning(sprintf("%s: no RMSE found for study '%s', skipping.", model_name, study))
        next
      }
      rows[[length(rows) + 1]] <- data.frame(
        study    = study,
        variable = names(vec),
        model    = model_name,
        rmse     = as.numeric(vec),
        stringsAsFactors = FALSE
      )
    }
  }

  out       <- do.call(rbind, rows)
  out$model <- factor(out$model, levels = c("mHMM", "mlVAR"))
  rownames(out) <- NULL
  return(out)
}


# ============================================================
# 2. summarise_rmse()
# ============================================================
# Pivots to wide format with mHMM and mlVAR side by side, plus
# the difference (positive = mHMM better, negative = mlVAR better)
# and which model won.
#
# Arguments:
#   rmse_long - output of rmse_to_long()
#
# Returns: data frame with columns:
#   study, variable, rmse_mHMM, rmse_mlVAR, diff, better_model

summarise_rmse <- function(rmse_long) {
  rmse_long %>%
    pivot_wider(names_from = model, values_from = rmse,
                names_prefix = "rmse_") %>%
    mutate(
      diff         = rmse_mlVAR - rmse_mHMM,
      better_model = case_when(
        is.na(diff) ~ NA_character_,
        diff < 0    ~ "mlVAR",
        diff > 0    ~ "mHMM",
        TRUE        ~ "tie"
      )
    ) %>%
    arrange(study, variable)
}


# ============================================================
# 3. plot_rmse_compare()
# ============================================================
# Faceted dot plot: one panel per study, variables on the y-axis,
# RMSE on the x-axis. A line connects the two model estimates per
# variable so differences are easy to read at a glance.
#
# Arguments:
#   rmse_long   - output of rmse_to_long()
#   output_file - optional path to save (PDF or PNG by extension)
#   width       - plot width in inches (default 8)
#   height      - plot height in inches (default: 4 per study)
#
# Returns: ggplot object (invisibly)

plot_rmse_compare <- function(rmse_long,
                              output_file = NULL,
                              width       = 8,
                              height      = NULL) {

  n_studies <- length(unique(rmse_long$study))
  if (is.null(height)) height <- max(4, n_studies * 4)

  p <- ggplot(rmse_long, aes(x = rmse, y = variable,
                              colour = model, shape = model)) +
    geom_line(aes(group = variable), colour = "grey70", linewidth = 0.4) +
    geom_point(size = 2.5) +
    facet_wrap(~ study, scales = "free_x") +
    scale_colour_manual(values = c(mHMM = "#E07B39", mlVAR = "#3A86C8")) +
    scale_shape_manual(values  = c(mHMM = 16, mlVAR = 17)) +
    labs(
      title  = "Prediction RMSE: mHMM vs mlVAR",
      x      = "RMSE (test set, pooled across individuals)",
      y      = NULL,
      colour = "Model", shape = "Model"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank()
    )

  if (!is.null(output_file)) {
    ext <- tolower(tools::file_ext(output_file))
    ggsave(output_file, p, width = width, height = height,
           device = if (ext == "pdf") "pdf" else NULL, dpi = 300)
    message("Plot saved to: ", output_file)
  }

  print(p)
  return(invisible(p))
}


# --- mHMM: build named list manually after predict_mhmm() calls ---
variables_list <- list(
  Bernstein_2019 = c("Cheerful", "Anxiety", "Content", "Sad"),
  Blanke_2017    = c("Nervous", "Downhearted", "Distressed", "Happy", "Relaxed", "Content"),
  Brans_2013_ds1 = c("Angry", "Depressed", "Anxiety", "Stressed", "Relaxed", "Happy"),
  Brans_2013_ds2 = c("Angry", "Depressed", "Anxiety", "Sad", "Relaxed", "Happy"),
  Bringmann_2013 = c("Excited", "Worried", "Anxiety", "Sad", "Relaxed"),
  Bringmann_2016 = c("Angry", "Depressed", "Sad", "Anxious", "Relaxed", "Happy"),
  Rowland_2020   = c("Happy", "Excited", "Relaxed", "Satisfied", "Angry", "Anxious", "Depressed", "Sad")
)


# --- Compare ---
rmse_long    <- rmse_to_long(mlhmm_rmse_list, mlvar_rmse_list)
rmse_summary <- summarise_rmse(rmse_long)

print(rmse_summary)
write.csv(rmse_summary, "Output/rmse_comparison.csv", row.names = FALSE)

plot_rmse_compare(rmse_long, output_file = "Output/rmse_comparison.pdf")

# --- Quick overall win count ---
rmse_summary %>% count(better_model)

# --- Mean RMSE per model across all studies and variables ---
rmse_long %>%
  group_by(model) %>%
  summarise(mean_rmse = mean(rmse, na.rm = TRUE),
            sd_rmse   =   sd(rmse, na.rm = TRUE))




  
list.files("/Documents and Settings/zenog/OneDrive/Documents/Psychologie/ResMas/Y3B4_Thesis/Thesis/", recursive = TRUE)
  

1 + 1
