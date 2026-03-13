
library(brms)

delta_rmse <- calc_delta_rmse(list("Output/N1_Results/Rescaled_results/VAR_kalman_filter_RMSE/", 
                                   "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/RMSE/"), 
                              dataset_names, 
                              c("var_kf", "hmm_kf"))
studies_info <- read.csv("studiesInfo2.csv") %>%
  rename(dataset = datasets)
data_path <- "Data/CleanRescaled/"
s_state_path <- "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/Nstates/"

get_covariates <- function(data_path, 
                           s_state_path,
                           dataset_names, 
                           studies_info, 
                           delta_rmse){
  
  n_obeservations_list <- list()
  s_states_found_list <- list()
  for(ds in dataset_names){
    n_obeservations_list[[ds]] <- read.csv(paste0(data_path, "rescaled_", ds, ".csv")) %>% 
      group_by(ID) %>%
      summarise(n_observations_per_pp = length(ID)) %>%
      mutate(dataset = ds)
    s_states_found_list[[ds]] <- read.csv(paste0(s_state_path, "Nstates_", ds, ".csv")) %>%
      select(-X) %>%
      mutate(dataset = ds)
  }
  n_obeservations <- do.call(rbind, n_obeservations_list)
  s_states_found <- do.call(rbind, s_states_found_list)
 
  
  delta_rmse %>%
    inner_join(n_obeservations, by = c("ID", "dataset")) %>%
    inner_join(s_states_found, by = c("ID", "dataset")) %>%
    inner_join(studies_info, "dataset") %>%
    select(-X) %>%
    mutate(var_n_paramaters = n_paramaters_var(nvars), 
           hmm_n_paramaters = n_paramaters_hmm(nvars, nstate)) %>%
    rename(scale_min = ScaleMin,
           scale_max = ScaleMax,
           s_state = nstate, 
           n_variables = nvars,
           n_pp = N, 
           n_beebs = nBeeps,
           population = Population) %>%
    mutate(delta_parameters = var_n_paramaters - hmm_n_paramaters, 
           scale_size = scale_max - scale_min) %>%
    return()
}

delta_rmse


delta_rmse %>%
  summarise(var = mean(mean_var_kf),
            hmm = mean(mean_hmm_kf),
            d = mean(delta_rmse))
ds <- "Bar_2020_ds1"

delta_rmse_covariates %>%
  summarise(var = mean(mean_var_kf),
            hmm = mean(mean_hmm_kf),
            d = mean(delta_rmse))


delta_rmse_covariates %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = 100)
  pull(dataset)


delta_rmse_covariates <- get_covariates("Data/CleanRescaled/", 
                                        "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/Nstates/",
                                        dataset_names, 
                                        studies_info,
                                        delta_rmse)
delta_rmse_covariates <- delta_rmse_covariates


delta_rmse_covariates <- get_covariates("Data/CleanRescaled/", 
               "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/Nstates/",
               dataset_names, 
               studies_info,
               delta_rmse) %>%
  mutate(n_observations_per_pp_center = scale(n_observations_per_pp, center = TRUE),
         n_variables_center = scale(n_variables, center = TRUE),
         n_beebs_center = scale(n_beebs, center = TRUE),
         delta_parameters_center = scale(delta_parameters, center = TRUE),
         scale_size_center = scale(scale_size, center = TRUE),
         population = factor(population, levels = c("General", "Clinical", "Students")),
         delta_rmse = delta_rmse*100
         )  %>%
  filter(mean_var_kf < 1)
# %>%
#    select(delta_rmse, n_observations_per_pp_center, n_variables_center, n_beebs_center, 
#         delta_parameters_center, scale_size_center, population, dataset
#          ) %>%
#   as.data.frame() %>%
#   head()


baseline_model <- brms::brm(
  delta_rmse ~ n_observations_per_pp_center + n_variables_center +
    n_beebs_center + delta_parameters_center + scale_size_center + population,
  data = delta_rmse_covariates, 
  family = gaussian(),
  # prior = c(brms::set_prior("normal(0, 5)", class = "b")),  # Weakly informative priors
  chains = 4, 
  cores = 4, 
  iter = 2000
)


random_intercept_model <- brms::brm(
  delta_rmse ~ n_observations_per_pp_center + n_variables_center + n_beebs_center + 
    delta_parameters_center + scale_size_center + population + 
    (1 | dataset),  
  data = delta_rmse_covariates, 
  family = gaussian(),
  # prior = c(brms::set_prior("normal(0, 5)", class = "b"),
  #           brms::set_prior("cauchy(0, 2)", class = "sd")),  
  chains = 4, 
  cores = 4, 
  iter = 2000
)


random_slope_model <- brms::brm(
  delta_rmse ~ n_observations_per_pp_center + n_variables_center + n_beebs_center + 
    delta_parameters_center + scale_size_center + population +
    (1 + n_observations_per_pp_center | dataset),  
  data = delta_rmse_covariates, 
  family = gaussian(),
  # prior = c(brms::set_prior("normal(0, 5)", class = "b"),
  #           brms::set_prior("cauchy(0, 2)", class = "sd")),  
  chains = 4, 
  cores = 4, 
  iter = 2000
)
0.225  -.183
brms::loo_moment_match(brms::add_criterion(random_slope_model, "loo"))

loo <- ?brms::loo(baseline_model, random_intercept_model, random_slope_model)


brms::loo_compare(brms::add_criterion(baseline_model, "loo"), brms::add_criterion(random_intercept_model, "loo"), brms::add_criterion(random_slope_model, "loo"), criterion = "loo")

?brms::bayes_R2(baseline_model)
brms::bayes_R2(random_intercept_model)
brms::bayes_R2(random_slope_model)


