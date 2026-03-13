
path <- "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/Nstates/"


percentage_best_states_list <- list()
for(ds in dataset_names){
  percentage_best_states_list[[ds]] <- read.csv(paste0(path, "Nstates_", ds,".csv")) %>%
    select(-X) %>%
    filter(!is.na(nstate)) %>%
    count(nstate) %>%
    mutate(nstate = gsub(pattern = "(\\d{1})", replacement = "S\\1", nstate)) %>%
    complete(nstate = c("S2", "S3", "S4", "S5"), fill = list(n = 0)) %>% 
    mutate(n = n/sum(n)*100) %>%
    column_to_rownames(var = "nstate") %>%
    t() %>%
    as.data.frame() %>%
    remove_rownames() 
}
  
percentage_best_states <- do.call(rbind, percentage_best_states_list) %>%
  arrange(desc(S2)) %>%
  t()




pdf(file = "Output/Figures/AppendixD_which_state_lowest_aic.pdf", 
    width = standard_figure_with,
    height = 5)


percentage_models_barplot(
  percantages_matrix = percentage_best_states,
  n_datasets = ncol(percentage_best_states), 
  n_models = 5,
  palette =  wesanderson::wes_palette("Zissou1", 5)[1:4],
  margin_bottom_graphs = 1.5,
  margin_dataset_labels = 8,
  dataset_labels_position = -1,
  labels_width = 1.4,
  barplot_width = 1,
  split = 3
)


dev.off()
