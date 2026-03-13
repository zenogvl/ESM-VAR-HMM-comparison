

hmm_path <- "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/Ntries/"
var_path <- "Output/N1_Results/Rescaled_results/VAR_kalman_filter_RMSE/"


percentage_valid_model_list <- list()
for(ds in dataset_names_all){
  rmse_var <- read.csv(paste0(var_path, "RMSE_", ds,".csv")) %>%
    select(-X) %>% 
    mutate(var_mean_rmse = rowMeans(across(-ID), na.rm = TRUE)) %>%
    mutate(var_mean_rmse = ifelse(var_mean_rmse > 1, NA, var_mean_rmse)) %>%
    select(ID, var_mean_rmse) 
  
  ntries_hmm <- read.csv(paste0(hmm_path, "Ntries_", ds,".csv")) %>%
    select(-X) %>%
    mutate(across(starts_with("X"), ~ replace(.x, .x == 100, NA)))
  
  percentage_valid_model_list[[ds]] <- full_join(rmse_var, ntries_hmm, by = "ID") %>%
    mutate(across(-ID, ~replace(.x, !is.na(.x), 1))) %>%
    mutate(across(-ID, ~replace(.x, is.na(.x), 0))) %>%
    summarise(across(-ID, sum)/ length(ID) * 100) %>% 
    mutate(dataset = ds) %>%
    rename(S2 = X2, S3 = X3, S4 = X4, S5 = X5, var = var_mean_rmse) %>%
    select(dataset, S2, S3, S4, S5, var)
}
percentage_valid_models <- do.call(rbind, percentage_valid_model_list) %>%
  arrange(desc(S2),desc(S3)) %>%
  select(-dataset) %>%
  t()


pdf(file = "Output/Figures/AppendixC_percentages_valid_models.pdf", 
    width = standard_figure_with,
    height = 4)


percentage_models_barplot(
  percantages_matrix = percentage_valid_models,
  n_datasets =  ncol(percentage_valid_models), 
  n_models = 6,
  palette = wesanderson::wes_palette("Zissou1", 5), 
  margin_bottom_graphs = 1.5,
  margin_dataset_labels = 8,
  dataset_labels_position = -1,
  labels_width = 1.4,
  barplot_width = 1,
  split = 3
)


dev.off()













layout(matrix(c(1:(split*2),rep(split*2+1, split*2)), nrow = 2, byrow = TRUE), 
       widths = rep(c(labels_width, barplot_width), split),
       heights = c(12,1))
par(oma = c(1, 1, 1, 2)) 

for(i in 1:split){

  input_selection <- percentage_valid_models[,sort((1:n_datasets)[ntile(1:n_datasets,3) == i],decreasing = TRUE)] 
  
  clean_labels <- clean_dataset_names(colnames(input_selection))
  max_n_datasets <- ntile(1:n_datasets,3)  %>%
    table() %>% 
    max()
  n_datasets_selection <- ncol(input_selection)
  par(mar = c(margin_bottom_graphs,margin_dataset_labels,0,0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0,0), ylim = c(0,max_n_datasets*6))
  axis(2, at = seq(3.5, n_datasets_selection * 6 - 2.5, by = 6), labels = clean_labels, 
       las = 1, tick = FALSE, line = dataset_labels_position)
  
  par(mar = c(margin_bottom_graphs,0.1,0,0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0,100), ylim = c(0,max_n_datasets*6))
  axis(1, at = seq(0, 100, length.out = 5), labels = paste0(seq(0, 100, length.out = 5), "%"))
  barplot(input_selection, 
          width = 1,
          #space = 0,
          horiz = TRUE,
          beside = TRUE, 
          axes = FALSE,
          axisnames = FALSE,
          col = palette, 
          add = TRUE
  )
}


par(mar = c(0,0,0,0))
plot.new()
legend("bottom", 
       legend = clean_model_names(row.names(input_selection)), 
       fill = palette, 
       horiz = TRUE, inset = c(0, 0), bty = "n")




