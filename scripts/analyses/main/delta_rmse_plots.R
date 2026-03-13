
palette_boxplots <- wesanderson::wes_palette("Moonrise2", 3)

delta_rmse <- calc_delta_rmse(list("Output/N1_Results/Rescaled_results/VAR_kalman_filter_RMSE/", 
                                   "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/RMSE/"), 
                              dataset_names, 
                              c("var_kf", "hmm_kf"))

lowest_rmse <- find_lowest_rmse(path_list = list("Output/N1_Results/Rescaled_results/VAR_kalman_filter_RMSE/", 
                                                 "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/RMSE/"),
                                dataset_names = dataset_names, 
                                model_names = c("var", "hmm"), 
                                show_equal_results = FALSE,
                                format_output = "matrix")



###################### Figure 2 ########################

pdf(file = "Output/Figures/Figure_2_histogram_rmse.pdf", 
    width = standard_figure_with,
    height = 3)


x_lim_start <- -.2
x_lim_end <- .3
y_lim_end <- 120

layout(1)
par(mar = c(3, 3, 1, 1))
plot.new()
plot.window(xlim = c(x_lim_start, x_lim_end), ylim = c(0, y_lim_end))
axis(1, at = seq(x_lim_start, x_lim_end, .1), labels = seq(x_lim_start, x_lim_end, .1))
axis(2, at = seq(0, y_lim_end, 20), labels =  seq(0, y_lim_end, 20))


delta_rmse %>% 
  filter(delta_rmse > x_lim_start & delta_rmse < x_lim_end) %>% 
  pull(delta_rmse) %>%
  hist(breaks = seq(x_lim_start, x_lim_end, by = .002),
       xlim = c(x_lim_start + .01, x_lim_end - .01), 
       add = TRUE, 
       axes = FALSE, 
       col = palette_boxplots[1])
abline(v = 0, lty = 2, col = "red")
text(x = x_lim_start + 0.05, y = y_lim_end * 0.9, labels = "VAR beter", col = "black", cex = 1.2, adj = 0)  
text(x = x_lim_end - 0.05, y = y_lim_end * 0.9, labels = "HMM beter", col = "black", cex = 1.2, adj = 1) 


dev.off()


  
###################### Figure 3 ########################

pdf(file = "Output/Figures/Figure_3_barplot_lowest_rmse.pdf", 
    width = standard_figure_with,
    height = 6)


lowest_rmse_percentage_barplot(lowest_rmse[,order(lowest_rmse[1,])], 
                               palette = palette_boxplots[2:3],
                               split = 2, 
                               margin_bottom_graphs = 1.5, 
                               margin_dataset_labels = 6,
                               dataset_labels_position = -3.5,
                               labels_width = 1,
                               barplot_width = 1.3,
                               vline_colour = "red")





dev.off()

###################### Figure 4 ########################



datasets_order <- delta_rmse %>%
  filter(mean_var_kf < 1) %>%
  group_by(dataset) %>%
  summarise(mean_delta_rmse = mean(delta_rmse)) %>%
  ungroup() %>%
  mutate(clean_study_names = clean_dataset_names(dataset)) %>%
  arrange(mean_delta_rmse)

n_studies <- length(dataset_names)

pdf(file = "Output/Figures/Figure_4_boxplots.pdf", 
    width = standard_figure_with,
    height = 8.5)

layout(matrix(c(1:3, rep(4,3)), nrow = 2, ncol = 3, byrow = TRUE), 
       widths = c(2,3,3), 
       heights = c(20,.5))

amount_jitter <- .3
amount_alpha <- .3
points_size <- .3
start_ylim <- 1

par(mar = c(3, 10, 1, 1))
plot.new()
plot.window(xlim = c(0,0), ylim = c(start_ylim, n_studies*2))
axis(2, at = seq(1, n_studies*2, by = 2), labels = datasets_order$clean_study_names, las = 1)

par(mar = c(3, 0, 1, 1))
plot.new()
plot.window(xlim = c(0,.5), ylim = c(start_ylim, n_studies*2))
axis(1, at = seq(0,.5,.1), labels = seq(0,.5,.1))

for(i in 1:n_studies){
  data_study <- delta_rmse %>%
    filter(dataset ==  datasets_order$dataset[i])
  location_y_axis <- seq(1, n_studies*2, by = 2)[i]
  n_observations <- nrow(data_study)
  
  boxplot(data_study$mean_var_kf, 
          at = location_y_axis - .4,
          horizontal = TRUE, 
          add = TRUE, 
          outline = FALSE,
          axes = FALSE,
          frame.plot = FALSE,
          show.names = FALSE, 
          col = palette_boxplots[2]
  )
  points(data_study$mean_var_kf,
         jitter(rep(location_y_axis - .4, n_observations), amount = amount_jitter),
         pch = 1,
         cex = points_size,
         col = adjustcolor(palette_boxplots[2], alpha.f = amount_alpha)
  )
  boxplot(data_study$mean_hmm_kf, 
          at = location_y_axis + .4,
          horizontal = TRUE, 
          add = TRUE, 
          outline = FALSE,
          axes = FALSE,
          frame.plot = FALSE,
          show.names = FALSE, 
          col = palette_boxplots[3]
  )
  points(data_study$mean_hmm_kf,
         jitter(rep(location_y_axis + .4, n_observations), amount = amount_jitter),
         pch = 1,
         cex = points_size,
         col = adjustcolor(palette_boxplots[3], alpha.f = amount_alpha)
  )
  
}

x_lim_start <- -.2
x_lim_end <- .3

plot.new()
plot.window(xlim = c(x_lim_start, x_lim_end), ylim = c(start_ylim, n_studies*2))
axis(1, at = seq(x_lim_start, x_lim_end, .1), labels = seq(x_lim_start, x_lim_end, .1))
abline(v = 0, lty = 2, col = "black")

for(i in 1:n_studies){
  data_study <- delta_rmse %>%
    filter(dataset ==  datasets_order$dataset[i])
  location_y_axis <- seq(1, n_studies*2, by = 2)[i]
  n_observations <- nrow(data_study)
  
  boxplot(data_study$delta_rmse, 
          at = location_y_axis,
          horizontal = TRUE, 
          add = TRUE, 
          outline = FALSE,
          axes = FALSE,
          frame.plot = FALSE,
          show.names = FALSE, 
          col = palette_boxplots[1]
  )
  points(data_study$delta_rmse,
         jitter(rep(location_y_axis, n_observations), amount = amount_jitter),
         pch = 1,
         cex = points_size,
         col = adjustcolor(palette_boxplots[1], alpha.f = amount_alpha)
  )
}


par(mar = c(0,0,0,0))
plot.new()
legend("bottom", 
       legend = c( "VAR", "HMM", expression(Delta ~ "RMSE")), 
       fill = palette_boxplots[c(2,3,1)], 
       horiz = TRUE, inset = c(0, 0), bty = "n")





dev.off()


