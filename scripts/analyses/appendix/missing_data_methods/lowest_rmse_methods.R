

var_lowest_rmse <- find_lowest_rmse(path_list = list("Output/N1_Results/Rescaled_results/VAR_kalman_filter_RMSE/",
                                                     "Output/N1_Results/Rescaled_results/VAR_listwise_deletion_RMSE/"),
                                    dataset_names = dataset_names, 
                                    model_names = c("kalman_filter", "listwise_deletion"), 
                                    show_equal_results = TRUE,
                                    format_output = "matrix")

palette <- wesanderson::wes_palette("Darjeeling1", 3)



pdf(file = "Output/Figures/AppendixA_var_methods_lowest_rmse.pdf", 
    width = standard_figure_with,
    height = 5)

lowest_rmse_percentage_barplot(var_lowest_rmse[,order(var_lowest_rmse[1,], decreasing = TRUE)], 
                               palette = wesanderson::wes_palette("Darjeeling1", 3), 
                               split = 2, 
                               margin_bottom_graphs = 1.5, 
                               margin_dataset_labels = 6,
                               dataset_labels_position = -3.5,
                               labels_width = 1,
                               barplot_width = 1.3,
                               vline_colour = "black")

dev.off()


hmm_lowest_rmse <- find_lowest_rmse(path_list = list("Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/RMSE/",
                                                     "Output/N1_Results/Rescaled_results/HMM_no_update_results/RMSE/",
                                                     "Output/N1_Results/Rescaled_results/HMM_seperate_series_results/RMSE/"),
                                    dataset_names = dataset_names, 
                                    model_names = c("kalman_filter", "continued_transitions", "seperate__time_series"), 
                                    show_equal_results = FALSE,
                                    format_output = "matrix")




pdf(file = "Output/Figures/AppendixA_hmm_methods_lowest_rmse.pdf", 
    width = standard_figure_with,
    height = 5)


lowest_rmse_percentage_barplot(hmm_lowest_rmse[,order(hmm_lowest_rmse[1,], decreasing = TRUE)], 
                               palette = wesanderson::wes_palette("Moonrise3", 3), 
                               split = 2, 
                               margin_bottom_graphs = 1.5, 
                               margin_dataset_labels = 6,
                               dataset_labels_position = -3.5,
                               labels_width = 1,
                               barplot_width = 1.3,
                               vline_colour = "black")




dev.off()


layout(matrix(1:20, nrow = 4, byrow = TRUE))

palette_names <- c("BottleRocket1", "BottleRocket2", "Rushmore1", "Royal1", "Royal2", 
                   "Zissou1", "Darjeeling1", "Darjeeling2", "Chevalier1" , "FantasticFox1" , "Moonrise1", 
                   "Moonrise2", "Moonrise3", "Cavalcanti1", "GrandBudapest1", "GrandBudapest2", "IsleofDogs1", 
                   "IsleofDogs2", "FrenchDispatch", "AsteroidCity2", "AsteroidCity2", "AsteroidCity3")


wesanderson::wes_palette(palette_names[1])
wesanderson::wes_palette(palette_names[2])
wesanderson::wes_palette(palette_names[3])
wesanderson::wes_palette(palette_names[4])
wesanderson::wes_palette(palette_names[5])
wesanderson::wes_palette(palette_names[6])
wesanderson::wes_palette(palette_names[7])
wesanderson::wes_palette(palette_names[8])
wesanderson::wes_palette(palette_names[9])
wesanderson::wes_palette(palette_names[10])
wesanderson::wes_palette(palette_names[11])
wesanderson::wes_palette(palette_names[12])
wesanderson::wes_palette(palette_names[13])
wesanderson::wes_palette(palette_names[14])
wesanderson::wes_palette(palette_names[15])
wesanderson::wes_palette(palette_names[16])
wesanderson::wes_palette(palette_names[17])
wesanderson::wes_palette(palette_names[18])
wesanderson::wes_palette(palette_names[19])
wesanderson::wes_palette(palette_names[20])

wesanderson::wes_palette(palette_names[22])

graphics.off() 