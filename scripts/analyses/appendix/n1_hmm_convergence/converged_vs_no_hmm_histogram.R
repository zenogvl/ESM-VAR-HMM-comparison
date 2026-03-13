

var_path <- "Output/N1_Results/Rescaled_results/VAR_kalman_filter_RMSE/"

delta_rmse <- calc_delta_rmse(list(var_path, 
                                   "Output/N1_Results/Rescaled_results/HMM_kalman_filter_results/RMSE/"), 
                              dataset_names_all, 
                              c("var_kf", "hmm_kf"))
ds <- dataset_names_all[1]
var_rmse_converged <- list()
for(ds in dataset_names_all){
  id_converged_hmm <- delta_rmse %>%
    filter(dataset == ds) %>% 
    pull(ID)
  
  var_rmse_converged[[ds]] <- read.csv(paste0(var_path, "RMSE_", ds,".csv")) %>%
    mutate(mean_rmse = rowMeans(select(., -ID, -X), na.rm = TRUE)) %>%
    select(ID, mean_rmse) %>%
    mutate(converged_hmm = ifelse(ID %in% id_converged_hmm, "converged", "missing"),
           dataset = ds)
}
var_rmse_converged <- do.call(rbind, var_rmse_converged)

var_rmse_converged %>%
  filter(converged_hmm == "converged") %>% 
  summarise(invalid = sum(mean_rmse > 1), 
            total = length(mean_rmse))

var_rmse_converged %>%
  na.omit() %>%
  filter(converged_hmm == "missing") %>% 
  summarise(invalid = sum(mean_rmse > 1), 
            total = sum(converged_hmm == "missing"))



var_rmse_converged  %>%
  na.omit()



var_rmse_converged %>%
  na.omit() %>%
  filter(mean_rmse < 1) %>%
  group_by(converged_hmm) %>%
  summarise(mean = mean(mean_rmse),
            median = median(mean_rmse),
            sd = sd(mean_rmse))


pdf(file = "Output/Figures/AppendixC_histogram_rmse.pdf", 
    width = standard_figure_with,
    height = 3)


x_lim_start <- 0
x_lim_end <- 1
y_lim_end <- 120

layout(1)
par(mar = c(4, 4, 1, 1), 
    mgp = c(3, 2, 1))
plot.new()
plot.window(xlim = c(x_lim_start, x_lim_end), ylim = c(0, y_lim_end))
axis(1, at = seq(x_lim_start, x_lim_end, .1), labels = seq(x_lim_start, x_lim_end, .1))
axis(2, at = seq(0, y_lim_end, 20), labels =  seq(0, y_lim_end, 20))



var_rmse_converged %>%
  filter(converged_hmm == "converged" & mean_rmse < 1) %>% 
  na.omit() %>%
  pull(mean_rmse) %>%
  hist(breaks = seq(x_lim_start, x_lim_end, by = .005),
       xlim = c(x_lim_start + .01, x_lim_end - .01), 
       add = TRUE, 
       axes = FALSE, 
       col = adjustcolor("darkorange2", alpha.f = .5),
       border = "darkorange2"
  )

var_rmse_converged %>%
  filter(converged_hmm == "missing" & mean_rmse < 1) %>% 
  na.omit() %>%
  pull(mean_rmse) %>%
  hist(breaks = seq(x_lim_start, x_lim_end, by = .005),
       xlim = c(x_lim_start + .01, x_lim_end - .01), 
       add = TRUE, 
       axes = FALSE, 
       col = adjustcolor("turquoise3", alpha.f = .5),
       border = "turquoise3"
         )

legend("topright", 
       legend = c("Converged HMM", "No HMM"), 
       fill = c(adjustcolor("darkorange2", alpha.f = 0.5), adjustcolor("turquoise3", alpha.f = 0.5)), 
       border = c("darkorange2", "turquoise3"),
       bty = "n")


dev.off()




x_lim_start <- 0
x_lim_end <- 1
y_lim_end <- 120

layout(1)
par(mar = c(3, 3, 1, 1))
plot.new()
plot.window(xlim = c(x_lim_start, x_lim_end), ylim = c(0, y_lim_end))

# Extract data beforehand to adjust histogram parameters
hist1_data <- var_rmse_converged %>%
  filter(converged_hmm == "converged" & mean_rmse < 1) %>% 
  na.omit() %>%
  pull(mean_rmse)

hist2_data <- var_rmse_converged %>%
  filter(converged_hmm == "missing" & mean_rmse < 1) %>% 
  na.omit() %>%
  pull(mean_rmse)

# Compute histograms without plotting
hist1 <- hist(hist1_data, breaks = seq(x_lim_start, x_lim_end, by = .005), plot = FALSE)
hist2 <- hist(hist2_data, breaks = seq(x_lim_start, x_lim_end, by = .005), plot = FALSE)

# Plot first histogram with semi-transparency
plot(hist1, col = adjustcolor("darkorange2", alpha.f = 0.5), 
     border = "darkorange2", 
     xlim = c(x_lim_start + .01, x_lim_end - .01), 
     ylim = c(0, y_lim_end), 
     axes = FALSE, 
     main = "", 
     xlab = "Mean RMSE", 
     ylab = "Frequency")

# Overlay second histogram with transparency
plot(hist2, col = adjustcolor("turquoise3", alpha.f = 0.5), 
     border = "turquoise3", 
     add = TRUE)

axis(1, at = seq(x_lim_start, x_lim_end, .1), labels = seq(x_lim_start, x_lim_end, .1))
axis(2, at = seq(0, y_lim_end, 10), labels = seq(0, y_lim_end, 10))

