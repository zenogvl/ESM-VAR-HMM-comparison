

clean_dataset_names <- function(dataset_names){
  dataset_names %>%
    gsub(pattern = "_", replacement = " ") %>%
    gsub(pattern = "(\\d{4})", replacement = "(\\1)") %>%
    gsub(pattern = "([a-z])([A-Z])", replacement = "\\1 \\2")  %>%
    gsub(pattern = ") ds1", replacement = "a)")%>%
    gsub(pattern = ") ds2", replacement = "b)")
}


clean_model_names <- function(model_names){
  model_names %>%
    gsub(pattern = "_", replacement = " ") %>% 
    gsub(pattern = "all equal", replacement = "Equal RMSE") %>%
    gsub(pattern = "\\b(\\w)", replacement = "\\U\\1", perl = TRUE) %>%
    gsub(pattern = "S(\\d{1})", replacement = "HMM S=\\1") %>%
    gsub(pattern = "var", replacement = "VAR", ignore.case = TRUE) %>%
    gsub(pattern = "hmm", replacement = "HMM", ignore.case = TRUE)
  
}

lowest_rmse_percentage_barplot <- function(matrix_input, 
                                           palette, 
                                           #clean_labels, 
                                           split = 1, 
                                           margin_bottom_graphs = 2, 
                                           margin_dataset_labels = 8, 
                                           dataset_labels_position = -1, 
                                           labels_width = 1,
                                           barplot_width = 3, 
                                           vline_colour = "red"){
  
  n_bars = ncol(matrix_input)
  
  if(split == 2){
    layout(matrix(c(1:4,rep(5, 4)), nrow = 2, byrow = TRUE), 
           widths = c(labels_width, barplot_width, labels_width, barplot_width),
           heights = c(12,1))
    par(oma = c(.5, .5, .5, 1.5)) 
    
    matrix_input_1 <- matrix_input[,as.integer(n_bars/2):1]
    
    
  } else {
    
    layout(matrix(c(1:2,rep(3, 2)), nrow = 2, byrow = TRUE), 
           widths = c(labels_width, barplot_width),
           heights = c(12,1))
    par(oma = c(1, 1, 1, 1)) 
    matrix_input_1 <- matrix_input
    
  }
  
  clean_labels <- clean_dataset_names(colnames(matrix_input_1))
  n_bars_1 <- ncol(matrix_input_1)
  par(mar = c(margin_bottom_graphs,margin_dataset_labels,0,0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0,0), ylim = c(0,n_bars_1))
  axis(2, at = (1:n_bars_1 -.5), labels = clean_labels, 
       las = 1, tick = FALSE, line = dataset_labels_position)
  
  par(mar = c(margin_bottom_graphs,0.1,0,0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0,1), ylim = c(0,n_bars_1))
  axis(1, at = seq(0, 1, length.out = 5), labels = paste0(seq(0, 100, length.out = 5), "%"))
  barplot(matrix_input_1, 
          width = 1,
          space = 0,
          horiz = TRUE,
          beside = FALSE, 
          axes = FALSE,
          axisnames = FALSE,
          col = palette, 
          add = TRUE
  )
  abline(v=.5, col = vline_colour)
  
  if(split == 2){
    matrix_input_2 <- matrix_input[,n_bars:as.integer(n_bars/2)]
    clean_labels <- clean_dataset_names(colnames(matrix_input_2))
    
    n_bars_2 <- ncol(matrix_input_2)
    par(mar = c(margin_bottom_graphs,margin_dataset_labels,0,0), xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(0,0), ylim = c(0,n_bars_2))
    axis(2, at = (1:n_bars_2 -.5), labels = clean_labels, 
         las = 1, tick = FALSE, line = dataset_labels_position)
    
    par(mar = c(margin_bottom_graphs,0.1,0,0), xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(0,1), ylim = c(0,n_bars_2))
    axis(1, at = seq(0, 1, length.out = 5), labels = paste0(seq(0, 100, length.out = 5), "%"))
    barplot(matrix_input_2, 
            width = 1,
            space = 0,
            horiz = TRUE,
            beside = FALSE, 
            axes = FALSE,
            axisnames = FALSE,
            col =palette, 
            add = TRUE
    )
    abline(v=.5, col = vline_colour)
  }
  par(mar = c(0,0,0,0))
  plot.new()
  legend("bottom", 
         legend = clean_model_names(row.names(matrix_input)), 
         fill = palette, 
         horiz = TRUE, inset = c(0, 0), bty = "n")
  
  
}




percentage_models_barplot <- function(
    percantages_matrix,
    n_datasets, 
    n_models,
    palette, 
    margin_bottom_graphs = 1.5,
    margin_dataset_labels = 8,
    dataset_labels_position = -1,
    labels_width = 1.4,
    barplot_width = 1,
    split = 3
) {
  layout(matrix(c(1:(split*2),rep(split*2+1, split*2)), nrow = 2, byrow = TRUE), 
         widths = rep(c(labels_width, barplot_width), split),
         heights = c(12,1))
  par(oma = c(1, 1, 1, 2)) 
  
  
  for(i in 1:split){
    input_selection <- percantages_matrix[,sort((1:n_datasets)[ntile(1:n_datasets,3) == i],decreasing = TRUE)] 
    
    clean_labels <- clean_dataset_names(colnames(input_selection))
    max_n_datasets <- ntile(1:n_datasets,3)  %>%
      table() %>% 
      max()
    n_datasets_selection <- ncol(input_selection)
    par(mar = c(margin_bottom_graphs,margin_dataset_labels,0,0), xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(0,0), ylim = c(0,max_n_datasets*n_models))
    axis(2, 
         at = seq((n_models + 1)/2, n_datasets_selection * n_models - (n_models - 1)/2, by = n_models), 
         labels = clean_labels, 
         las = 1, 
         tick = FALSE, 
         line = dataset_labels_position)
    
    par(mar = c(margin_bottom_graphs,0.1,0,0), xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(0,100), ylim = c(0,max_n_datasets*n_models))
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
  
}


