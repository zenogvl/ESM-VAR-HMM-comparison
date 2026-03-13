

max_vars <- 19


pdf(file = "Output/Figures/AppendixB_variables_vs_paramaters.pdf", 
    width = 5.9,
    height = 5)

layout(1)
par(mar = c(3,3,1,1))
plot.new()
plot.window(xlim = c(1,max_vars), ylim = c(0,400))
axis(1, 
     at = 1:max_vars, 
     labels = ifelse(1:max_vars %% 2 == 1, 1:max_vars, NA),
     cex.axis = 1
     ) 
axis(2, at = seq(0,400, 50), labels = seq(0,400, 50))
lines(x = 1:max_vars, y = n_paramaters_var(1:max_vars), 
      type = "l", lwd = 2, col = wesanderson::wes_palette("Zissou1", 5)[5])
lines(x = 1:max_vars, y = n_paramaters_hmm(1:max_vars,2), 
      type = "l", lwd = 2, col = wesanderson::wes_palette("Zissou1", 5)[1])
lines(x = 1:max_vars, y = n_paramaters_hmm(1:max_vars,3), 
      type = "l", lwd = 2, col = wesanderson::wes_palette("Zissou1", 5)[2])
lines(x = 1:max_vars, y = n_paramaters_hmm(1:max_vars,4), 
      type = "l", lwd = 2, col = wesanderson::wes_palette("Zissou1", 5)[3])
lines(x = 1:max_vars, y = n_paramaters_hmm(1:max_vars,5), 
      type = "l", lwd = 2, col = wesanderson::wes_palette("Zissou1", 5)[4])
legend("topleft", 
       legend = c("HMM S=2", "HMM S=3", "HMM S=4","HMM S=5", "VAR"), 
       fill =  wesanderson::wes_palette("Zissou1", 5), 
       bty = "n", 
       inset = c(0, 0.01)
       )

dev.off()


