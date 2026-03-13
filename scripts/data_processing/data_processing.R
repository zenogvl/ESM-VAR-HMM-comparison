


f <- list.files("Scripts/data_processing/")

metaData <- vector("list", length(f))
names(metaData) <- gsub(".R", "", gsub("Import_", "", f))
variables_list <- vector("list", length(f))
names(variables_list) <- gsub(".R", "", gsub("Import_", "", f))

for(i_proc in 1:length(f)){
  print(f[i_proc])
  #Run scripts to create clean data
  source(paste0("Scripts/data_processing/", f[i_proc]))
  
  
  #Get variables vector 
  #metaData[[i_proc]]$variables <- Emotions
  #metaData[[i_proc]]$dsName <- gsub(".R", "", gsub("Import_", "", f[i_proc]))
  variables_list[[i_proc]] <- Emotions
  
}
saveRDS(variables_list, "Data/variables_list.rds")



saveRDS(metaData, "Data/metaDateRFile.rds")

