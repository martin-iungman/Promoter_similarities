setwd("/home/tinchis/Documents/Promoter architecture")
selected_samples <- CAGEr::FANTOM5humanSamples[-grep("cell line", 
                                                     CAGEr::FANTOM5humanSamples[,"type"]),"data_url"]
selected_samples<-stringr::str_replace(selected_samples, "http","https")
writeLines(selected_samples[10:50],"Data/FANTOM5_files.txt")
