selected_samples <- CAGEr::FANTOM5humanSamples[-grep("cell line", 
                                                     CAGEr::FANTOM5humanSamples[,"type"]),]
selected_samples[["data_url"]]<-stringr::str_replace(selected_samples[["data_url"]], "http","https")
write.table(selected_samples,"Data/FANTOM5_files_info.tsv")
writeLines(selected_samples[["data_url"]],"Data/FANTOM5_files.txt")
