library(tidyverse)
library(CAGEr)
library(furrr)
library(rtracklayer)
data(FANTOM5humanSamples)
View(FANTOM5humanSamples)

selected_samples <- FANTOM5humanSamples[-grep("cell line", 
                                             FANTOM5humanSamples[,"type"]),]
# astrocyteCAGEset <- importPublicData(source = "FANTOM5", dataset = "human", 
#                                      sample = selected_samples[1:3,"sample"])
selected_samples[[5]]<-str_replace(selected_samples[[5]], "^http","https")
n.cores=availableCores()-2
plan(multisession, workers=n.cores)
data<-future_map(seq(1,n.cores), ~import(selected_samples[.x,5]))

seq_depth<-map_dfr(data,~sum(.x$score))%>%t()%>%as.data.frame()%>%rownames_to_column(var="tissue")%>%rename(seq_depth=V1)
ggplot(seq_depth, aes(y=seq_depth, x=tissue,fill=tissue))+geom_bar(stat="identity")+scale_y_log10()+theme(axis.text.x = element_text(size=5,angle=90),legend.position = "none")
