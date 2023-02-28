library(plyranges)
library(rtracklayer)
library(tidyverse)
library(foreach)
library(future)
dir.create("Tables/CAGE_filtered_bed")
df_genes<-read_csv("Tables/alt_prom_dataset.csv")
inputFiles = list.files("Data/FANTOM5_files", full.names = T)
library_size<-read_delim("Data/library_size.txt",delim=" ", col_names=F)
names(library_size)<-c("lib_size","file")
sample_names<-gsub(".CNhs.+","",library_size$file)%>%gsub( "(_-_)|(,_)", "_", .)
tissue<-gsub("(donor|pool|treated|rep|response).+","",sample_names)%>%gsub("_$","",.)%>%tolower()
library_size$file<-paste0("Data/FANTOM5_files/",library_size$file, ".gz")
library_size_merged<-library_size%>%mutate(tissue=tissue)%>%group_by(tissue)%>%summarise(lib_size=sum(lib_size),samples=n())
library_size_merged%>%arrange(desc(lib_size))

#selected<-"mast_cell"

myCluster <- makeCluster(30, # number of cores to use
                         type = "FORK") # type of cluster
registerDoParallel(myCluster)
foreach::foreach(selected=unique(library_size_merged$tissue)[7:9])%dopar%
{bed_files<-purrr::map(library_size$file[which(tissue==selected)], ~rtracklayer::import(.x,format="bed"))
if(length(bed_files)){
  u=unique(do.call("c",bed_files)) 
  mcols(u)=NULL
  for(k in 1:length(bed_files)){ 
    mcols(u)[match(bed_files[[k]], u), paste0("score_",k)] = bed_files[[k]]$score
    mcols(u)[[paste0("score_",k)]]<-replace_na( mcols(u)[[paste0("score_",k)]],0)
  }
  mcols(u)$score<-rowSums(mcols(u)%>%as.matrix()) 
}else{
  
  u<-bed_files[[1]]
} 
rm(bed_files)
gc()
prom_bed<-import.bed("Tables/final_database.bed")
mcols(prom_bed)<-mcols(prom_bed)%>%as_tibble()%>%select(name)
intersect<-join_overlap_intersect(u,prom_bed)
export.bed(intersect,paste0("Tables/CAGE_filtered_bed/",selected,"_filt.bed"))
rm(u)
gc()}
stopCluster(myCluster)
