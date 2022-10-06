packages<-c("tidyverse","rtracklayer","BSgenome.Hsapiens.UCSC.hg19","plyr","furrr","tictoc","foreach")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
find.tss.selected<-function(data,true.tss){
    a<-findOverlaps(true.tss, data, ignore.strand = F)
     data.frame(name= true.tss$name[queryHits(a)],
                Raw_count = mcols(data)$score[subjectHits(a)])%>%
      group_by(name) %>% dplyr::summarise(across(Raw_count,sum,.names="sum"))
}
h37_tss_selected_gr<-import.bed("Tables/h37_tss_selected_gr.bed")
inputFiles = list.files("Data/FANTOM5_files", full.names = T)
sample_names<-gsub(".CNhs.+","",basename(inputFiles))%>%gsub( "( - )|(, )|( )", "_", .)
tissue<-gsub("(donor|pool|treated|rep|response).+","",sample_names)%>%gsub("_$","",.)%>%tolower()
#donor<-str_extract(sample_names,"donor.")
n.cores=5
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK",
  outfile=""
)
print("initializing analysis")
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()



ol_list<-foreach(i=unique(tissue),.packages=c("magrittr","dplyr","tidyr","rtracklayer"))%dopar%{
  ol<-list()
  bed_files<-purrr::map(inputFiles[which(tissue==i)], ~rtracklayer::import(.x,format="bed"))
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
  ol[[i]]<-find.tss.selected(u, h37_tss_selected_gr)
  rm(bed_files,u)
  gc()
  return(ol)
}
parallel::stopCluster(my.cluster)
ol_list<-flatten(ol_list)
print("writing output")
write_rds(ol_list,"Tables/true_tss_activ_tissue.rds")

	

