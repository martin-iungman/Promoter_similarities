packages<-c("tidyverse","rtacklayer","BSgenome.Hsapiens.UCSC.hg19","plyr","furrr","tictoc")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

n.cores=parallel::detectCores()-2
inputFiles = list.files("Data/FANTOM5_files", full.names = T)[1:20]
plan(multisession,workers=n.cores)
bed_files<-future_map(inputFiles[1:20], ~import(.x,format="bed"))
plan("default")
sample_names<-gsub(".CNhs.+","",basename(inputFiles))%>%gsub( "( - )|(, )|( )", "_", .)
tissue<-gsub("_donor.","",sample_names)
donor<-str_extract(sample_names,"donor.")


tissue_beds<-vector(mode="list",length=length(unique(tissue)))
names(tissue_beds)<-unique(tissue)
for(i in unique(tissue)){
  group<-bed_files[grep(i,sample_names)] #split the bed_files in list according to the tissue
  samples=length(group)
  if(samples>1){
    u=unique(do.call("c",group)) #concatenate the GRanges and collapse
    mcols(u)=NULL
    for(k in 1:samples){ #for each one of the samples, the score will be add as metadata for the matching range, and when a range is missing in one sample, replace with 0
      mcols(u)[match(group[[k]], u), paste0("score_",k)] = group[[k]]$score
      mcols(u)[[paste0("score_",k)]]<-replace_na( mcols(u)[[paste0("score_",k)]],0)
    }
    mcols(u)$score<-rowSums(mcols(u)%>%as.matrix())  #sum the scores rowwise
    tissue_beds[[i]]<-u  #add the summarised GRange to the list
  }else{
    v<-group[[1]]
    # names(mcols(v))[2]<-"sum"
    tissue_beds[[i]]<-v
    } #if there was just one sample for the tissue, don't do anything
}

names(bed_files)<-sample_names
map_dfr(bed_files,~sum(mcols(.x)$score))%>%t()%>%as.data.frame()%>%rownames_to_column(var="sample_name")%>%dplyr::rename(seq_depth=V1)%>%mutate(tissue=str_remove(sample_name,"_donor."))%>%
  ggplot(aes(y=seq_depth, x=sample_name,fill=tissue))+geom_bar(stat="identity")+scale_y_log10()+theme(axis.text.x = element_text(size=5,angle=90),legend.position = "none")

map_dfr(tissue_beds,~sum(.x$score))%>%t()%>%as.data.frame()%>%rownames_to_column(var="sample_name")%>%dplyr::rename(seq_depth=V1)%>%mutate(tissue=str_remove(sample_name,"_donor."))%>%
  ggplot(aes(y=seq_depth, x=sample_name,fill=tissue))+geom_bar(stat="identity")+scale_y_log10()+theme(axis.text.x = element_text(size=5,angle=90),legend.position = "none")


h37_tss_selected_gr<-import.bed("Tables/h37_tss_selected_gr.bed")
find.tss.selected<-function(data,true.tss){
    a<-findOverlaps(true.tss, data, ignore.strand = F)
     data.frame(name= true.tss$name[queryHits(a)],
                Raw_count = mcols(data)$score[subjectHits(a)])%>%
      group_by(name) %>% dplyr::summarise(across(Raw_count,sum,.names="sum"))
}
plan(multisession,workers=n.cores)
tissue_beds_ol<-future_map(tissue_beds,~find.tss.selected(.x, h37_tss_selected_gr))
plan("default")
write_rds(tissue_beds_ol, "Tables/true_tss_activ_tissue.rds")
tissue_beds_ol2<-tissue_beds_ol%>%map2(names(tissue_beds_ol),~.x%>%dplyr::rename({{.y}}:=sum))
wide_table<-join_all(tissue_beds_ol2,by="name",type = "full")

#deberia hacer un heatmap aqui

