library(tidyverse)
library(rtracklayer)
library(plyranges)

inputFiles<-list.files("Tables/CAGE_filtered_bed")
samples<-sub("_filt.bed","",inputFiles)%>%sub("^.+bed//","",.)

bed_files<-map(inputFiles,~import.bed(paste0("Tables/CAGE_filtered_bed/",.x)))

u=unique(do.call("c",bed_files)) 
mcols(u)=NULL
for(k in 1:length(bed_files)){ 
  mcols(u)[match(bed_files[[k]], u), paste0("score_",samples[k])] = bed_files[[k]]$score
  mcols(u)[[paste0("score_",samples[k])]]<-replace_na( mcols(u)[[paste0("score_",samples[k])]],0)
}
mcols(u)$score<-rowSums(mcols(u)%>%as.matrix()) 
prom_bed<-import.bed("Tables/final_database.bed")
mcols(prom_bed)<-mcols(prom_bed)%>%as_tibble()%>%select(name)
intersect<-join_overlap_intersect(u,prom_bed)
mcols(intersect)$n.samples<-apply(mcols(intersect)%>%as_tibble()%>%select(starts_with("score_")), 1,function(x) length(which(x>0)))
write_tsv(intersect,"Tables/full_all_samples_CAGE.tsv")
mcols(intersect)<-mcols(intersect)%>%as_tibble()%>%select(name,score,n.samples)
export.bed(intersect,"Tables/CAGE_filtered_bed/all_samples_filt.bed")
write_tsv(intersect,"Tables/all_samples_CAGE.tsv")
