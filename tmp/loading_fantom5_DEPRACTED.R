library(tidyverse)
library(rtracklayer)
library(CAGEr)
library(BSgenome.Hsapiens.UCSC.hg19)

join_prom_tables<-function(gene_promoter_table, data_table,identifier){
  data_table$name<-as.character(data_table$name)
  names(data_table)[2]<-"count"
  df<-left_join(gene_promoter_table,data_table,by=c("cluster_id"="name")) %>%
    mutate(count = replace_na(count, 0))  
  df<-group_by(df,by=gene_symbol) %>%
    mutate(sum_gen = sum(count,na.rm=T),
           prop =count/sum_gen)
  df <- df %>%
    ungroup() %>%
    select(-by)
  names(df)[length(names(df)) - (2:0)] <- paste0(identifier,"_",names(df)[length(names(df)) - (2:0)])
  return(df)}

inputFiles = list.files("Data/FANTOM5_files", full.names = T)
# ce <- CAGEexp( genomeName     = "BSgenome.Hsapiens.UCSC.hg19"
#                , inputFiles     = inputFiles
#                , inputFilesType = "ctss"
#                , sampleLabels   = gsub( ".CNhs.+", "", basename(inputFiles))%>%gsub( "( - )|(, )|( )", "_", .)
# )

bed_files<-map(inputFiles, ~import(.x,format="bed"))
sample_names<-gsub(".CNhs.+","",basename(inputFiles))%>%gsub( "( - )|(, )|( )", "_", .)
tissue<-gsub("_donor.","",sample_names)
donor<-str_extract(sample_names,"donor.")

bed_files_df<-bed_files%>%map(~as.data.frame(.x))
names(bed_files)<-sample_names
# bed_files_df<-pmap(list(bed_files_df,sample_names,tissue),~..1%>%mutate(sample_name=..2, tissue=..3))
seq_depth<-map_dfr(bed_files_df,~sum(.x$score))%>%t()%>%as.data.frame()%>%rownames_to_column(var="sample_name")%>%dplyr::rename(seq_depth=V1)%>%mutate(tissue=str_remove(sample_name,"_donor."))
ggplot(seq_depth, aes(y=seq_depth, x=sample_name,fill=tissue))+geom_bar(stat="identity")+scale_y_log10()+theme(axis.text.x = element_text(size=5,angle=90),legend.position = "none")

h37_tss_selected_gr<-import.bed("Tables/h37_tss_selected_gr.bed")
find.tss.selected<-function(data,true.tss,colname){
  a<-findOverlaps(true.tss, data, ignore.strand = F)
  data.frame(name= true.tss$name[queryHits(a)],
             Raw_count = mcols(data)$score[subjectHits(a)])%>%
    group_by(name) %>% summarise(across(Raw_count,sum,.names="sum_{colname}"))
}
bed_files_sum<-map2(bed_files, names(bed_files),~find.tss.selected(.x, h37_tss_selected_gr,.y))

clusters_selected_genes<-read_csv("Tables/clusters_selected_genes.csv")
wide.table_full<-clusters_selected_genes
for (i in seq(1,length(sample_names))) {
  wide.table_full=join_prom_tables(wide.table_full,bed_files_sum[[i]],sample_names[i])
  i=i+1
}
for(i in length(unique(tissue))){
  
}
test<-wide.table_full%>%rowwise()%>%mutate(test=across(matches(paste0("sum_",tissue[1])),sum))
¨