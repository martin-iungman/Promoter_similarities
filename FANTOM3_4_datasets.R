library(CAGEr)
library(FANTOM3and4CAGE)
library(rtracklayer)
library(tidyverse)
library(glue)

#cargo data de cage de FANTOM3 y 4
data(FANTOMtissueCAGEhuman)
names(FANTOMtissueCAGEhuman)
#excluyo los que estan asociados a enfermedad o no son control. 
data<-map2(FANTOMtissueCAGEhuman, names(FANTOMtissueCAGEhuman), ~.x%>%select(matches(str_glue("chr|pos|strand|{.y}|control"))))
for(i in 1:length(data)){
  if(any(names(data[[i]])=="control")){
    var=names(data)[i]
    data[[i]]<-data[[i]]%>%dplyr::rename({{var}}:=control)
  }
}
data<-keep(data,~ncol(.x)==4)
#veo profundidad de secuenciacion de cada tejido
seq_depth<-map_dfr(data,~sum(.x[[4]]))%>%t()%>%as.data.frame()%>%rownames_to_column(var="tissue")%>%rename(seq_depth=V1)
ggplot(seq_depth, aes(y=seq_depth, x=tissue,fill=tissue))+geom_bar(stat="identity")+scale_y_log10()+theme(axis.text.x = element_text(size=5,angle=90),legend.position = "none")

#hago overlap con la data de true tss
h37_tss_selected_gr<-import.bed("Tables/h37_tss_selected_gr.bed")
#hago un liftOver delhg19 al hg18
h18_tss_selected<-unlist(liftOver(h37_tss_selected_gr, chain=import.chain("Data/hg19ToHg18.over.chain")))

find.tss.selected<-function(data,true.tss,colname){
  a<-findOverlaps(true.tss, data, ignore.strand = F)
  data.frame(name= true.tss$name[queryHits(a)],
             Raw_count = mcols(data)$score[subjectHits(a)])%>%
    group_by(name) %>% summarise(across(Raw_count,sum,.names="sum_{colname}"))
}
data_gr<-map(data,~GRanges(seqnames=.x[[1]],range=IRanges(start = .x[[2]],width=1),strand=.x[[3]],score=.x[[4]]))
tss_sel<-map2(data_gr,names(data_gr),~find.tss.selected(.x,h18_tss_selected,.y))

#cant de tss distintos
map_dfr(tss_sel, nrow)%>%t()%>%as.data.frame()%>%rownames_to_column(var="tissue")%>%dplyr::rename(tss_n=V1)%>%
  ggplot(aes(tissue,tss_n, fill=tissue))+geom_bar(stat="identity")+scale_y_log10()+theme(axis.text.x = element_text(size=5,angle=90),legend.position = "none")
map_dfr(tss_sel, ~sum(.x[[2]]))%>%t()%>%as.data.frame()%>%rownames_to_column(var="tissue")%>%dplyr::rename(tss_n=V1)%>%
  ggplot(aes(tissue,tss_n, fill=tissue))+geom_bar(stat="identity")+scale_y_log10()+theme(axis.text.x = element_text(size=5,angle=90),legend.position = "none")

## unificar todos los df
prom_df<-read_csv("Tables/clusters_selected_genes.csv")%>%select(cluster_id,gene_symbol,distances)
tss_sel_df<-purrr::reduce(tss_sel, full_join, by="name")%>%
  mutate(across(.cols=where(is.numeric),.fns = ~replace_na(.x,replace=0)))
tss_sel_df<-left_join(tss_sel_df,prom_df, by=c("name"="cluster_id"))
tss_sel_df<-tss_sel_df%>%rowwise()%>%mutate(counts_tot=sum(c_across(starts_with("sum"))), n_samples=sum(c_across(starts_with("sum"))>1))%>%ungroup()
tss_sel_df%>%ggplot(aes(counts_tot))+geom_histogram()+scale_x_log10()+xlab("Counts totales")+geom_vline(xintercept = 10)
tss_sel_df%>%ggplot(aes(n_samples))+geom_histogram()+geom_vline(xintercept = 3,linetype="dashed")
tss_sel_df<-tss_sel_df%>%filter(counts_tot>=10, n_samples>=3)
#me quedo con los casos de n_prom=2 y reordeno
prom_list<-tss_sel_df%>%select(-n_samples,-counts_tot,-distances)%>%group_split(gene_symbol)%>%map(unique)%>%keep(~nrow(.x)==2)%>%
  map(~.x%>%pivot_longer(starts_with("sum"),names_to="tissue",values_to = "counts"))
prom_list<-map(prom_list, ~.x%>%pivot_wider(names_from = name,values_from = counts))
#calculo correlacion de spearman para cada uno
cor_vector<-map_dbl(prom_list,~cor(.x[[grep("^p", names(.x))[1]]],
    .x[[grep("^p", names(.x))[2]]],method = "spearman"))
names(cor_vector)<-map_chr(prom_list,~unique(.x$gene_symbol))

as.data.frame(cor_vector)%>%ggplot(aes(cor_vector))+geom_density()+xlab("correlacion de Spearman")
#cant de prom con los dos vectores con expresion en al menos un tejido (sino no se computa la correlacion)
length(which(!is.na(cor_vector)))

cor_df<-as.data.frame(cor_vector)%>%rownames_to_column(var="gene_symbol")%>%left_join(tss_sel_df, by="gene_symbol")
cor_df%>%ggplot(aes(distances, cor_vector))+geom_point()+scale_x_log10()+geom_smooth()
cor_df%>%ggplot(aes(counts_tot,cor_vector))+geom_point()+scale_x_log10()+geom_smooth(method = "gam")
cor_df%>%ggplot(aes(x=as_factor(n_samples),y=counts_tot))+geom_boxplot(outlier.alpha = 0)+geom_jitter(width=0.1,size=0.2)+scale_y_log10()
cor_sum<-cor_df%>%group_by(gene_symbol)%>%summarise(cor_vector, counts=sum(counts_tot))%>%unique()
cor_sum%>%  ggplot(aes(counts,cor_vector))+geom_point()+scale_x_log10()+geom_smooth(method = "gam")
cor(cor_sum$counts,cor_sum$cor_vector,method="spearman")

h18_filt<-subset(h18_tss_selected, name%in%unlist(prom_list%>%map(~names(.x)[-c(1,2)])))
mcols(h18_filt)<-mcols(h18_filt)%>%as_tibble()%>%mutate(gene=str_remove_all(name,"(p.@)|(,.+)"))
gr_list_genes<-split(h18_filt, h18_filt$gene)
dist_to_near_intra<-map(as.list(gr_list_genes),distanceToNearest)
dist_df<-data.frame(
  cluster_id=unlist(map2(as.list(gr_list_genes),dist_to_near_intra,~.x[queryHits(.y)]$name)),
  nearest_id=unlist(map2(as.list(gr_list_genes),dist_to_near_intra,~.x[subjectHits(.y)]$name)),
  distances=unlist(map(dist_to_near_intra, ~as_vector(mcols(.x)))))
cor_df%>%left_join(dist_df, by=c("name"="cluster_id"))%>%select(distances,cor_vector, gene_symbol)%>%unique()%>%
  ggplot(aes(distances, cor_vector))+geom_point(size=0.5)+geom_smooth()+scale_x_log10()

tss_sel_df%>%group_by(gene_symbol)%>%summarise(n_prom=n(),counts_tot,n_samples)%>%
  ggplot(aes(as_factor(n_prom),counts_tot))+geom_boxplot()+scale_y_log10()
tss_sel_df%>%group_by(gene_symbol)%>%summarise(n_prom=n(),counts_tot,n_samples)%>%
  ggplot(aes(as_factor(n_prom),n_samples))+geom_boxplot()

prom_df2<-prom_list2%>%map_dfr(~.x)
ggplot(prom_df2,aes(as_factor(n_samples),sum))+geom_boxplot()+scale_y_log10()
prom_list2<-tss_sel_df%>%select(-n_samples,-counts_tot,-distances)%>%group_split(gene_symbol)%>%map(unique)%>%
  map(~.x%>%pivot_longer(starts_with("sum"),names_to="tissue",values_to = "counts"))
prom_list2<-map(prom_list2,~.x%>%rowwise()%>%mutate(sum=sum(c_across(starts_with("p"))), n_prom_act=sum(c_across(starts_with("p"))>1)))


cor_vector2<-map_dbl(prom_list2,~cor(.x$sum,
                                   .x$n_prom_act,method = "spearman"))
names(cor_vector2)<-map_chr(prom_list2,~unique(.x$gene_symbol))
as.data.frame(cor_vector2)%>%ggplot(aes(cor_vector2))+geom_density()+xlab("correlacion de Spearman")
prom_list3<-keep(prom_list2,~ncol(.x)>=6)
cor_vector3<-cor_vector2[which(names(cor_vector2)%in%map_chr(prom_list3,~unique(.x$gene_symbol)))]
