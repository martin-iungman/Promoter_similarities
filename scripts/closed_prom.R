library(rtracklayer)
library(tidyverse)
h37_tss_selected_gr<-import("Tables/h37_tss_selected_gr.bed")
#mcols(h37_tss_selected_gr)<-cbind(mcols(h37_tss_selected_gr), gene=str_remove_all(h37_tss_selected_gr$name,"(p.+@)|(,.+)"),width=width(h37_tss_selected_gr))

clusters_selected_genes<-read_csv("Tables/clusters_selected_genes.csv")
#clusters_selected_genes<-mcols(h37_tss_selected_gr)%>%as_tibble()%>%select(name,width)%>%left_join(clusters_selected_genes,.,by=c("cluster_id"="name"))
width_ic_inf= sort(clusters_selected_genes$width)[nrow(clusters_selected_genes)/40]
width_ic_sup=sort(clusters_selected_genes$width)[nrow(clusters_selected_genes)-nrow(clusters_selected_genes)/40]
ggplot(clusters_selected_genes, aes(width))+geom_density()+
  geom_vline(xintercept = c(width_ic_inf, width_ic_sup), linetype="dashed")

closed_prom<-list()
for(i in 1:nrow(clusters_selected_genes)){
  if(!is.na(clusters_selected_genes[i,"distances"])){
  dist=clusters_selected_genes[i,"width"]+clusters_selected_genes[i,"distances"]+
    clusters_selected_genes[which(clusters_selected_genes$cluster_id==clusters_selected_genes$nearest_id[i])[1],"width"]
  if(dist<width_ic_sup){
    closed_prom[[i]]<-c(unname(sort(unlist(c(clusters_selected_genes[i,"cluster_id"],
               clusters_selected_genes[i,"nearest_id"])))),
               dist)
    names(closed_prom[[i]])<-c("prom1","prom2","dist")
    }

  }
}

# De los promotores totales seleccionaddos, cuantos tienen otro del mismo gen a menos de 71pb?
compact(closed_prom)%>%length()
#closed_prom_df<-compact(closed_prom%>%map(~.x%>%unlist()))%>%map_dfr(~.x)#%>%mutate(proms=sort(c(cluster_id,nearest_id)))
closed_prom_df<-closed_prom%>%map_dfr(~.x)
closed_prom_df<-closed_prom_df%>%unique()
write_rds(closed_prom_df, "Tables/proximal_promoters.rds")
