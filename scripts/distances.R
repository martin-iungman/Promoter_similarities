packages<-c("rtracklayer","stringi","tidyverse")
invisible(lapply(packages, library, character.only = TRUE))

h37_tss_selected_gr<-import("Tables/h37_tss_selected_gr.bed")
mcols(h37_tss_selected_gr)<-cbind(mcols(h37_tss_selected_gr), gene=str_remove_all(h37_tss_selected_gr$name,"(p.+@)|(,.+)"))

clusters_selected_genes<-read_csv("Tables/clusters_selected_genes.csv")
active_prom<-read_tsv("Tables/active_prom.tsv")
clusters_selected_genes<-clusters_selected_genes%>%filter(cluster_id%in%active_prom[[1]])

h37_tss_selected_gr <- subset(h37_tss_selected_gr, name %in% active_prom[[1]])

# get cluster width
mcols(h37_tss_selected_gr)<-cbind(mcols(h37_tss_selected_gr),width=width(h37_tss_selected_gr))
clusters_selected_genes<-mcols(h37_tss_selected_gr)%>%as_tibble()%>%select(name,width)%>%left_join(clusters_selected_genes,.,by=c("cluster_id"="name"))

##width distribution
ggplot(clusters_selected_genes, aes(width))+geom_density()
ggsave("Plots/width_distribution.png")

#distancia al prom mas cercano
print("Calculating distances")
gr_list_genes<-split(h37_tss_selected_gr, h37_tss_selected_gr$gene)
dist_to_near_intra<-map(as.list(gr_list_genes),distanceToNearest)
dist_df<-data.frame(
  cluster_id=unlist(map2(as.list(gr_list_genes),dist_to_near_intra,~.x[queryHits(.y)]$name)),
  nearest_id=unlist(map2(as.list(gr_list_genes),dist_to_near_intra,~.x[subjectHits(.y)]$name)),
  distances=unlist(map(dist_to_near_intra, ~as_vector(mcols(.x)))))

clusters_selected_genes<-left_join(clusters_selected_genes, dist_df, by=c("cluster_id"))

ggplot(clusters_selected_genes, aes(distances))+geom_density()+scale_x_log10()
ggsave("Plots/distance_distribution.png")

#dist_to_near<-distanceToNearest(h37_tss_selected_gr)
#dist_df2<-data.frame(cluster_id=h37_tss_selected_gr[queryHits(dist_to_near),]$name,
#                     nearest_id_extragen=h37_tss_selected_gr[subjectHits(dist_to_near),]$name,
#                     distance_extragen=as_vector(mcols(dist_to_near)))
#clusters_selected_genes<-left_join(clusters_selected_genes, dist_df2, by="cluster_id")

#ggplot(clusters_selected_genes, aes(distances_extragen))+geom_density()+scale_x_log10()

write_csv(clusters_selected_genes,"Tables/active_prom_df.csv")
