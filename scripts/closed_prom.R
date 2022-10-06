library(rtracklayer)
library(tidyverse)
h37_tss_selected_gr<-import("Tables/h37_tss_selected_gr.bed")
#mcols(h37_tss_selected_gr)<-cbind(mcols(h37_tss_selected_gr), gene=str_remove_all(h37_tss_selected_gr$name,"(p.+@)|(,.+)"),width=width(h37_tss_selected_gr))

active_prom<-read_csv("Tables/active_prom_df.csv")

width_ic_inf= sort(active_prom$width)[nrow(active_prom)/40]
width_ic_sup=sort(active_prom$width)[nrow(active_prom)-nrow(active_prom)/40]
ggplot(active_prom, aes(width))+geom_density()+
  geom_vline(xintercept = c(width_ic_inf, width_ic_sup), linetype="dashed")
ggsave("Plots/width_distribution_IC.png")

closed_prom<-list()
for(i in 1:nrow(active_prom)){
  if(!is.na(active_prom[i,"distances"])){
  dist=active_prom[i,"width"]+active_prom[i,"distances"]+
    active_prom[which(active_prom$cluster_id==active_prom$nearest_id[i])[1],"width"]
  if(dist<width_ic_sup){
    closed_prom[[i]]<-c(unname(sort(unlist(c(active_prom[i,"cluster_id"],
               active_prom[i,"nearest_id"])))),
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
