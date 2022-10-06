packages<-c("tidyverse","furrr","doParallel")
lapply(packages,library,character.only = TRUE)
tissue_rds<-readRDS("Tables/true_tss_activ_tissue.rds")
clusters_selected_genes<-read_csv("Tables/clusters_selected_genes.csv")
n.cores <- 20

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK",
  outfile=""
)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

tissue_rds<-future_map(tissue_rds,~full_join(.x,data.frame(name=unique(clusters_selected_genes$cluster_id)))%>%dplyr::mutate(sum=tidyr::replace_na(sum,0)),.progress=T)



#density pĺot de cantidad de reads por tissue (agrupadas las replicas)
totcounts_per_tissue<-map_dbl(tissue_rds,~sum(.x$sum))
data.frame(total_counts=totcounts_per_tissue)%>%
  ggplot(aes(total_counts))+
  geom_density()+scale_x_log10()+
  geom_point(aes(x=total_counts,y=runif(length(totcounts_per_tissue),min=0,max=0.2)),size=0.5)+
  ylab("density")
ggsave("Plots/total_counts_per_tissue.png")

#La normalizacion de los counts la hago dividiendo los counts por promotor/"library size" del tissue. Onda cpm.
#Density plot de la suma de normalizados para cadea promotor (cada punto es un promotor)
map2(tissue_rds,totcounts_per_tissue, ~mutate(.x,sum=10E6*sum/.y))%>%
  map2_dfr(names(tissue_rds),~mutate(.x,tissue=.y))%>%
  group_by(name)%>%summarise(total_prom=sum(sum))%>%ggplot(aes(total_prom+0.5))+geom_density()+scale_x_log10()+xlab("Sum of normalized counts across all samples")
ggsave("Plots/total_normalized_counts_per_prom.png")

tissue_long<-map2_dfr(tissue_rds,names(tissue_rds),~mutate(.x,tissue=.y))
#Lo mismo de antes pero sin normalizacion
tissue_long%>%group_by(name)%>%summarise(total_prom=sum(sum))%>%ggplot(aes(total_prom+0.5))+geom_density()+scale_x_log10()
ggsave("Plots/total_counts_per_prom.png")

#tissue_long%>%group_by(name)%>%count(vars = sum>10)%>%arrange(desc(n))
#tissue_long_n.act<-map_dfr(seq(0,25,by=5), ~tissue_long%>%group_by(name)%>%count(vars = sum>.x)%>%mutate(threshold=as_factor(.x)))
tissue_long_n.act<-map_dfr(seq(0,25,by=5), ~tissue_long%>%filter(sum>=.x)%>%group_by(name)%>%count(.drop=F)%>%mutate(threshold=as_factor(.x)))
tissue_long_n.act%>%ggplot(aes(n,col=threshold))+geom_histogram(binwidth=1)+xlab("Number of tissue where promoter is active (count>threshold)")
ggsave("Plots/threshold_prom_activity_definition.png")

# que tenga al menos 5 counts en 5 tejidos
tissue_activity_filt<-tissue_long_n.act%>%filter(threshold==5&n>=5)

tissue_activity_filt%>%group_by(threshold)%>%count()
tissue_activity_filt_name<-unique(tissue_activity_filt$name)
length(tissue_activity_filt_name)
write_tsv(tissue_activity_filt[,1],"Tables/active_prom.tsv")

#por otro lado, tomo los que tienen alta actividad (>25counts) en un prom pero no en mas de 5. Estos los voy a considerar promotores pero no los voy a intentar correlacionar con sus vecinos (daria todo NA y la correlacion seria medio invalida)
tissue_restricted_prom<-tissue_long_n.act%>%filter(threshold==25)%>%filter(!name%in%tissue_activity_filt$name)
nrow(tissue_restricted_prom)
write_tsv(tissue_restricted_prom,"Tables/tissue_restricted_prom.tsv")
