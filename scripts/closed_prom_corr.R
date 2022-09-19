packages<-c("tidyverse","rtacklayer","furrr","tictoc","doParallel","foreach")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))
proximal_promoters<-readRDS("Tables/proximal_promoters.rds")
clusters_selected_genes<-read_csv("Tables/clusters_selected_genes.csv")
tissue_beds_ol<-read_rds("Tables/true_tss_activ_tissue.rds")
n.cores <- parallel::detectCores() - 2

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK",
  outfile=""
)

tissue_beds_ol<-tissue_beds_ol%>%map(~full_join(.x,clusters_selected_genes[,c(1,2)],by=c("name"="cluster_id"))%>%mutate(sum=replace_na(sum,0)))
tissue_genes<-tissue_beds_ol%>%map(~.x%>%group_by(gene_symbol)%>%dplyr::summarise(gene_sum=sum(sum)))

#ggplot()+geom_density(data = tissue_genes[[1]],mapping=aes(gene_sum+0.9),col="blue")+
#  geom_density(data = tissue_genes[[2]],mapping=aes(gene_sum+0.9),col="red")+
#  geom_density(data = tissue_genes[[3]],mapping=aes(gene_sum+0.9),col="black")+
#  geom_density(data = tissue_genes[[4]],mapping=aes(gene_sum+0.9),col="green")+
#  geom_density(data = tissue_genes[[5]],mapping=aes(gene_sum+0.9),col="darkred")+
#  geom_density(data = tissue_genes[[6]],mapping=aes(gene_sum+0.9),col="darkred")+
#  geom_density(data = tissue_genes[[7]],mapping=aes(gene_sum+0.9),col="darkred")+
#  geom_density(data = tissue_genes[[8]],mapping=aes(gene_sum+0.9),col="darkred")+
#  scale_x_log10()

#disribucion de los counts de los promotores

tissue_beds_ol%>%map_dfr(~.x[sample(1:nrow(.x),size=10000),])%>%mutate(sum2=sample(sum))%>%
  ggplot(aes(log10(sum+sum2)))+geom_density()+xlim(0,4)
ggsave("prom_activ_density.png")

###### deberia agrupar en listas los pares y nombrarlos como a los genes. Asi ya tengo facil para hacer otra ronda
# prox_prom_list<-vector(mode="list",length=nrow(split_proximal)) #reemplazar el length a  nrow(proximal_promoters)
# tic()
# for(k in 1:nrow(proximal_promoters)){ #for each one of the pairs, I'll repeat the sequence
#   df<-data.frame(tissue=unique(tissue), p1=NA,p2=NA) #preallocation
#   p1=proximal_promoters[[k,1]]
#   p2=proximal_promoters[[k,2]]
#   sym<-tissue_beds_ol[[1]][[grep(paste0("^",p1,"$"),tissue_beds_ol[[1]]$name)[1],3]] #find the name of the respective gene symbol
#   for(i in 1:length(unique(tissue))){
#     expr_gene<-tissue_genes[[i]][[grep(paste0("^",sym,"$"),tissue_genes[[i]]$gene_symbol),2]] #find the expression of the whole gene
#     if(expr_gene>10){ #only consider samples where the gene expression is relevant (value to define...how?)
#       df[i,2]<-tissue_beds_ol[[i]][[grep(p1,tissue_beds_ol[[i]]$name)[1],2]]  #counts for one of the promoters in tissue i
#       df[i,3]<-tissue_beds_ol[[i]][[grep(p2,tissue_beds_ol[[i]]$name)[1],2]]  #counts for the other one in tissue i
#       if(df[i,2]==0&df[i,3]==0){
#         df[i,2]=df[i,3]=NaN}else if((df[i,2]==0)|(df[i,3]==0)){
#         df[i,2]=df[i,2]+0.5
#         df[i,3]=df[i,3]+0.5
#         }
#       }else{
#         df[i,2]<-NaN
#         df[i,3]<-NaN
#       }
#   }
#     prox_prom_list[[k]][[1]]<-df
#     prox_prom_list[[k]][[2]]<-c(sym,p1,p2)
#     prox_prom_list[[k]][[3]]<-log2(df$p1/df$p2)
#     prox_prom_list[[k]][[4]]<-var(prox_prom_list[[k]][[3]],na.rm=T)
# }
# toc()
# prox_prom_list<-compact(prox_prom_list)
# prox_prom_list%>%map_dbl(~.x[[4]])%>%as_data_frame()%>%filter(!is.na(value),abs(value)!=100)%>%nrow()
# prox_prom_list%>%map_dbl(~.x[[4]])%>%as_data_frame()%>%filter(!is.na(value),abs(value)!=100)%>%ggplot(aes(value))+geom_density()+scale_x_log10()

#### le agrgaria otra condicion para que los counts de p1+p2 son menores a un umbral, no considere (ademas de los del gen total). Principalmente para no distorsionar tanto con grandes cocientes entre valores bajos
# 
# n.cores <- parallel::detectCores() - 2
# 
# my.cluster <- parallel::makeCluster(
#   n.cores, 
#   type = "PSOCK"
# )
# split_proximal<-left_join(proximal_promoters,clusters_selected_genes[,c(1,2)],by=c("prom1"="cluster_id"))%>%distinct()%>%group_split(gene_symbol)
# names(split_proximal)<-map_chr(split_proximal,~unique(.x$gene_symbol) ) 
# 
# prox_prom_list<-vector(mode="list",length=length(split_proximal)) #reemplazar el length a  nrow(proximal_promoters)
# names(prox_prom_list)<-names(split_proximal)
# doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
# foreach::getDoParWorkers()
# tic()
# prox_prom_list<-foreach(gene=1:length(split_proximal))%dopar%{
#   prox_prom_list[[gene]]<-vector(mode="list",length=nrow(split_proximal[[gene]]))
# 
# for(k in 1:nrow(split_proximal[[gene]])){ 
#   df<-data.frame(tissue=unique(tissue), p1=NA,p2=NA) #preallocation
#   p1=split_proximal[[gene]][[k,1]]
#   p2=split_proximal[[gene]][[k,2]]
#   sym<-names(split_proximal)[gene] #find the name of the respective gene symbol
#   for(i in 1:length(unique(tissue))){
#     expr_gene<-tissue_genes[[i]][[grep(paste0("^",sym,"$"),tissue_genes[[i]]$gene_symbol),2]] #find the expression of the whole gene
#     if(expr_gene>10){ #only consider samples where the gene expression is relevant (value to define...how?)
#       df[i,2]<-tissue_beds_ol[[i]][[grep(p1,tissue_beds_ol[[i]]$name)[1],2]]  #counts for one of the promoters in tissue i
#       df[i,3]<-tissue_beds_ol[[i]][[grep(p2,tissue_beds_ol[[i]]$name)[1],2]]  #counts for the other one in tissue i
#       if(df[i,2]==0&df[i,3]==0){ #if no counts for either promoter, return NA
#         df[i,2]=df[i,3]=NaN}else if((df[i,2]==0)|(df[i,3]==0)){
#           df[i,2]=df[i,2]+0.5
#           df[i,3]=df[i,3]+0.5
#         } #if any of the promoters has 0 counts, add a pseudocount. The choice of 0.5 was to add the unity (in this case 1, as they are absolute counts) divided by the numbers of groups (2 in this case. actually they must be the n_prom of the gene, but it was more computing...) It will make a difference if the other promoter has just 1 count (1/0.5=2) but won't be of great magnitude if the value is rather small (5/0.5=10). The  
#     }else{
#       df[i,2]<-NaN
#       df[i,3]<-NaN
#     }
#   }
#   prox_prom_list[[gene]][[k]][[1]]<-df
#   prox_prom_list[[gene]][[k]][[2]]<-c(sym,p1,p2)
#   prox_prom_list[[gene]][[k]][[3]]<-log2(df$p1/df$p2)
#   prox_prom_list[[gene]][[k]][[4]]<-var(prox_prom_list[[gene]][[k]][[3]],na.rm=T)
#   }
# return(prox_prom_list[[gene]])
# }
# toc()
# prox_prom_list<-compact(prox_prom_list)
# prox_prom_list%>%map(~map_dbl(.x,~.x[[4]]))%>%unlist()%>%as_data_frame()%>%filter(!is.na(value))%>%nrow()
# prox_prom_list%>%map(~map_dbl(.x,~.x[[4]]))%>%unlist()%>%as_data_frame()%>%filter(!is.na(value))%>%ggplot(aes(value))+geom_density()+scale_x_log10()

#######
split_proximal<-left_join(proximal_promoters,clusters_selected_genes[,c(1,2)],by=c("prom1"="cluster_id"))%>%distinct()%>%group_split(gene_symbol)
names(split_proximal)<-map_chr(split_proximal,~unique(.x$gene_symbol) ) 

prox_prom_list<-vector(mode="list",length=length(split_proximal)) 
names(prox_prom_list)<-names(split_proximal)

doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

print("Initiate pairwise correlations on the true dataset")
tic()
prox_prom_list<-foreach(gene=1:length(split_proximal))%dopar%{
  prox_prom_list[[gene]]<-vector(mode="list",length=nrow(split_proximal[[gene]]))

for(k in 1:nrow(split_proximal[[gene]])){ 
  df<-data.frame(tissue=unique(tissue), p1=NA,p2=NA) #preallocation
  p1=split_proximal[[gene]][[k,1]]
  p2=split_proximal[[gene]][[k,2]]
  sym<-names(split_proximal)[gene] #find the name of the respective gene symbol
  for(i in 1:length(unique(tissue))){
      df[i,2]<-tissue_beds_ol[[i]][[grep(p1,tissue_beds_ol[[i]]$name)[1],2]]  #counts for one of the promoters in tissue i
      df[i,3]<-tissue_beds_ol[[i]][[grep(p2,tissue_beds_ol[[i]]$name)[1],2]]  #counts for the other one in tissue i
      if(df[i,2]+df[i,3]<20){df[i,2]=df[i,3]=NaN}else if(df[i,2]==0|df[i,3]==0){
        df[i,2]=df[i,2]+0.5
        df[i,3]=df[i,3]+0.5
      }
  }
  prox_prom_list[[gene]][[k]]$df<-df
  prox_prom_list[[gene]][[k]]$gene_symbol<-sym
  prox_prom_list[[gene]][[k]]$p1<-p1
  prox_prom_list[[gene]][[k]]$p2<-p2
  prox_prom_list[[gene]][[k]]$values<-log2(df$p1/df$p2)
  prox_prom_list[[gene]][[k]]$variance<-var(prox_prom_list[[gene]][[k]]$values,na.rm=T)
}
return(prox_prom_list[[gene]])
}
toc()
prox_prom_list<-compact(prox_prom_list)
prox_prom_list%>%map(~map_dbl(.x,~.x$variance))%>%unlist()%>%as_data_frame()%>%filter(!is.na(value))%>%nrow()
prox_prom_list%>%map(~map_dbl(.x,~.x$variance))%>%unlist()%>%as_data_frame()%>%filter(!is.na(value))%>%ggplot(aes(value))+geom_density()+scale_x_log10()
ggsave("variance_density.png")
writeRDS(prox_prom_list,"prox_prom_list.rds")

### PARA PERMUTACIONES. 
#permuted_proximal_promoters<-arrange(proximal_promoters,)[1:length(flatten(prox_prom_list)),c(1,2)]
# permuted_proximal_promoters[,2]<-sample(permuted_proximal_promoters[[2]])
permuted_proximal_promoters<-data.frame(p1=map_chr(flatten(prox_prom_list),~.x$p1),p2=sample(map_chr(flatten(prox_prom_list),~.x$p2)))
prox_prom_perm<-vector(mode="list",length=nrow(permuted_proximal_promoters))
doParallel::registerDoParallel(cl = my.cluster)

print("Initiate pairwise correlations on the permuted dataset")
tic()
prox_prom_perm<-foreach(k=1:length(prox_prom_perm))%dopar%{
#for(k in 1:length(prox_prom_perm)){
  df<-data.frame(tissue=unique(tissue), p1=NA,p2=NA) 
  p1=permuted_proximal_promoters[[k,1]]
  p2=permuted_proximal_promoters[[k,2]]
  sym1<-tissue_beds_ol[[1]][[grep(paste0("^",p1,"$"),tissue_beds_ol[[1]]$name)[1],3]] 
  sym2<-tissue_beds_ol[[1]][[grep(paste0("^",p2,"$"),tissue_beds_ol[[1]]$name)[1],3]]
  for(i in 1:length(unique(tissue))){
    expr_gene1<-tissue_genes[[i]][[grep(paste0("^",sym1,"$"),tissue_genes[[i]]$gene_symbol)[1],2]]
    expr_gene2<-tissue_genes[[i]][[grep(paste0("^",sym2,"$"),tissue_genes[[i]]$gene_symbol)[1],2]]
    #if(is.na(!prox_prom_list[[sym1]][[grep(p1,map_chr(prox_prom_list[[sym1]],~paste(.x$p1,.x$p2)))]]$df[i,2])&is.na(!prox_prom_list[[sym2]][[grep(p2,map_chr(prox_prom_list[[sym2]],~paste(.x$p1,.x$p2)))]]$df[i,2])){
    df[i,2]<-tissue_beds_ol[[i]][[grep(p1,tissue_beds_ol[[i]]$name)[1],2]]  #counts for one of the promoters in tissue i
    df[i,3]<-tissue_beds_ol[[i]][[grep(p2,tissue_beds_ol[[i]]$name)[1],2]]  #counts for the other one in tissue i
    if(df[i,2]+df[i,3]<20){df[i,2]=df[i,3]=NaN}else if(df[i,2]==0|df[i,3]==0){
      df[i,2]=df[i,2]+0.5
      df[i,3]=df[i,3]+0.5
    }
    #}else{df[i,2]=df[i,3]=NaN}
  }
  prox_prom_perm[[k]]$df<-df
  prox_prom_perm[[k]]$gene_symbol1<-sym1
  prox_prom_perm[[k]]$gene_symbol2<-sym2
  prox_prom_perm[[k]]$p1<-p1
  prox_prom_perm[[k]]$p2<-p2
  prox_prom_perm[[k]]$values<-log2((df$p1/expr_gene1)/(df$p2/expr_gene2))
  prox_prom_perm[[k]]$variance<-var(prox_prom_perm[[k]]$values,na.rm=T)
  
  return(prox_prom_perm[[k]])
}
toc()
prox_prom_perm<-compact(prox_prom_perm)
prox_prom_perm%>%map_dbl(~.x$variance)%>%as_data_frame()%>%filter(!is.na(value))%>%nrow()
prox_prom_perm%>%map_dbl(~.x$variance)%>%as_data_frame()%>%filter(!is.na(value))%>%ggplot(aes(value))+geom_density()+scale_x_log10()
ggsave("variance_permuted_density.png")
writeRDS(prox_prom_perm,"prox_prom_perm.rds")

ggplot()+
  geom_density(data=prox_prom_list%>%map(~map_dbl(.x,~.x$variance))%>%unlist()%>%as_data_frame()%>%filter(!is.na(value)),mapping=aes(value))+
  geom_density(data=prox_prom_perm%>%map_dbl(~.x$variance)%>%as_data_frame()%>%filter(!is.na(value)),mapping=aes(value),col="blue")
ggsave("variance_compared_density.png")
