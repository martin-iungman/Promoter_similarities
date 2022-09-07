library(rtracklayer)
library(stringi)
library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38<-BSgenome.Hsapiens.UCSC.hg38

epd_all_prom<-import.bed("Data/All_promoters.bed")
mcols(epd_all_prom)$gene<-mcols(epd_all_prom)$name%>%str_remove("_.+")
df_genes<-mcols(epd_all_prom)%>%as_tibble()%>%group_by(gene)%>%mutate(n_prom=n())

#histograma de numero de prom por gen
df_genes%>%select(-name)%>%distinct()%>%ggplot(aes(n_prom))+geom_histogram(stat="count")
#numeor de genes con mas de un promotor:
df_genes%>%filter(n_prom>1)%>%select(gene)%>%distinct()%>%nrow()

epd_list_genes<-split(epd_all_prom, epd_all_prom$gene)
dist_to_near_intra<-map(as.list(epd_list_genes),distanceToNearest)
dist_epd<-data.frame(
  cluster_id=unlist(map2(as.list(epd_list_genes),dist_to_near_intra,~.x[queryHits(.y)]$name)),
  nearest_id=unlist(map2(as.list(epd_list_genes),dist_to_near_intra,~.x[subjectHits(.y)]$name)),
  distances=unlist(map(dist_to_near_intra, ~as_vector(mcols(.x))))
  )

df_genes<-left_join(df_genes, dist_epd, by=c("name"="cluster_id"))
#density de distancia entre promotores del mismo gen
df_genes%>%filter(n_prom>1)%>%select(gene,distances)%>%distinct()%>%ggplot( aes(distances))+geom_density()+scale_x_log10()

df_genes%>%filter(n_prom>1)%>%group_by(gene)%>%summarise(mean_distance=mean(distances))%>%
  ggplot( aes(mean_distance))+geom_density()+scale_x_log10()+geom_vline(xintercept = 100, linetype="dashed")

#Ahora armo promotores para buscar motivos core
epd_down100_prom<-epd_all_prom%>%resize(100, fix="end")

#Busco presencia de TATA
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38<-BSgenome.Hsapiens.UCSC.hg38

pfm_tbp <- getMatrixByName(JASPAR2020, "TBP")
seqLogo(toICM(pfm_tbp))
pwm_tbp<-toPWM(pfm_tbp, pseudocounts=0.8)
dna<-getSeq(Hsapiens,epd_down100_prom)

siteset<- searchSeq(pwm_tbp,dna, min.score="80%", strand="*")
gff3_tata<-writeGFF3(siteset)
tata_prom<-epd_down100_prom[as.numeric(unique(gff3_tata$seqname))]$name

df_genes<-cbind(df_genes, strand=as.vector(strand(epd_down100_prom)))
df_genes<-df_genes%>%mutate(TATA=ifelse(name%in%tata_prom, T,F))
df_genes_tata<-df_genes%>%filter(name%in%tata_prom)%>%cbind(pos=unique(gff3_tata$seqname))
df_genes_tata%>%group_by(gene)%>%summarise(n_prom_tata=n())
gff3_tata_ext<-gff3_tata%>%rename(strand_TATA=strand)%>%left_join(df_genes_tata, by=c("seqname"="pos"))%>%filter(strand_TATA==strand)%>%group_by(name)%>%filter(score.x==max(score.x))
gff3_tata_ext<-gff3_tata_ext%>%mutate(start=100-start,end=100-end)
  
df_genes<-df_genes%>%group_by(gene)%>%mutate(prop_TATA=round(sum(TATA)/n_prom,digits=3))
#Density de promotores con TATA
df_genes%>%filter( n_prom>1)%>%ggplot(aes(prop_TATA))+geom_density(adjust=2)

gff3_tata_ext%>%ggplot(aes(start))+geom_density()+xlim(0,100)+geom_vline(xintercept = 31,linetype="dashed")

# Voy a comparar con los genes TAT& que me da EPD 
epd_tata<-import.bed("Data/TATA.bed")
length(epd_tata)
sum(epd_tata$name%in%gff3_tata_ext$name)/length(epd_tata)

#relacion numero de promotores con presencia de TATA en uno. 
#Es un poco engañoso dado que un gen con mas cantidad de promo0tores es mas probable que tenga alguno que sea TAT (categoria rara)
df_genes%>%ggplot(aes(n_prom, fill=prop_TATA>0))+geom_histogram(adjust=2,alpha=0.5, position = "dodge")+scale_x_log10()
#hago test
df_genes<-df_genes%>%mutate(mult_prom=ifelse(n_prom>1,TRUE,FALSE))
table(df_genes$mult_prom,df_genes$TATA)
fisher.test(table(df_genes$mult_prom,df_genes$TATA))
library(ggmosaic)
df_genes%>%ggplot()+geom_mosaic(aes(x=product(TATA),fill=mult_prom))

#esperados vs observados para distribucion de TATA en genes con dos promotores
n_tata2<-df_genes%>%ungroup()%>%filter(n_prom==2)%>%nrow()
frec_tata2<-df_genes%>%ungroup()%>%filter(n_prom==2)%>%select(TATA)%>%as_vector()%>%sum(.)
sin_tata2_esp<-((n_tata2-frec_tata2)/n_tata2)^2
sin_tata2_obs<-df_genes%>%ungroup()%>%filter(n_prom==2,prop_TATA==0)%>%nrow()/n_tata2
con_tata2_esp<-(frec_tata2/n_tata2)^2
con_tata2_obs<-df_genes%>%ungroup()%>%filter(n_prom==2,prop_TATA==1)%>%nrow()/n_tata2
het_tata2_esp<-2*(frec_tata2/n_tata2)*  ((n_tata2-frec_tata2)/n_tata2)          
het_tata2_obs<-df_genes%>%ungroup()%>%filter(n_prom==2,prop_TATA==0.5)%>%nrow()/n_tata2

# Aplico Hardy Weinberg
library(HardyWeinberg)
HWChisq(c(AA=df_genes%>%ungroup()%>%filter(n_prom==2,prop_TATA==0)%>%nrow()/2,
          AB=df_genes%>%ungroup()%>%filter(n_prom==2,prop_TATA==0.5)%>%nrow()/2,
          BB=df_genes%>%ungroup()%>%filter(n_prom==2,prop_TATA==1)%>%nrow()/2))

# La D=0.5(obs_het -esp_het). Osea que D<0 es que hay menos "heterocigotas" de lo esperado.

# Aplico modelo de bondad de ajuste a todas las distribuciones hasta 6 promotores alternativos. 
# Uso Binomio de Newton para determinar la distribucion de los esperados


TATA_modelo<-data.frame()
for (n in 1:6) { #numero de promotores alternativos del gen
  N=df_genes%>%filter(n_prom==n)%>%nrow() #numero de promotores totales en la subpoblacion
  p=df_genes%>%filter(n_prom==n,TATA==T)%>%nrow()/N #frecuencia de promotores TATA en la subpoblacion
  q=df_genes%>%filter(n_prom==n,TATA==F)%>%nrow()/N ##frecuencia de promotores sin TATA en la subpoblacion
  for(k in 0:n){ #numero de promotores TATA de los n que hay en los genes
    vector<-c(n,N,N/n,k=k,
              esp=choose(n,k)*p^(n-k)*q^k, #binomio de newton
              obs=df_genes%>%filter(n_prom==n,prop_TATA==round(1-k/n,3))%>%nrow()/N,
              obs_raw=df_genes%>%filter(n_prom==n,prop_TATA==round(1-k/n,3))%>%nrow()/n,
              prop=k/n)
    TATA_modelo<-rbind(TATA_modelo,vector)
  }
}
names(TATA_modelo)<-c("n_prom","N_prom","N","k","esp","obs","obs_raw","prop_TATA")
chisq_list<-map(1:6, ~chisq.test(
  TATA_modelo%>%filter(n_prom==.x)%>%select(obs_raw)%>%as_vector(), 
  p=TATA_modelo%>%filter(n_prom==.x)%>%select(esp)%>%as_vector()))
