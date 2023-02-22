library(tidyverse)
shape_entropy<-function(df,norm=F,m=NULL){
  tot_sum<-sum(df$score)
  df$relscore<-df$score/tot_sum
  shape_entropy=2+sum(df$relscore*log2(df$relscore))
  if(norm){
    if(is.null(m)){
      m=abs(max(df$start)-min(df$start))+1}
    shape_entropy=(-1/log(m))*sum(df$relscore*log(df$relscore))
  }
  return(shape_entropy)
}
shape_extropy<-function(df,norm=F,m=NULL){
  tot_sum<-sum(df$score)
  df$relscore<-df$score/tot_sum
  shape_extropy=sum(df$relscore*log2(df$relscore))
  if(norm){
    if(is.null(m)){
      m=abs(max(df$start)-min(df$start))+1}
    shape_extropy=(-1/((m-1)*log(m/(m-1))))*sum((1-df$relscore)*log(1-df$relscore))
  }
  return(shape_extropy)
}
gini_coefficient <- function(x, weights=rep(1,length=length(x))){ ## copied form reldist package
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox]/sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights*x)
  n <- length(nu)
  nu <- nu / nu[n]
  sum(nu[-1]*p[-n]) - sum(nu[-n]*p[-1])
}
calculate_gini<-function(df, include_zeros=F, alpha=NULL){
  if(!is.null(alpha)){
    width_ic<-calculate_width(df,alpha)
    df<-df%>%filter(start>=width_ic[["low"]],start<=width_ic[["up"]])
  }
  if(include_zeros){
    df<-tibble(start=seq(min(df$start), max(df$start)), score=0)%>%filter(!start%in%df$start)%>%rbind(df%>%select(start, score))
  }
  gini_coefficient(df$score)
}
get_stats_samples<-function(pop.vctr,n.samples,counts.sample){
  sample.vctr<-map(1:n.samples,~sample(pop.vctr,counts.sample))
  sample.df<-list()
  sample.stats<-list()
  for(i in 1:n.samples){
    df<-data.frame(table(sample.vctr[[i]]))
    names(df)<-c("start","score")
    df$start<-df$start%>%as.integer()
    sample.df[[i]]<-df
    sample.stats[[i]]<-data.frame(entropy=shape_entropy(df),norm_entropy=shape_entropy(df,norm=T),
                                  extropy=shape_extropy(df),norm_extropy=shape_extropy(df,norm=T),
                                  gini=calculate_gini(df))
  }
  sample.stats<-sample.stats%>%list_rbind()
  return(sample.stats)}

# max.len=30
# min.counts=200
# max.counts=1000
#len=sample.int(30)[1]
#counts=sample.int(1000)
#counts=counts[counts%inseq(min.counts,max.counts)]
#pop.vctr=sample.int(len, size=counts, replace=T, prob = min(c(pnorm(1:len, mean=len/2),pnorm(1:len, mean=len/2, lower.tail = F))))

#norm distribution
n.samples=100
counts_vctr=seq(100,1000,by=100)
len_vctr=seq(5,50,by=10)
p_vctr=seq(0.1,0.9,by=0.1)
params<-expand.grid(counts=counts_vctr, len=len_vctr,p=p_vctr)
perc.sample=0.1
pop.vctr=list()
pop.stats=tibble(pop.entropy=NULL,pop.norm_entropy=NULL,pop.extropy=NULL,pop.norm_extropy=NULL,pop.gini=NULL,pop.range=NULL,counts=NULL,len=NULL,p=NULL)
sample.stats=list()
for(j in 1:nrow(params)){
  if(j%%5==0){print(j/nrow(params))}
  counts=params[j,1]
  len=params[j,2]
  p=params[j,3]
pop.vctr[[j]]=rbinom(counts, len, p=p)+1 #el +1 es porque el 0 esta incluido en la distribucion binomial
pop.range=max(pop.vctr[[j]])-min(pop.vctr[[j]])
# hist(pop.vctr[[j]],breaks=seq(min(pop.vctr[[j]]),max(pop.vctr[[j]])))
pop.df<-data.frame(table(pop.vctr[[j]]))
names(pop.df)<-c("start","score")
pop.df$start<-pop.df$start%>%as.integer()
pop.stats<-bind_rows(pop.stats,c(pop.entropy=shape_entropy(pop.df),pop.norm_entropy=shape_entropy(pop.df,norm=T),
                      pop.extropy=shape_extropy(pop.df),pop.norm_extropy=shape_extropy(pop.df,norm=T),
                      pop.gini=calculate_gini(pop.df), pop.range=pop.range, counts=counts, len=len, p=p))

counts.sample=perc.sample*counts
sample.stats[[j]]<-get_stats_samples(pop.vctr[[j]],n.samples,counts.sample)%>%mutate(n.samples=n.samples,perc.sample=perc.sample, counts.sample=counts.sample)%>%bind_cols(pop.stats[j,])
}
# sample.stats[[j]]%>%ggplot()+geom_density(aes(entropy,fill="entropy"))+geom_density(aes(norm_entropy,,fill="norm_entropy"))+
#   geom_density(aes(extropy,fill="extropy"))+geom_density(aes(norm_extropy, fill="norm_extropy"))+geom_density(aes(gini, fill="gini"))+xlab("Measure")
sample.stats.sum<-map(sample.stats,~summarise(.x,counts=unique(counts),len=unique(len),p=unique(p),
                                              exact_entropy=unique(abs((median(entropy)-pop.entropy)/pop.entropy)),
                                              exact_norm_entropy=unique(abs((median(norm_entropy)-pop.norm_entropy)/pop.norm_entropy)),
                                              exact_extropy=unique(abs((median(extropy)-pop.extropy)/pop.extropy)),
                                              exact_norm_extropy=unique(abs((median(norm_extropy)-pop.norm_extropy)/pop.norm_extropy)),
                                              exact_gini=unique(abs((median(gini)-pop.gini)/pop.gini)),
                                              across(c("entropy","norm_entropy","extropy","norm_extropy","gini"),.fns=~abs((sd(.x)/mean(.x))))))%>%list_rbind()


hist(c(rbinom(100, 20, 0.4),rbinom(30,50,0.5)))

