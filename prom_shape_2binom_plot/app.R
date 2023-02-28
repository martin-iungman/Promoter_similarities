invisible(library(shiny))
invisible(library(tidyverse))
shape_entropy<-function(df,norm=F,m=NULL,pseudocount=NULL){
  if(!is.null(pseudocount)){
    df<-tibble(start=seq(min(df$start), max(df$start)), score=pseudocount)%>%filter(!start%in%df$start)%>%rbind(df%>%select(start, score))
  }
  tot_sum<-sum(df$score)
  df$relscore<-df$score/tot_sum
  shape_entropy=2+sum(df$relscore*log2(df$relscore))
  if(norm){
    if(is.null(m)){
      m=abs(max(df$start)-min(df$start))}
    shape_entropy=(-1/log(m))*sum(df$relscore*log(df$relscore))
  }
  return(shape_entropy)
}

shape_extropy<-function(df,include_zeros=F,norm=F,m=NULL){
  if(include_zeros){
    df<-tibble(start=seq(min(df$start), max(df$start)), score=0)%>%filter(!start%in%df$start)%>%rbind(df%>%select(start, score))
  }
  tot_sum<-sum(df$score)
  df$relscore<-df$score/tot_sum
  shape_extropy=sum(df$relscore*log2(df$relscore))
  if(norm){
    if(is.null(m)){
      m=abs(max(df$start)-min(df$start))}
    shape_extropy=(-1/((m-1)*log(m/(m-1))))*sum(1-df$relscore*log(1-df$relscore))
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
  sample.vctr<-map(seq(n.samples),~sample(pop.vctr,counts.sample))
  sample.df<-list()
  sample.stats<-list()
  for(i in 1:n.samples){
    df<-data.frame(table(sample.vctr[[i]]))
    names(df)<-c("start","score")
    df$start<-df$start%>%as.character()%>%as.integer()
    sample.df[[i]]<-df
    sample.stats[[i]]<-data.frame(entropy=shape_entropy(df),norm_entropy=shape_entropy(df,norm=T),
                                  extropy=shape_extropy(df),norm_extropy=shape_extropy(df,norm=T),
                                  gini=calculate_gini(df))
  }  
  sample.stats<-sample.stats%>%list_rbind()
  return(sample.stats)
}

accuracy<-function(observed,real){
  return(median(abs((observed-real)/real)))
}

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Promoter Shape Metrics (double binomial)"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      sliderInput("counts",
                  "Promoter counts: ",
                  min = 1,
                  max = 1000,
                  value = 50),
      sliderInput("count.ratio",
                  "Ratio of counts between the two curves",
                  min=0,
                  max=1,
                  value=0.5),
      
      sliderInput("len1",
                  "Binomial size (1): ",
                  min = 1,
                  max = 100,
                  value = 30),
      sliderInput("len2",
                  "Binomial size (2): ",
                  min = 1,
                  max = 100,
                  value = 30),
      sliderInput("p1",
                  "Probability shape (1): ",
                  min = 0,
                  max = 1,
                  value = 0.5),
      sliderInput("p2",
                  "Probability shape (2): ",
                  min = 0,
                  max = 1,
                  value = 0.5),

      sliderInput("perc.sample",
                  "Percentage of counts sampled: ",
                  min = 1,
                  max=100,
                  value=10),
      sliderInput("n.samples",
                  "Number of samples: ",
                  min = 1,
                  max=1000,
                  value=100)
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot"),
      tableOutput("tableSumm"),
      tableOutput("tableSampleP"),
      tableOutput("tableSampleA"),
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  pop.vctr<-reactive({ 
    pop.vctr=c(rbinom(floor(input$counts*input$count.ratio), input$len1, input$p1),
               rbinom(floor(input$counts*(1-input$count.ratio)), input$len2, input$p2))+1  
    return(pop.vctr)
  })
  
  pop.df<-reactive({
    pop.df<-data.frame(table(pop.vctr()))
    names(pop.df)<-c("start","score")
    pop.df$start<-pop.df$start%>%as.character()%>%as.integer()
    return(pop.df)
  })
  
  stats<-reactive({
    pop.df=pop.df()
    pop.stats<-data.frame(entropy=shape_entropy(pop.df),norm_entropy=shape_entropy(pop.df,norm=T),
                          extropy=shape_extropy(pop.df),norm_extropy=shape_extropy(pop.df,norm=T),
                          gini=calculate_gini(pop.df), range=max(pop.df$start)-min(pop.df$start))
    return(pop.stats) 
  })
  
  samples<-reactive({
    counts.sample=floor(input$perc.sample*input$counts/100)
    return(get_stats_samples(pop.vctr(),input$n.samples,counts.sample))
  })
  
  precision_samples<-reactive({
    samples<-samples()
    summ_table<-data.frame(CV_entropy=abs(mad(samples$entropy)/median(samples$entropy)),
                           CV_norm_entropy=abs(mad(samples$norm_entropy)/median(samples$norm_entropy)),
                           CV_extropy=abs(mad(samples$extropy)/median(samples$extropy)),
                           CV_norm_extropy=abs(mad(samples$norm_extropy)/median(samples$norm_extropy)),
                           CV_gini=abs(mad(samples$gini)/median(samples$gini)))
    return(summ_table)
  })
  
  accuracy_samples<-reactive({
    samples<-samples()
    pop<-stats()
    summ_table<-data.frame(acc_entropy=accuracy(samples$entropy,pop$entropy),
                           acc_norm_entropy=accuracy(samples$norm_entropy,pop$norm_entropy),
                           acc_extropy=accuracy(samples$extropy,pop$extropy),
                           acc_norm_extropy=accuracy(samples$norm_extropy,pop$norm_extropy),
                           acc_gini=accuracy(samples$gini,pop$gini))
    
    return(summ_table)
  })
  
  output$distPlot <- renderPlot({
    pop.df()%>%ggplot(aes(start,score))+geom_col()+theme_bw()+xlab("Position")
  })
  output$tableSumm<-renderTable(stats())
  output$tableSampleP<-renderTable(precision_samples())
  output$tableSampleA<-renderTable(accuracy_samples())
  
}

# Run the application 
shinyApp(ui = ui, server = server)
