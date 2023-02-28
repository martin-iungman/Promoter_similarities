invisible(library(shiny))
invisible(library(tidyverse))
invisible(library(data.table))
invisible(library(DT))
library(rtracklayer)
files<-list.files("../Tables/CAGE_filtered_bed/",full.names = T)
samples<-sub("_filt.bed","",files)%>%sub("^.+bed//","",.)
calculate_width<-function(df,alpha=0){
  df<-arrange(df, start)
  df<-data.frame(cumsum=cumsum(df$score),pos=df$start)
  if(alpha==0){
    up=max(df$start)
    low=min(df$start)
    width=abs(up-low)
  }else{ 
  suppressWarnings(low<-df$pos[min(which(df$cumsum<(alpha/2)*max(df$cumsum)))+1])
  if(is.na(low)){low=min(df$pos)}
  up<-df$pos[min(which(df$cumsum>(1-alpha/2)*max(df$cumsum)))]
  width<-abs(up-low)+1}
  return(list(width=width,low=low,up=up))
}
plot_width<-function(df, alpha, type=c("histogram", "cumsum"), name=NULL){
  df<-arrange(df, start)
  df_sum<-data.frame(cumsum=cumsum(df$score),pos=df$start)
  suppressWarnings(low<-df_sum$pos[min(which(df_sum$cumsum<(alpha/2)*max(df_sum$cumsum)))+1])
  if(is.na(low)){low=min(df_sum$pos)}
  up<-df_sum$pos[min(which(df_sum$cumsum>(1-alpha/2)*max(df_sum$cumsum)))]
  if(type=="histogram"){
    plot<-ggplot(df, aes(start,score))+geom_col()+geom_vline(xintercept = c(up+0.5,low-0.5),linetype="dashed")+xlab("Position")+ylab("Counts")
  }
  if(type=="cumsum"){
    plot<-df_sum%>%ggplot(aes(pos,cumsum))+geom_line()+geom_point()+geom_vline(xintercept = c(up,low),linetype="dashed")+xlab("Position")+ylab("Cumulated sum")
  }
  if(!is.null(name)){
    plot<-plot+ggtitle(name)
  }
  print(plot)
}
shape_entropy<-function(df,alpha=0,norm=F,m=NULL,pseudocount=NULL){
    if(!is.null(pseudocount)){
    df<-tibble(start=seq(min(df$start), max(df$start)), score=pseudocount)%>%filter(!start%in%df$start)%>%rbind(df%>%select(start, score))
    }
  if(alpha!=0){
    width_ic<-calculate_width(df,alpha)
    df<-df%>%filter(start>=width_ic[["low"]],start<=width_ic[["up"]])
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

shape_extropy<-function(df,alpha=0,include_zeros=F,norm=F,m=NULL){
  if(include_zeros){
    df<-tibble(start=seq(min(df$start), max(df$start)), score=0)%>%filter(!start%in%df$start)%>%rbind(df%>%select(start, score))
  }
  if(alpha!=0){
    width_ic<-calculate_width(df,alpha)
    df<-df%>%filter(start>=width_ic[["low"]],start<=width_ic[["up"]])
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
calculate_gini<-function(df, include_zeros=F, alpha=0){
  if(alpha!=0){
    width_ic<-calculate_width(df,alpha)
    df<-df%>%filter(start>=width_ic[["low"]],start<=width_ic[["up"]])
  }
  if(include_zeros){
    df<-tibble(start=seq(min(df$start), max(df$start)), score=0)%>%filter(!start%in%df$start)%>%rbind(df%>%select(start, score))
  }
  gini_coefficient(df$score)
}
ui <- fluidPage(

    titlePanel("Promoter Shape"),

    sidebarLayout(
        sidebarPanel(
          selectizeInput("tissue","Select a CAGE sample: ",choices=samples),
          numericInput("min.activ", "Minimum number of counts per promoter:",10),
          selectizeInput("prom", "Select a promoter: ", choices=""),
          sliderInput("alpha","Alpha", 0,1,0.1, step=0.01),
          checkboxInput("complete_table","Show complete table"),
          checkboxInput("apply_alpha", "Restrict promoters length prior to statistics?"),
          actionButton("refresh_table","Refresh Table")
        ),

        mainPanel(
          tabsetPanel(
            tabPanel("Plots",
                     plotOutput("hist"),
                     plotOutput("cumsum")),
            tabPanel("Table", tags$br(),tags$br(), dataTableOutput("stats"))
          )
        )
    )
)

server <- function(input, output,session) {
  import_bed<-reactive({
    import.bed(files[which(samples==input$tissue)])})
  active_prom<-reactive({
    import_bed()%>%as_tibble()%>%group_by(name)%>%summarise(sum=sum(score))%>%filter(sum>=input$min.activ)%>%select(name)%>%as_vector()
  })  
  observeEvent(req(input$tissue, input$min.activ),{
    updateSelectizeInput(session,inputId="prom",choices=unname(active_prom()),server=T)})
  list_proms<-reactive({
    bed_filt<-subset(import_bed(), name%in%active_prom())
    list_prom<-group_split(as_tibble(bed_filt), name)
    names(list_prom)<-group_keys(as_tibble(bed_filt)%>%group_by(name))%>%as_vector()
    return(list_prom)
  })
  val<-reactive({
    if(input$apply_alpha){
      input$alpha
    }else{0}
  })
  
  plot_hist<-reactive({
    req(input$prom)
    plot_width(list_proms()[[input$prom]], alpha=input$alpha, "histogram", name=input$prom)})
  plot_cumsum<-reactive({
    req(input$prom)
    plot_width(list_proms()[[input$prom]], alpha=input$alpha, "cumsum", name=input$prom)})
  observeEvent(input$refresh_table,{
    stats<-list()
    if(input$complete_table){
      all_prom=active_prom()
    }else{all_prom=input$prom}
    for (prom in all_prom) {
      df<-list_proms()[[prom]]
      stats[[which(active_prom()==prom)]]<-data.frame(name=prom, counts=sum(df$score),width=calculate_width(df,val())[["width"]],entropy=shape_entropy(df,alpha=val()),norm_entropy=shape_entropy(df,alpha=val(), norm=T),
                                                          extropy=shape_extropy(df,alpha=val()),norm_extropy=shape_extropy(df,alpha=val(),norm=T),
                                                          gini=calculate_gini(df,alpha=val()), gini_zero=calculate_gini(df,alpha=val(),include_zero=T))%>%mutate(across(where(is.numeric),round,3))
      
    }
    stats%>%list_rbind()
    output$stats<-DT::renderDataTable(data.table::as.data.table(stats))
  })
  output$hist<-renderPlot(plot_hist())
  output$cumsum<-renderPlot(plot_cumsum())

}

shinyApp(ui = ui, server = server)
