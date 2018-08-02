# Required Packages and Functions -----------------------------------------
list.of.packages <- c("shiny", "magrittr","ggplot2", "fields", "data.table", "stringr", "ape", "DT", 
                      "shinyjs", "cosinor2")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
}
library(shiny)
library(magrittr)
library(ggplot2)
library(fields)
library(data.table)
library(stringr)
library(ape)
library(DT)
library(shinyjs)
library(cosinor2)
source("functions.R")
###=== end of packages and functions loading ===###
ui<-fluidPage(
  useShinyjs(),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        'input.dataset === "Cosinor Analysis of Rhythms"',
        h4("Cosinor Analysis of Rhythms"),
        verticalLayout(
          fileInput(
            'raw',
            label = h4(popover(title = "?", 
                               content = a("Please click here to see an example",
                                           href="./Example/results.csv",
                                           download="results.csv",
                                           target="_blank"), 
                               trigger = "focus", html = TRUE),
                       "Upload your wave signal in .csv file:"),
            accept = c('csv', 'comma-separated-values', '.csv')
          ),
          checkboxInput('style', 'Scatter plot'),
          numericInput('TimeTagsCol', 'Cloumn number of time', 1),
          numericInput('ValueCols', 'Start cloumn number of sample', 2),
          numericInput('ValueCole', 'End Cloumn number of sample', 4),
          numericInput('period', 'Period assumption', 24),
          textInput('xtitle', 'Label of x-axis', "Time(h)"),
          textInput('ytitle', 'Label of y-axis', "Expression Level"),
          actionButton("OK", "OK", class="btn-primary")
        )
      ),
      conditionalPanel(
        'input.dataset === "About & Help"',
        h4("About & Help")
      )
    ),
    mainPanel(
      tabsetPanel(
        id='dataset',
        tabPanel(
          "Cosinor Analysis of Rhythms",
          verticalLayout(
            plotOutput("plot.out"),
            textOutput("F.statistic"),
            textOutput("rhythm.p"),
            htmlOutput("CurveFun"),
            textOutput("R2"),
            textOutput("P.value")
          )
        ),
        tabPanel(
          "About & Help"
        )
      )
    )
  )
)
server<-function(input, output, session) {
  res<-eventReactive(input$OK, {
    raw<-read.csv(input$raw$datapath, header = TRUE, sep = ",")
    #Parameters setting
    TimeTagsCol<-input$TimeTagsCol
    period<-input$period
    xtitle<-input$xtitle
    ytitle<-input$ytitle
    ValueCol<-c(input$ValueCols:input$ValueCole)
    temp<-data.frame(Time=raw[, TimeTagsCol], Value=rowMeans(raw[, ValueCol]), SD=rowSds(as.matrix(raw[, ValueCol])))
    fit_per_est<-cosinor.lm(Value~time(Time), period = period, data = temp)
    FitCurve<-data.frame(x=temp$Time, y=fit_per_est$coefficients[1]+fit_per_est$coefficients[2]*cos(2*pi*temp$Time/period+2*fit_per_est$coefficients[3]))
    FitCurve<-data.frame(x=temp$Time, y=fit_per_est$fit$fitted.values)
    res<-cosinor.detect(fit_per_est)
    ForScatter<-melt(raw, id.vars=TimeTagsCol)
    #Results Part
    CurveFun<-paste0(round(fit_per_est$coefficients[1],2),
                     " + ", round(fit_per_est$coefficients[2],2),
                     "cos(2\u03C0t/", period, " + ", round(pi-fit_per_est$coefficients[3],2),
                     ")")
    F.statistic<-res[1]
    rhythm.p<-res[4]
    R2<-cor(temp$Value, FitCurve$y)^2
    P.value<-cor.test(temp$Value, FitCurve$y)$p.value
    res<-list(raw, temp, FitCurve, ForScatter, CurveFun, F.statistic, rhythm.p, R2, P.value)
    return(res)
  })
  plot.out<-reactive({
    if (is.null(res())){
      return(NULL)
    } else {
      ForScatter<-res()[4][[1]]
      FitCurve<-res()[3][[1]]
      temp<-res()[2][[1]]
      xtitle<-input$xtitle
      ytitle<-input$ytitle
      if (input$style == TRUE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=value), ForScatter)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none")+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)
      } else if (input$style == FALSE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=Value), temp)+
          geom_errorbar(aes(x=Time, ymax=Value+SD, ymin=Value-SD), temp)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none")+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)
      }
      return(gp)
    }
  })
  output$plot.out<-renderPlot(plot.out())
  output$F.statistic<-renderText(paste0("F statitics: ", res()[6]))
  output$rhythm.p<-renderText(paste0("P-value for rhythms: ", res()[7]))
  output$CurveFun<-renderUI(HTML(paste0("Function of fitted curve: ", res()[5])))
  output$R2<-renderText(paste0("R squared (Goodness of fit): ", res()[8]))
  output$P.value<-renderText(paste0("P-value of fit: ", res()[9]))
  session$onSessionEnded(stopApp)
}
shinyApp(ui=ui, server=server)