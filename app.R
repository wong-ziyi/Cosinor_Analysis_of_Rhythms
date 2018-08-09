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
# Shiny UI part -----------------------------------------
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
          numericInput('ValueCols', 'Start cloumn number of sample', 2, min = 1, step = 1),
          numericInput('ValueCole', 'End Cloumn number of sample', 4, min = 1, step = 1),
          numericInput('XInterval', 'X-axis Display Interval', 1, min=1),
          textInput('xtitle', 'Label of x-axis', "Time(h)"),
          textInput('ytitle', 'Label of y-axis', "Expression Level"),
          splitLayout(
            withBusyIndicatorUI(
              actionButton("OK", "OK", class="btn-primary")
            ),
            downloadButton('Results.csv', 'Result Download')
          )
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
            textOutput("P.value"),
            plotOutput("plot.period")
          )
        ),
        tabPanel(
          "About & Help"
        )
      )
    )
  )
)
# Shiny server function -----------------------------------------
server<-function(input, output, session) {
  res<-eventReactive(input$OK, {
    withBusyIndicatorServer("OK", {
      if (is.null(input$raw$datapath)) {
        stop("Please submit your data first")
        return(NULL)
      } else {
        raw<-read.csv(input$raw$datapath, header = TRUE, sep = ",")
        #Parameters setting (Passed from UI input)
        TimeTagsCol<-input$TimeTagsCol #Column number of the time series
        ValueCol<-c(input$ValueCols:input$ValueCole) #Column number of the results
        xtitle<-input$xtitle #Title of x-axis
        ytitle<-input$ytitle #Title of y-axis
        #Set up x-axis label of ticks by iterative method
        TimeSeq<-c()
        TimeS<-raw[, TimeTagsCol]
        for (i in 1:length(TimeS)) {
          if(TimeS[i]<24){
            TimeSeq0<-TimeS[i]
          } else if (TimeS[i]>=24){
            TimeSeq0<-TimeS[i]-24
            TimeS<-TimeS-24
          }
          TimeSeq<-c(TimeSeq, TimeSeq0)
        }
        ticks<-TimeSeq # Make final sequence for labling the x-axis ticks
        #++///////////////////////////////////////////////////////////////////////////++#
        #++==In follow processing, only mean value was used to fit a cosine function==++#
        #++///////////////////////////////////////////////////////////////////////////++#
        #Make temporal data that contains the calculated mean and standrad deviation
        if(input$ValueCols != input$ValueCole){
          temp<-data.frame(Time=raw[, TimeTagsCol], Value=rowMeans(raw[, ValueCol]), SD=rowSds(as.matrix(raw[, ValueCol])))
        } else if(input$ValueCols == input$ValueCole){
          temp<-data.frame(Time=raw[, TimeTagsCol], Value=raw[, ValueCol], SD=0)
        }
        #Estimate the period by modified iterative function from cosinor2
        period0<-periodogram_wzy(data = temp, timecol = 1, firstsubj = 2, lastsubj = 2)
        period<-period0$plot_env$best # Pass the best results
        #Get best fitted cosinor model
        temp0<-temp
        temp0[,1]<-temp0[,1]*60 #convert hours to minutes
        fit_per_est<-cosinor.lm(Value~time(Time), period = period, data = temp0)
        period<-period/60 #convert minutes to hours
        #Calculate the estimated fitted value
        FitCurve<-data.frame(x=seq(from=first(raw[, TimeTagsCol]), to=last(raw[, TimeTagsCol]), by=0.5), y=fit_per_est$coefficients[1]+fit_per_est$coefficients[2]*cos(2*pi*seq(from=first(raw[, TimeTagsCol]), to=last(raw[, TimeTagsCol]), by=0.5)/period+pi-fit_per_est$coefficients[3]))
        #Statistically detect the rhythms
        res<-cosinor.detect(fit_per_est)
        #Make temporal for scatter plot
        ForScatter<-melt(raw[, c(TimeTagsCol, ValueCol)], id.vars=TimeTagsCol)
        colnames(ForScatter)<-c("Time", "variable","value")
        ForScatter<-cbind(ForScatter, test=fit_per_est$coefficients[1]+fit_per_est$coefficients[2]*cos(2*pi*ForScatter$Time/period+pi-fit_per_est$coefficients[3]))
        #Results Part (combine all results and output)
        CurveFun<-paste0(format(round(fit_per_est$coefficients[1],3), nsmall=3),
                         " + ", format(round(fit_per_est$coefficients[2],3), nsmall=3),
                         "cos(2\u03C0t/", format(round(period, 3), nsmall=3), " + ", format(round(pi-fit_per_est$coefficients[3],3), nsmall=3),
                         ")") # cosinor model
        F.statistic<-res[1] # F statistics
        rhythm.p<-res[4] # P-value for F statistics
        R2<-cor(ForScatter$value, ForScatter$test)^2 # Coefficient of Determination (GOF, goodness of fit) which be calculated from Perason's correlation coefficient
        P.value<-cor.test(ForScatter$value, ForScatter$test)$p.value # Significance for this Perason's correlation coefficient
        res.out<-list(raw, temp, FitCurve, ForScatter, CurveFun, F.statistic, rhythm.p, R2, P.value, ticks, fit_per_est, period, period0) # Built output results
      } #UI effect: error indicator. end.
    }) #UI effect: busy indicator. end. 
    return(res.out)
  })
  plot.out<-reactive({
    if (is.null(res())){
      return(NULL)
    } else {
      #Pass value to local
      ticks<-res()[10][[1]]
      ForScatter<-res()[4][[1]]
      FitCurve<-res()[3][[1]]
      temp<-res()[2][[1]]
      xtitle<-input$xtitle
      ytitle<-input$ytitle
      #Scatter or Mean SD?
      if (input$style == TRUE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=value), ForScatter)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=input$XInterval)], labels=ticks[seq(from=1, to=length(temp$Time), by=input$XInterval)])
      } else if (input$style == FALSE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=Value), temp)+
          geom_errorbar(aes(x=Time, ymax=Value+SD, ymin=Value-SD), temp)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=input$XInterval)], labels=ticks[seq(from=1, to=length(temp$Time), by=input$XInterval)])
      }
      return(gp)
    }
  })
  output$plot.out<-renderPlot(plot.out())
  output$plot.period<-renderPlot(print(res()[13][[1]]))
  output$F.statistic<-renderText(paste0("F statitics: ", res()[6]))
  output$rhythm.p<-renderText(paste0("P-value for rhythms: ", res()[7]))
  output$CurveFun<-renderUI(HTML(paste0("Function of fitted curve: ", res()[5])))
  output$R2<-renderText(paste0("R squared (Goodness of fit): ", res()[8]))
  output$P.value<-renderText(paste0("P-value of fit: ", res()[9]))
  output$Results.csv<-downloadHandler(
    filename = function(){
      paste("Results", "csv", sep = ".")
    },
    content = function(fname){
      res.out<-data.frame(terms=c("MESOR:", 
                                  "Amplitude:", 
                                  "Period:", 
                                  "Phase:", 
                                  "F statistics:", 
                                  "P-value for rhythm detection:", 
                                  "Fitted Curve Function:", 
                                  "R squared (Goodness of fit):", 
                                  "P-value for fit"), 
                          results=c(res()[11][[1]]$coefficients[1],
                                    res()[11][[1]]$coefficients[2],
                                    res()[12][[1]],
                                    res()[11][[1]]$coefficients[3],
                                    res()[6][[1]],
                                    res()[7][[1]],
                                    res()[5][[1]],
                                    res()[8][[1]],
                                    res()[9][[1]]
                                    )
                          )
      Sys.sleep(2)
      write.csv(x=res.out, file = fname, sep = ",", row.names = FALSE)
    }
  )
  session$onSessionEnded(stopApp)
}
# Shiny app finishing function -----------------------------------------
shinyApp(ui=ui, server=server)