# Required Packages and Functions -----------------------------------------
list.of.packages <- c("shiny", "magrittr","ggplot2", "fields", "data.table", "stringr", "ape", "DT", 
                      "shinyjs", "cosinor2", "xlsx")
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
          numericInput('Interval', 'Interval of Data', 4, min = 0),
          numericInput('XInterval', 'X-axis Display Interval', 4, min=1),
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
        #Set up x-axis label of ticks
        if(input$Interval != 0){
          TimeSeq<-seq(from=first(raw$Time), to=24, by=input$Interval) # Make general time sequence 0-24
          if(last(TimeSeq)==24){
            TimeSeq<-TimeSeq[-length(TimeSeq)]
          } # Make general time sequence 0->24
        } else {
          TimeSeq<-(0:23) # If interval un defined, then make defaut general time sequence 0-23
        }
        position<-match(raw[1, TimeTagsCol], TimeSeq) # Identify the start position of time in upload data
        ticks<-rep(TimeSeq, ceil(last(raw[, TimeTagsCol])/24)+1) # Make final sequence for labling the x-axis ticks
        #Make temporal data that contains the calculated mean and standrad deviation
        if(input$ValueCols != input$ValueCole){
          temp<-data.frame(Time=raw[, TimeTagsCol], Value=rowMeans(raw[, ValueCol]), SD=rowSds(as.matrix(raw[, ValueCol])))
        } else if(input$ValueCols == input$ValueCole){
          temp<-data.frame(Time=raw[, TimeTagsCol], Value=raw[, ValueCol], SD=0)
        }
        #Estimate the period by modified iterative function from cosinor2
        period<-periodogram_wzy(data = temp, timecol = 1, firstsubj = 2, lastsubj = 2)
        period<-period$plot_env$best # Pass the best results
        #Get best fitted cosinor model
        fit_per_est<-cosinor.lm(Value~time(Time), period = period, data = temp)
        #Calculate the estimated fitted value
        FitCurve<-data.frame(x=seq(from=first(raw[, TimeTagsCol]), to=last(raw[, TimeTagsCol]), by=0.5), y=fit_per_est$coefficients[1]+fit_per_est$coefficients[2]*cos(2*pi*seq(from=first(raw[, TimeTagsCol]), to=last(raw[, TimeTagsCol]), by=0.5)/period+pi-fit_per_est$coefficients[3]))
        #Statistically detect the rhythms
        res<-cosinor.detect(fit_per_est)
        #Make temporal for scatter plot
        ForScatter<-melt(raw[, c(TimeTagsCol, ValueCol)], id.vars=TimeTagsCol)
        ForScatter<-cbind(ForScatter, test=fit_per_est$coefficients[1]+fit_per_est$coefficients[2]*cos(2*pi*ForScatter$Time/period+pi-fit_per_est$coefficients[3]))
        #Results Part (combine all results and output)
        CurveFun<-paste0(round(fit_per_est$coefficients[1],2),
                         " + ", round(fit_per_est$coefficients[2],2),
                         "cos(2\u03C0t/", period, " + ", round(pi-fit_per_est$coefficients[3],2),
                         ")") # cosinor model
        F.statistic<-res[1] # F statistics
        rhythm.p<-res[4] # P-value for F statistics
        R2<-cor(ForScatter$value, ForScatter$test)^2 # Coefficient of Determination (GOF, goodness of fit) which be calculated from Perason's correlation coefficient
        P.value<-cor.test(ForScatter$value, ForScatter$test)$p.value # Significance for this Perason's correlation coefficient
        res<-list(raw, temp, FitCurve, ForScatter, CurveFun, F.statistic, rhythm.p, R2, P.value, position, ticks, fit_per_est, period) # Built output results
      } #UI effect: error indicator. end.
    }) #UI effect: busy indicator. end. 
    return(res)
  })
  plot.out<-reactive({
    if (is.null(res())){
      return(NULL)
    } else {
      #Pass value to local
      position<-res()[10][[1]]
      ticks<-res()[11][[1]]
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
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=input$XInterval)], labels=ticks[position:(length(temp$Time)-1+position)][seq(from=1, to=length(temp$Time), by=input$XInterval)])
      } else if (input$style == FALSE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=Value), temp)+
          geom_errorbar(aes(x=Time, ymax=Value+SD, ymin=Value-SD), temp)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=input$XInterval)], labels=ticks[position:(length(temp$Time)-1+position)][seq(from=1, to=length(temp$Time), by=input$XInterval)])
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
                          results=c(res()[12][[1]]$coefficients[1],
                                    res()[12][[1]]$coefficients[2],
                                    res()[13][[1]],
                                    res()[12][[1]]$coefficients[3],
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