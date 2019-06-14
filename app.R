# Required Packages and Functions -----------------------------------------
list.of.packages <- c("shiny", "BH","magrittr","ggplot2", "fields", "data.table", "stringr", "ape", "htmlwidgets", "DT", "Hmisc",
                      "shinyjs", "cosinor2", "spam", "matrixStats","gtable", "munsell")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages)
}
library(shiny)
library(BH)
library(magrittr)
library(matrixStats)
library(htmlwidgets)
library(fields)
library(data.table)
library(stringr)
library(ape)
library(DT)
library(shinyjs)
library(Hmisc)
library(gtable)
library(munsell)
library(ggplot2)
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
          checkboxInput('compare', "Test two groups"),
          radioButtons("PeriodSet", "Period setting:",
                       c("Set the maximal interation" = "AutoPer",
                         "A fixed value" = "PFix")),
          uiOutput("compare.in"),
          checkboxInput('style', 'Scatter plot'),
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
          uiOutput("compare.out")
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
  output$compare.in<-renderUI({
    if(input$compare==FALSE){
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
        radioButtons("sample", "Experiment design:",
                     c("Repeat measurement" = "Rep",
                       "Independent/single measurement" = "Ind")),
        if(input$PeriodSet=="AutoPer"){
          numericInput('interation', 'Maximal interation for period', 48, min = 3, step = 1)
        }else if(input$PeriodSet=="PFix"){
          numericInput('PeriodFix', 'Fixed Period', 24, min = 3, step = 1)
        },
        numericInput('TimeTagsCol', 'Cloumn number of time', 1),
        numericInput('ValueCols', 'Start cloumn number of sample', 2, min = 1, step = 1),
        numericInput('ValueCole', 'End Cloumn number of sample', 4, min = 1, step = 1),
        numericInput('XInterval', 'X-axis Display Interval', 1, min=1),
        textInput('xtitle', 'Label of x-axis', "Time(h)"),
        textInput('ytitle', 'Label of y-axis', "Expression Level")
      )
    } else if(input$compare==TRUE){
      splitLayout(
        verticalLayout(
          fileInput(
            'raw1',
            label = "Upload your Group1 in .csv file:",
            accept = c('csv', 'comma-separated-values', '.csv')
          ),
          radioButtons("sample1", "Experiment design:",
                       c("Repeat measurement" = "Rep",
                         "Independent/single measurement" = "Ind")),
          if(input$PeriodSet=="AutoPer"){
            numericInput('interation1', 'Maximal interation for period', 48, min = 3, step = 1)
          }else if(input$PeriodSet=="PFix"){
            numericInput('PeriodFix1', 'Fixed Period', 24, min = 3, step = 1)
          },
          numericInput('TimeTagsCol1', 'Cloumn number of time', 1),
          numericInput('ValueCols1', 'Start cloumn number of sample', 2, min = 1, step = 1),
          numericInput('ValueCole1', 'End Cloumn number of sample', 4, min = 1, step = 1),
          numericInput('XInterval1', 'X-axis Display Interval', 1, min=1),
          textInput('xtitle1', 'Label of x-axis', "Time(h)"),
          textInput('ytitle1', 'Label of y-axis', "Expression Level")
        ),
        verticalLayout(
          fileInput(
            'raw2',
            label = "Upload your Group2 in .csv file:",
            accept = c('csv', 'comma-separated-values', '.csv')
          ),
          radioButtons("sample2", "Experiment design:",
                       c("Repeat measurement" = "Rep",
                         "Independent/single measurement" = "Ind")),
          if(input$PeriodSet=="AutoPer"){
            numericInput('interation2', 'Maximal interation for period', 48, min = 3, step = 1)
          }else if(input$PeriodSet=="PFix"){
            numericInput('PeriodFix2', 'Fixed Period', 24, min = 3, step = 1)
          },
          numericInput('TimeTagsCol2', 'Cloumn number of time', 1),
          numericInput('ValueCols2', 'Start cloumn number of sample', 2, min = 1, step = 1),
          numericInput('ValueCole2', 'End Cloumn number of sample', 4, min = 1, step = 1),
          numericInput('XInterval2', 'X-axis Display Interval', 1, min=1),
          textInput('xtitle2', 'Label of x-axis', "Time(h)"),
          textInput('ytitle2', 'Label of y-axis', "Expression Level")
        )
      )
    }
  })
  output$compare.out<-renderUI({
    if(input$compare==FALSE){
      verticalLayout(
        plotOutput("plot.out"),
        textOutput("F.statistic"),
        textOutput("rhythm.p"),
        htmlOutput("CurveFun"),
        textOutput("R2"),
        textOutput("P.value"),
        plotOutput("plot.period")
      )
    }else if(input$compare==TRUE){
      verticalLayout(
        tags$h4("F-test (Bingham et al. 1982)"),
        dataTableOutput("table"),
        tags$h4("Student t-test (independently calculate parameters for each subject)"),
        dataTableOutput("table.t"),
        tags$hr(),
        splitLayout(
          verticalLayout(
            tags$h4("Group 1"),
            plotOutput("plot.out1"),
            textOutput("F.statistic"),
            textOutput("rhythm.p"),
            htmlOutput("CurveFun"),
            textOutput("R2"),
            textOutput("P.value"),
            dataTableOutput("table.G1"),
            dataTableOutput("table.CI.G1"),
            plotOutput("plot.period")
          ),
          verticalLayout(
            tags$h4("Group 2"),
            plotOutput("plot.out2"),
            textOutput("F.statistic2"),
            textOutput("rhythm.p2"),
            htmlOutput("CurveFun2"),
            textOutput("R22"),
            textOutput("P.value2"),
            dataTableOutput("table.G2"),
            dataTableOutput("table.CI.G2"),
            plotOutput("plot.period2")
          )
        )
      )
    }
  })
  res<-eventReactive(input$OK, {
    withBusyIndicatorServer("OK", {
      if (is.null(input$raw$datapath) & (is.null(input$raw1$datapath)|is.null(input$raw2$datapath))) {
        stop("Please submit your data first")
        return(NULL)
      } else {
        if(input$compare==FALSE){
          raw<-read.csv(input$raw$datapath, header = TRUE, sep = ",")
          #Parameters setting (Passed from UI input)
          TimeTagsCol<-input$TimeTagsCol #Column number of the time series
          Cols<-input$ValueCols
          Cole<-input$ValueCole
          ValueCol<-c(Cols:Cole) #Column number of the results
          xtitle<-input$xtitle #Title of x-axis
          ytitle<-input$ytitle #Title of y-axis
          design<-switch(input$sample,
                         Rep=1,
                         Ind=2)
          if(input$PeriodSet=="PFix"){
            res.out<-rhythms_wzy(raw, TimeTagsCol, Cols, Cole, ValueCol, xtitle, ytitle, design, iteration = "NA", PFix = input$PeriodFix)
          } else if(input$PeriodSet=="AutoPer"){
            res.out<-rhythms_wzy(raw, TimeTagsCol, Cols, Cole, ValueCol, xtitle, ytitle, design, iteration = input$interation, PFix = "NA")
          }
        } else if(input$compare==TRUE){
          ###Group1
          raw1<-read.csv(input$raw1$datapath, header = TRUE, sep = ",")
          #Parameters setting (Passed from UI input)
          TimeTagsCol<-input$TimeTagsCol1 #Column number of the time series
          Cols<-input$ValueCols1
          Cole<-input$ValueCole1
          ValueCol<-c(Cols:Cole) #Column number of the results
          xtitle<-input$xtitle1 #Title of x-axis
          ytitle<-input$ytitle1 #Title of y-axis
          design<-switch(input$sample1,
                         Rep=1,
                         Ind=2)
          if(input$PeriodSet=="PFix"){
            res1<-rhythms_wzy(raw1, TimeTagsCol, Cols, Cole, ValueCol, xtitle, ytitle, design, iteration = "NA", PFix = input$PeriodFix1)
          } else if(input$PeriodSet=="AutoPer"){
            res1<-rhythms_wzy(raw1, TimeTagsCol, Cols, Cole, ValueCol, xtitle, ytitle, design, iteration = input$interation1, PFix = "NA")
          }
          if(Cols!=Cole){
            ME1<-c()
            Am1<-c()
            Pe1<-c()
            Ph1<-c()
            for (i in Cole:Cols) {
              if(input$PeriodSet=="AutoPer"){
                temp<-rhythms_wzy(raw1, TimeTagsCol, i, i, i, xtitle, ytitle, design, iteration = input$interation1)
              } else if(input$PeriodSet=="PFix"){
                temp<-rhythms_wzy(raw1, TimeTagsCol, i, i, i, xtitle, ytitle, design, iteration = input$PeriodFix1)
              }
              ME1<-c(ME1, temp[[11]][["coefficients"]][1])
              Am1<-c(Am1, temp[[11]][["coefficients"]][2])
              Pe1<-c(Pe1, temp[[12]][1])
              Ph1<-c(Ph1, temp[[14]][1])
            }
            G1<-data.frame(
              "MESOR"=ME1,
              "Amplitude"=Am1,
              "Period"=Pe1,
              "Phase"=Ph1
            )
          }
          ###Group2
          raw2<-read.csv(input$raw2$datapath, header = TRUE, sep = ",")
          #Parameters setting (Passed from UI input)
          TimeTagsCol<-input$TimeTagsCol2 #Column number of the time series
          Cols<-input$ValueCols2
          Cole<-input$ValueCole2
          ValueCol<-c(Cols:Cole) #Column number of the results
          xtitle<-input$xtitle2 #Title of x-axis
          ytitle<-input$ytitle2 #Title of y-axis
          design<-switch(input$sample2,
                         Rep=1,
                         Ind=2)
          if(input$PeriodSet=="PFix"){
            res2<-rhythms_wzy(raw2, TimeTagsCol, Cols, Cole, ValueCol, xtitle, ytitle, design, iteration = "NA", PFix = input$PeriodFix2)
          } else if(input$PeriodSet=="AutoPer"){
            res2<-rhythms_wzy(raw2, TimeTagsCol, Cols, Cole, ValueCol, xtitle, ytitle, design, iteration = input$interation2, PFix = "NA")
          }
          if(Cols!=Cole){
            ME2<-c()
            Am2<-c()
            Pe2<-c()
            Ph2<-c()
            for (i in Cole:Cols) {
              if(input$PeriodSet=="AutoPer"){
                temp<-rhythms_wzy(raw2, TimeTagsCol, i, i, i, xtitle, ytitle, design, iteration = input$interation2)
              } else if(input$PeriodSet=="PFix"){
                temp<-rhythms_wzy(raw2, TimeTagsCol, i, i, i, xtitle, ytitle, design, iteration = input$PeriodFix2)
              }
              ME2<-c(ME2, temp[[11]][["coefficients"]][1])
              Am2<-c(Am2, temp[[11]][["coefficients"]][2])
              Pe2<-c(Pe2, temp[[12]][1])
              Ph2<-c(Ph2, temp[[14]][1])
            }
            G2<-data.frame(
              "MESOR"=ME2,
              "Amplitude"=Am2,
              "Period"=Pe2,
              "Phase"=Ph2
            )
            ME<-t.test(ME1, ME2)
            Am<-t.test(Am1, Am2)
            if(sum(Pe1-Pe2)!=0){
              Pe<-t.test(Pe1, Pe2)
            } else {
              Pe<-list()
              Pe$statistic<-0
              Pe$p.value<-0
            }
            Ph<-t.test(Ph1, Ph2)
            Student.t<-rbind("MESOR"=c(ME$statistic, ME$p.value),
                             "Amplitude"=c(Am$statistic, Am$p.value),
                             "Period"=c(Pe$statistic, Pe$p.value),
                             "Phase"=c(Ph$statistic, Ph$p.value))
            colnames(Student.t)<-c("t", "p-value")
            Student.t<-as.data.frame(Student.t)
          } else {
            Student.t<-rbind("MESOR"=c("null", "null"),
                             "Amplitude"=c("null", "null"),
                             "Period"=c("null", "null"),
                             "Phase"=c("null", "null"))
            colnames(Student.t)<-c("t", "p-value")
            Student.t<-as.data.frame(Student.t)
          }
          res0<-cosinor.poptests(res1[11][[1]], res2[11][[1]])
          res0[3,5]<-res1[[14]][1]
          res0[3,6]<-res2[[14]][1]
          rownames(res0)[3]<-"Phase"
          res.out<-list(res1,
                        res2,
                        as.data.frame(res0), 
                        G1, 
                        G2, 
                        Student.t)
          res.out[[1]][[11]]$conf.ints<-rbind(
            res.out[[1]][[11]]$conf.ints,
            variance=c(
              abs(res.out[[1]][[11]]$conf.ints[1,1]-res.out[[1]][[11]]$conf.ints[2,1])/2,
              abs(res.out[[1]][[11]]$conf.ints[1,2]-res.out[[1]][[11]]$conf.ints[2,2])/2,
              abs(res.out[[1]][[11]]$conf.ints[1,3]-res.out[[1]][[11]]$conf.ints[2,3])/2
            )
          )
          res.out[[2]][[11]]$conf.ints<-rbind(
            res.out[[2]][[11]]$conf.ints,
            variance=c(
              abs(res.out[[2]][[11]]$conf.ints[1,1]-res.out[[2]][[11]]$conf.ints[2,1])/2,
              abs(res.out[[2]][[11]]$conf.ints[1,2]-res.out[[2]][[11]]$conf.ints[2,2])/2,
              abs(res.out[[2]][[11]]$conf.ints[1,3]-res.out[[2]][[11]]$conf.ints[2,3])/2
            )
          )
        }
      } #UI effect: error indicator. end.
    }) #UI effect: busy indicator. end. 
    return(res.out)
  })
  plot.out<-reactive({
    if (is.null(res())){
      return(NULL)
    } else {
      #Pass value to local
      if(input$compare==FALSE){
        ticks<-res()[10][[1]]
        ForScatter<-res()[4][[1]]
        FitCurve<-res()[3][[1]]
        temp<-res()[2][[1]]
        xtitle<-input$xtitle
        ytitle<-input$ytitle
        intervalx<-input$XInterval
      }else if(input$compre==TRUE){
        ticks<-res()[1][[1]][10][[1]]
        ForScatter<-res()[1][[1]][4][[1]]
        FitCurve<-res()[1][[1]][3][[1]]
        temp<-res()[1][[1]][2][[1]]
        xtitle<-input$xtitle
        ytitle<-input$ytitle
        intervalx<-input$XInterval
      }
      Sys.sleep(1)
      Sys.sleep(1)
      #Scatter or Mean SD?
      if (input$style == TRUE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=value), size=2, shape=1, fill="white", color="black", ForScatter)+
          geom_line(aes(x=x, y=y, colour="red"), size=1, FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=intervalx)], labels=ticks[seq(from=1, to=length(temp$Time), by=intervalx)])
      } else if (input$style == FALSE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=Value), shape=1, fill="white", color="black", size=1.5, temp)+
          geom_errorbar(aes(x=Time, ymax=Value+SD, ymin=Value-SD), width=0.5, size=0.5, position = position_dodge(0.3), temp)+
          geom_line(aes(x=x, y=y, colour="red"), size=1, FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=intervalx)], labels=ticks[seq(from=1, to=length(temp$Time), by=intervalx)])
      }
      return(gp)
    }
  })
  plot.out1<-reactive({
    if (is.null(res()) & input$compare==FALSE){
      return(NULL)
    } else {
      #Pass value to local
      ticks<-res()[1][[1]][10][[1]]
      ForScatter<-res()[1][[1]][4][[1]]
      FitCurve<-res()[1][[1]][3][[1]]
      temp<-res()[1][[1]][2][[1]]
      xtitle<-input$xtitle
      ytitle<-input$ytitle
      intervalx<-input$XInterval
      #Scatter or Mean SD?
      if (input$style == TRUE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=value), ForScatter)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=intervalx)], labels=ticks[seq(from=1, to=length(temp$Time), by=intervalx)])
      } else if (input$style == FALSE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=Value), temp)+
          geom_errorbar(aes(x=Time, ymax=Value+SD, ymin=Value-SD), temp)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=intervalx)], labels=ticks[seq(from=1, to=length(temp$Time), by=intervalx)])
      }
      return(gp)
    }
  })
  plot.out2<-reactive({
    if (is.null(res()) & input$compare==FALSE){
      return(NULL)
    } else {
      #Pass value to local
      ticks<-res()[2][[1]][10][[1]]
      ForScatter<-res()[2][[1]][4][[1]]
      FitCurve<-res()[2][[1]][3][[1]]
      temp<-res()[2][[1]][2][[1]]
      xtitle<-input$xtitle2
      ytitle<-input$ytitle2
      intervalx<-input$XInterval2
      #Scatter or Mean SD?
      if (input$style == TRUE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=value), ForScatter)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=intervalx)], labels=ticks[seq(from=1, to=length(temp$Time), by=intervalx)])
      } else if (input$style == FALSE){
        gp <- ggplot()+
          geom_point(aes(x=Time, y=Value), temp)+
          geom_errorbar(aes(x=Time, ymax=Value+SD, ymin=Value-SD), temp)+
          geom_line(aes(x=x, y=y, colour="red"), FitCurve)+
          theme_classic()+
          theme(legend.position = "none", text = element_text(size=26))+
          labs(x=xtitle, y=ytitle)+
          expand_limits(y = 0)+
          scale_x_continuous(breaks = c(temp$Time)[seq(from=1, to=length(temp$Time), by=intervalx)], labels=ticks[seq(from=1, to=length(temp$Time), by=intervalx)])
      }
      return(gp)
    }
  })
  output$plot.out<-renderPlot(plot.out())
  output$plot.out1<-renderPlot(plot.out1())
  output$plot.out2<-renderPlot(plot.out2())
  ###Group1 output
  output$plot.period<-renderPlot({
    if(input$compare==FALSE){
      print(res()[13][[1]])
    }else if(input$compare==TRUE){
      print(res()[1][[1]][[13]])
    }
  })
  output$F.statistic<-renderText({
    if(input$compare==FALSE){
      paste0("F statistics: ", res()[6])
    }else if(input$compare==TRUE){
      paste0("F statistics: ", res()[1][[1]][[6]])
    }
  })
  output$rhythm.p<-renderText({
    if(input$compare==FALSE){
      paste0("P-value for rhythms: ", res()[7])
    }else if(input$compare==TRUE){
      paste0("P-value for rhythms: ", res()[1][[1]][[7]])
    }
  })
  output$CurveFun<-renderUI({
    if(input$compare==FALSE){
      HTML(paste0("Function of fitted curve: ", res()[5]))
    }else if(input$compare==TRUE){
      HTML(paste0("Function of fitted curve: ", res()[1][[1]][[5]]))
    }
  })
  output$R2<-renderText({
    if(input$compare==FALSE){
      paste0("R squared (Goodness of fit): ", res()[8])
    }else if(input$compare==TRUE){
      paste0("R squared (Goodness of fit): ", res()[1][[1]][[8]])
    }
  })
  output$P.value<-renderText({
    if(input$compare==FALSE){
      paste0("P-value of fit: ", res()[9])
    }else if(input$compare==TRUE){
      paste0("P-value of fit: ", res()[1][[1]][[9]])
    }
  })
  output$table.G1<-renderDataTable(res()[[4]])
  output$table.CI.G1<-renderDataTable(as.data.frame(res()[[1]][[11]]$conf.ints))
  ###Group2 output
  output$plot.period2<-renderPlot({
    if(input$compare==FALSE){
      print(res()[13][[1]])
    }else if(input$compare==TRUE){
      print(res()[2][[1]][[13]])
    }
  })
  output$F.statistic2<-renderText({
    if(input$compare==FALSE){
      paste0("F statistics: ", res()[6])
    }else if(input$compare==TRUE){
      paste0("F statistics: ", res()[2][[1]][[6]])
    }
  })
  output$rhythm.p2<-renderText({
    if(input$compare==FALSE){
      paste0("P-value for rhythms: ", res()[7])
    }else if(input$compare==TRUE){
      paste0("P-value for rhythms: ", res()[2][[1]][[7]])
    }
  })
  output$CurveFun2<-renderUI({
    if(input$compare==FALSE){
      HTML(paste0("Function of fitted curve: ", res()[5]))
    }else if(input$compare==TRUE){
      HTML(paste0("Function of fitted curve: ", res()[2][[1]][[5]]))
    }
  })
  output$R22<-renderText({
    if(input$compare==FALSE){
      paste0("R squared (Goodness of fit): ", res()[8])
    }else if(input$compare==TRUE){
      paste0("R squared (Goodness of fit): ", res()[2][[1]][[8]])
    }
  })
  output$P.value2<-renderText({
    if(input$compare==FALSE){
      paste0("P-value of fit: ", res()[9])
    }else if(input$compare==TRUE){
      paste0("P-value of fit: ", res()[2][[1]][[9]])
    }
  })
  output$table.CI.G2<-renderDataTable(as.data.frame(res()[[2]][[11]]$conf.ints))
  output$table<-renderDataTable(as.data.frame(res()[[3]]))
  output$table.G2<-renderDataTable(res()[[5]])
  output$table.t<-renderDataTable(res()[[6]])
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
                                    res()[14][[1]],
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