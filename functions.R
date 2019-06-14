# Estimate period by iterative algorithm (modified from the fucntion, periodogram, in cosinor2 package)
periodogram_wzy<-function(data, timecol, firstsubj, lastsubj, iteration="NA", PFix = "NA"){
  periodogram<-c()
  periods<-c()
  #convert hour to minutes
  if(iteration!="NA"){
    end<-iteration
    start<-c(as.matrix(data[, timecol]))[4]-c(as.matrix(data[, timecol]))[1]
    if (lastsubj - firstsubj == 0) {
      colnames(data)[timecol]<-"Time"
      colnames(data)[firstsubj]<-"Subjy"
      for (i in seq(from=start, to=end, by=0.01)) {
        tryCatch({
          cosinor<-cosinor.lm(Subjy~time(Time),data=data,period=i)
          periodogram<-c(periodogram, cosinor.PR(cosinor)[[2]])
          periods<-c(periods,i)
        }, error=function(e){})
      }
    } else {
      for (i in seq(from=start, to=end, by=0.01)){
        tryCatch({
          TempData<-data[, c(timecol, firstsubj:lastsubj)]
          cosinor<-population.cosinor.lm(data = t(TempData[,c(firstsubj:lastsubj)]), time = TempData[,timecol], period = i, plot = FALSE)
          periodogram<-c(periodogram, cosinor.PR(cosinor)[[2]])
          periods<-c(periods,i)
        }, error=function(e){})
      }
    }
    df<-as.data.frame(cbind(period=periods, fit=periodogram))
    plot<-ggplot(df,aes(x=period, y=fit))+
      geom_point(aes(y=fit)) +
      geom_line(aes(y=fit)) +
      labs(x = "Period (hour)", y = "Coefficient of determination")
    best<-(periods[which(periodogram == max(periodogram,na.rm=T))])
    print(paste("The best fitting period is",best))
    return(plot)
  } else if(PFix!="NA"){
    df<-as.data.frame(cbind(period=c(1:12), fit=rep(0,12)))
    plot<-ggplot(df,aes(x=period, y=fit))+
      geom_point(aes(y=fit)) +
      geom_line(aes(y=fit)) +
      labs(x = "Period (hour)", y = "Coefficient of determination")
    best<-PFix
    print(paste("The best fitting period is",best))
    return(plot)
  }
}
rhythms_wzy<-function(raw, TimeTagsCol, Cols, Cole, ValueCol, xtitle, ytitle, design, iteration = "NA", PFix = "NA"){
  if(design==1 & Cols!=Cole){
    raw[,ValueCol]<-t(apply(raw[, ValueCol], 1, function(x)sort(x)))
  }
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
  #Make temporal data that contains the calculated mean and standrad deviation
  if(Cols != Cole){
    temp<-data.frame(Time=raw[, TimeTagsCol], Value=rowMeans(raw[, ValueCol]), SD=rowSds(as.matrix(raw[, ValueCol])))
  } else if(Cols == Cole){
    temp<-data.frame(Time=raw[, TimeTagsCol], Value=raw[, ValueCol], SD=0)
  }
  #Estimate the period by modified iterative function from cosinor2
  period0<-periodogram_wzy(data = raw, timecol = TimeTagsCol, firstsubj = Cols, lastsubj = Cole, iteration = iteration, PFix = PFix)
  period<-period0$plot_env$best # Pass the best results
  #Get best fitted cosinor model
  if(Cols==Cole){
    fit_per_est<-cosinor.lm(Value~time(Time), period = period, data = temp)
    phase<-correct.acrophase(fit_per_est)
  } else if(Cols!=Cole) {
    TempRaw<-raw[, c(TimeTagsCol, Cols:Cole)]
    fit_per_est<-population.cosinor.lm(data = t(TempRaw[,c(Cols:Cole)]), time = TempRaw[,TimeTagsCol], period = period, plot = FALSE)
    phase<-0
    for (j in 1:(Cole-Cols+1)) {
      phase<-phase+correct.acrophase(fit_per_est[["single.cos"]][[j]])
    }
    phase<-phase/(Cole-Cols+1)
  }
  #Calculate the estimated fitted value
  FitCurve<-data.frame(x=seq(from=first(raw[, TimeTagsCol]), to=last(raw[, TimeTagsCol]), by=0.01), y=as.numeric(fit_per_est$coefficients[1])+as.numeric(fit_per_est$coefficients[2])*cos(2*pi*seq(from=first(raw[, TimeTagsCol]), to=last(raw[, TimeTagsCol]), by=0.01)/period+phase))
  #Statistically detect the rhythms
  res<-cosinor.detect(fit_per_est)
  #Make temporal for scatter plot
  ForScatter<-melt(raw[, c(TimeTagsCol, ValueCol)], id.vars=TimeTagsCol)
  colnames(ForScatter)<-c("Time", "variable","value")
  ForScatter<-cbind(ForScatter, test=as.numeric(fit_per_est$coefficients[1])+as.numeric(fit_per_est$coefficients[2])*cos(2*pi*ForScatter$Time/period+phase))
  #Results Part (combine all results and output)
  CurveFun<-paste0(format(round(fit_per_est$coefficients[1],2), nsmall=2),
                   " + ", format(round(fit_per_est$coefficients[2],2), nsmall=2),
                   "cos(2\u03C0t/", format(round(period, 2), nsmall=2), with_plus(format(round(phase,2)), nsmall=2),
                   ")") # cosinor model
  F.statistic<-res[1] # F statistics
  rhythm.p<-res[4] # P-value for F statistics
  R2<-cor(ForScatter$value, ForScatter$test)^2 # Coefficient of Determination (GOF, goodness of fit) which be calculated from Perason's correlation coefficient
  P.value<-cor.test(ForScatter$value, ForScatter$test)$p.value # Significance for this Perason's correlation coefficient
  res.out<-list(raw,           #1
                temp,          #2
                FitCurve,      #3
                ForScatter,    #4
                CurveFun,      #5
                F.statistic,   #6
                rhythm.p,      #7
                R2,            #8
                P.value,       #9
                ticks,         #10
                fit_per_est,   #11
                period,        #12
                period0,       #13
                phase          #14
                ) # Built output results
  return(res.out)
}
# Some UI functions
popover <- function(title, content, header = NULL, html = TRUE, class = "btn-link", placement = c('right', 'top', 'left', 'bottom'), trigger = c('click', 'hover', 'focus', 'manual')) {
  
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
      )
    ),
    tags$a(
      tabindex = "0", href = NULL, role = "button", class = class, `data-toggle` = "popover",
      title = header, `data-content` = content, `data-animation` = TRUE, 'data-html' = html,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      title
    )
  )
}
# Need to use with the corresponding `withBusyIndicator` server function
withBusyIndicatorUI <- function(button) {
  id <- button[['attribs']][['id']]
  div(
    `data-for-btn` = id,
    button,
    span(
      class = "btn-loading-container",
      hidden(
        img(src = "ajax-loader-bar.gif", class = "btn-loading-indicator"),
        icon("check", class = "btn-done-indicator")
      )
    ),
    hidden(
      div(class = "btn-err",
          div(icon("exclamation-circle"),
              tags$b("Error: "),
              span(class = "btn-err-msg")
          )
      )
    )
  )
}
# Call this function from the server with the button id that is clicked and the
# expression to run when the button is clicked
withBusyIndicatorServer <- function(buttonId, expr) {
  # UX stuff: show the "busy" message, hide the other messages, disable the button
  loadingEl <- sprintf("[data-for-btn=%s] .btn-loading-indicator", buttonId)
  doneEl <- sprintf("[data-for-btn=%s] .btn-done-indicator", buttonId)
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  shinyjs::disable(buttonId)
  shinyjs::show(selector = loadingEl)
  shinyjs::hide(selector = doneEl)
  shinyjs::hide(selector = errEl)
  on.exit({
    shinyjs::enable(buttonId)
    shinyjs::hide(selector = loadingEl)
  })
  
  # Try to run the code when the button is clicked and show an error message if
  # an error occurs or a success message if it completes
  tryCatch({
    value <- expr
    shinyjs::show(selector = doneEl)
    shinyjs::delay(2000, shinyjs::hide(selector = doneEl, anim = TRUE, animType = "fade",
                                       time = 0.5))
    value
  }, error = function(err) { errorFunc(err, buttonId) })
}
# When an error happens after a button click, show the error
errorFunc <- function(err, buttonId) {
  errEl <- sprintf("[data-for-btn=%s] .btn-err", buttonId)
  errElMsg <- sprintf("[data-for-btn=%s] .btn-err-msg", buttonId)
  errMessage <- gsub("^ddpcr: (.*)", "\\1", err$message)
  shinyjs::html(html = errMessage, selector = errElMsg)
  shinyjs::show(selector = errEl, anim = TRUE, animType = "fade")
}
with_plus <- function(x, ...)
{
  if (x > 0)
  {
    sprintf(
      fmt = "+ %s", 
      format(x, ...)
    )
  }
  else
  {
    x
  }
}