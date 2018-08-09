# Estimate period by iterative algorithm (modified from the fucntion, periodogram, in cosinor2 package)
periodogram_wzy<-function(data, timecol, firstsubj, lastsubj, na.action){
  periodogram<-matrix()
  #convert hour to minutes
  data[, timecol]<-data[, timecol]*60
  end<-ceiling(last(data[,timecol]))
  start<-data[, timecol][3]-data[, timecol][1]
  if (lastsubj - firstsubj == 0) {
    colnames(data)[timecol]<-"Time"
    colnames(data)[firstsubj]<-"Subjy"
    for (i in seq(from=start, to=end, by=1)) {
      cosinor<-cosinor.lm(Subjy~time(Time),data=data,na.action=na.action,period=i)
      periodogram[[i]]<-cosinor.PR(cosinor)[[2]]
    }
  }
  else {
    for (i in seq(from=start, to=end, by=1)){
      cosinor<-population.cosinor.lm(data = data, timecol = timecol, firstsubj = firstsubj, lastsubj = lastsubj, na.action = na.action, period = i)
      periodogram[[i]]<-cosinor.PR(cosinor)[[2]]
    }
  }
  periods<-as.numeric(periodogram)
  periodogram<-data.frame(periodogram)
  rows<-(nrow(periodogram))
  plot<-ggplot(periodogram,aes(x=1:rows))+
    geom_point(aes(y=periodogram)) +
    geom_line(aes(y=periodogram)) +
    labs(x = "Period (minute)", y = "Coefficient of determination")
  best<-(which(periods == max(periods,na.rm=T)))
  print(paste("The best fitting period is",best))
  return(plot)
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