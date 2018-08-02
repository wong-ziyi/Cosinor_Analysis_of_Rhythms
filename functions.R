popover <- function(
  title, 
  content, 
  header = NULL, 
  html = TRUE,
  class = "btn-link",                    
  placement = c('right', 'top', 'left', 'bottom'),
  trigger = c('click', 'hover', 'focus', 'manual')) {
  
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