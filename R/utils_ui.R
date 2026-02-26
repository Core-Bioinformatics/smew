utils::globalVariables(c(".data"))

# UI utility helpers for SMEW Shiny app

##' Update selectize inputs for peaks with sensible defaults
##'
##' Updates a Shiny selectize input for peaks, preselecting a specified number of initial choices if available.
##'
##' @param session Shiny session object.
##' @param inputId Character; input id of the selectize control.
##' @param choices Character vector of choices.
##' @param selected_n Integer; how many initial choices to preselect (default 2).
##' @keywords shiny, selectize, peaks
update_peak_selectize <- function(session, inputId, choices, selected_n = 2){
  shiny::validate(shiny::need(length(choices) > 0, paste0('No choices available for ', inputId)))
  selected_vals <- utils::head(choices, selected_n)
  shiny::updateSelectizeInput(session, inputId, choices = choices, server = TRUE, selected = selected_vals)
}

##' Convenience wrapper for a consistent info drop menu
##'
##' Creates a consistent info drop menu using shinyWidgets::dropMenu.
##'
##' @param ns Namespace function.
##' @param id Character; base id for the info button.
##' @param title Character; title string for the menu.
##' @param items Character vector of bullet items.
##' @keywords shiny, info, menu
info_dropmenu <- function(ns, id, title, items){
  shinyWidgets::dropMenu(
    shinyWidgets::circleButton(ns(id), icon = shiny::icon('info-circle'), status = 'success'),
    shiny::tags$div(
      shiny::tags$h3(title),
      shiny::tags$ul(lapply(items, function(x) shiny::tags$li(x)))
    ),
    theme = 'light-border',
    placement = 'right',
    arrow = FALSE
  )
}
