
#' Presents metabolite annotations in an interactive table
#'
#' @description UI and server logic for browsing the peak/metabolite annotation table. Displays full annotation details with options to download as CSV.
#'
#' @details
#' \itemize{
#'  \item{Provides read-only access to annotation data with an interactive DataTable interface.}
#'  \item{Columns displayed: m_z, display_name, name, etc.}
#'  \item{Download handler for exporting the annotation table as CSV.}
#'}
#' @name IntroPanel_AnnoTab
#' @rdname IntroPanel_AnnoTab
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame with bulk sample metadata
#' @param show Logical; whether to render the panel (default: TRUE)
#' @return A shiny::tabPanel containing the UI elements for the annotation table panel
#' @export
IntroPanel_AnnoTabUI <- function(id, bulk.metadata, show = TRUE){
  ns <- shiny::NS(id)
  if(show){
    shiny::tabPanel(
      'Annotation Table',
      shinyWidgets::dropMenu(
        shinyWidgets::circleButton(ns("downloads"), icon = shiny::icon("download"),status = "success"),
        shiny::tags$div(
          shiny::tags$h3("Downloads"),
          shiny::fluidRow(
            shiny::column(10,offset=0,
                   shiny::tags$h4("Annotation table"),
                   shiny::textInput(ns('annoFileName'),'File name for download', value ='anno.csv', placeholder = 'anno.csv'),
                   shiny::downloadButton(ns('downloadAnno'), 'Download annotation table'))
          )),
        theme = "light-border",
        placement = "right",
        arrow = FALSE
      ),
      DT::dataTableOutput(ns('anno'))
)
  }else{
    NULL
  }
}

# ============ SERVER FUNCTIONS ============

#' @rdname IntroPanel_AnnoTab
#' @param bulk.intensity.matrix Numeric matrix of bulk sample intensities (features × samples)
#' @param anno Data frame with annotation information (must include columns: m_z, display_name, name, etc.)
#' @export
IntroPanel_AnnoTabServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno){
  ns <- shiny::NS(id)
  # check whether inputs (other than id) are reactive or not
  shiny::moduleServer(id, function(input, output, session){
    output[['anno']] <- DT::renderDT(anno[,1:4])

    output[['downloadAnno']] <- shiny::downloadHandler(
      filename = function() {
        paste(input[['annoFileName']])
      },
      content = function(file) {
        utils::write.csv(x = anno[,1:4], file = file, row.names = FALSE)
      }
    )
  })
}
