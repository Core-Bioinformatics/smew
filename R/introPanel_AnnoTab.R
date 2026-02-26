#' Introduction Panel - Annotation Module
#' @description Simple browser for the peak/metabolite annotation table.
#' Displays full annotation details with options to download as CSV.
#'
#' @details Provides read-only access to annotation data with interactive
#' DataTable interface. Columns displayed: m_z, display_name, name, etc.
#'
#' @keywords internal
#' @name intro_anno_panel

#' @rdname intro_anno_panel
#' @export
IntroAnnoUI <- function(id, bulk.metadata, show = TRUE){
  ns <- shiny::NS(id)
  if(show){
    shiny::tabPanel(
      'Annotation table',
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

#' @rdname intro_anno_panel
#' @export
IntroAnnoServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno){
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

# QCpanelApp <- function(){
#   shinyApp(
#     ui = fluidPage(QCpanelMetabUI('qc', bulk.metadata)),
#     server = function(input, output, session){
#       QCpanelMetabServer('qc', bulk.intensity.matrix[[1]], bulk.metadata[[1]])
#     }
#   )
# }
