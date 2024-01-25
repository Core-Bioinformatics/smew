#' @rdname IntroAnnoPanel
#' @export
IntroAnnopanelUI <- function(id, bulk.metadata, show = TRUE){
  ns <- NS(id)
  if(show){
    tabPanel(
      'Annotation table',
      DT::dataTableOutput(ns('anno')),
      dropMenu(
        circleButton(ns("downloads"), icon = icon("download"),status = "success"),
        tags$div(
          tags$h3("Downloads"),
          fluidRow(
            column(10,offset=0,
                   tags$h4("Annotation table"),
                   textInput(ns('annoFileName'),'File name for download', value ='anno.csv', placeholder = 'anno.csv'),
                   downloadButton(ns('downloadAnno'), 'Download annotation table'))
          )),
        theme = "light-border",
        placement = "right",
        arrow = FALSE
      ))
  }else{
    NULL
  }
}

#' @rdname IntroAnnoPanel
#' @export
IntroAnnopanelServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno){
  ns <- NS(id)
  # check whether inputs (other than id) are reactive or not
  moduleServer(id, function(input, output, session){
    output[['anno']] <- DT::renderDT(anno[,1:4])

    output[['downloadAnno']] <- downloadHandler(
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
