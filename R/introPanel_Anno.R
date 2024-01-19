#' @rdname IntroAnnoPanel
#' @export
IntroAnnopanelUI <- function(id, bulk.metadata, show = TRUE){
  ns <- NS(id)
  if(show){
    tabPanel(
      'Annotation table',
      DT::dataTableOutput(ns('anno'))
    )
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
