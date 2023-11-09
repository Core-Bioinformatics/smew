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

#' @rdname QCpanel
#' @export
IntroAnnopanelServer <- function(id, bulk.expression.matrix, bulk.metadata, anno){
  ns <- NS(id)
  # check whether inputs (other than id) are reactive or not
  moduleServer(id, function(input, output, session){
    print(anno)
    output[['anno']] <- DT::renderDT(anno[,1:4])
  })
}

# QCpanelApp <- function(){
#   shinyApp(
#     ui = fluidPage(QCpanelMetabUI('qc', bulk.metadata)),
#     server = function(input, output, session){
#       QCpanelMetabServer('qc', bulk.expression.matrix[[1]], bulk.metadata[[1]])
#     }
#   )
# }
