IntroSpatialVisPanelUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Spatial visualisation',
      selectInput(
            inputId = ns("samplesToShow"),
            label = "Select samples for clustering",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[,1])[1],
            multiple = TRUE
          ),
        selectInput(ns("peakName"), "Peaks to include:", multiple = FALSE, choices = character(0)),
        checkboxInput(ns("log2_intensity"),label = 'Show log2 intensity',value = FALSE),
        plotOutput(ns('plotPeak')),
    
        selectInput(ns("metadataName"), "Metadata to display:", multiple = FALSE, choices = colnames(full.metadata),selected=colnames(full.metadata)[length(colnames(full.metadata))]),
      
        plotOutput(ns('plotMetadata')),
      )
  }else{
    NULL
  }
}

#' @rdname DEsummaryPanel
#' @export
IntroSpatialVisPanelServer <- function(id, bulk.metadata, full.metadata, full.expression.matrix, anno){

  moduleServer(id, function(input, output, session){
    updateSelectizeInput(session, "peakName", choices = anno$display_name, server = TRUE, selected = anno$display_name[1])

    show_peak <- reactive({
      print(input[['log2_intensity']])
      my_peak = anno[anno$display_name==input[['peakName']],]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
      current.expression.matrix <- t(full.expression.matrix)[,full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']]]
      print(str(current.metadata))
#      current.metadata$peak = log2(current.expression.matrix[my_peak$m_z,]+1)
      if (input[['log2_intensity']]){
        current.metadata$peak = log2(current.expression.matrix[my_peak$m_z,]+1)
      } else {
        current.metadata$peak = current.expression.matrix[my_peak$m_z,]
      }
      print(str(current.metadata))
      print(summary(current.metadata))
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      legend_title = ifelse(input[['log2_intensity']],paste0('log2 ',my_peak$m_z,'\n intensity'),paste0(my_peak$m_z,'\n intensity'))
      print(legend_title)
      return(ggplot(current.metadata,aes(x=x,y=y,color=peak,fill=peak))+geom_tile()+
               facet_wrap(~current.metadata$Sample, nrow = floor(sqrt(length(input[['samplesToShow']]))), scales = 'free')  +                  
               theme_classic() +
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.line.y = element_blank(),
                     legend.title = element_text(legend_title))+ 
        scale_fill_gradient(name=legend_title,low = "lightgrey", high = "brown")+
        scale_color_gradient(name=legend_title,low = "lightgrey", high = "brown")+ 
        theme(aspect.ratio = 1))
    })
    
    show_metadata <- reactive({
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
      current.metadata$metadata = current.metadata[,input[['metadataName']]]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      return(ggplot(current.metadata,aes(x=x,y=y,color=metadata,fill=metadata))+geom_tile()+
               facet_wrap(~current.metadata$Sample, nrow = floor(sqrt(length(input[['samplesToShow']]))), scales = 'free')  +                  
               theme_classic() +
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank(),
                     axis.line.x = element_blank(),
                     axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank(),
                     axis.line.y = element_blank())+ 
               theme(aspect.ratio = 1))
    })
    

    output[['plotPeak']] <- renderPlot({
      show_peak()
    })
    output[['plotMetadata']] <- renderPlot({
      show_metadata()
    })
  })
}