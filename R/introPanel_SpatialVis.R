#' @rdname IntroSpatialVisPanel
#' @export
IntroSpatialVisPanelUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Spatial visualisation',
      selectInput(
            inputId = ns("samplesToShow"),
            label = "Select samples to show",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[,1])[1],
            multiple = TRUE
          ),
      sidebarLayout(
        sidebarPanel(
          selectInput(ns("peakName"), "Peaks to include:", multiple = FALSE, choices = character(0)),
          checkboxInput(ns("log2_intensity"),label = 'Show log2 intensity',value = FALSE),
          sliderInput(ns("capRange"), "Cap scale on percentiles:",
                      min = 0, max = 100,
                      value = c(5,95)),
          checkboxInput(ns('splitDensity'),value = T,label = 'Split density plot by sample'),
          actionButton(ns("go_plot_peak"),'Show spatial visualisation'),
        ),
      mainPanel(
        plotOutput(ns('peakDensity')))),
      fluidRow(column=10,plotOutput(ns('plotPeak'),height = 600)),
      sidebarLayout(
        sidebarPanel(
            selectInput(ns("metadataName"), "Metadata to display:", multiple = FALSE, choices = colnames(full.metadata),selected=colnames(full.metadata)[length(colnames(full.metadata))]),
            actionButton(ns("go_plot_metadata"),'Show spatial visualisation')),
        mainPanel(
          plotOutput(ns('plotMetadata'))),

      ))
  }else{
    NULL
  }
}

#' @rdname IntroSpatialVisPanel
#' @export
IntroSpatialVisPanelServer <- function(id, bulk.metadata, full.metadata, full.expression.matrix, anno){

  moduleServer(id, function(input, output, session){
    updateSelectizeInput(session, "peakName", choices = anno$display_name, server = TRUE, selected = anno$display_name[1])

    show_peak <- reactive({
      my_peak = anno[anno$display_name==input[['peakName']],]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
      current.expression.matrix <- t(full.expression.matrix)[,full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']]]
      print(input[['capRange']])
      caps = quantile(current.expression.matrix[my_peak$m_z,],probs=input[['capRange']]/100)
      print(caps)
      print(caps[1])
      print(caps[2])
      print(max(current.expression.matrix[my_peak$m_z,]))
      print(min(current.expression.matrix[my_peak$m_z,]))
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      current.metadata$peak = current.expression.matrix[my_peak$m_z,]
      if (input[['splitDensity']]){
        density.plot = ggplot2::ggplot(current.metadata,ggplot2::aes(x=peak,color=Sample))+geom_density()+theme_classic()+geom_vline(xintercept = caps[1])+geom_vline(xintercept = caps[2])
      } else {
        density.plot = ggplot2::ggplot(current.metadata,ggplot2::aes(x=peak))+geom_density()+theme_classic()+geom_vline(xintercept = caps[1])+geom_vline(xintercept = caps[2])
      }
      current.metadata$peak = pmin(caps[2],current.metadata$peak)
      current.metadata$peak = pmax(caps[1],current.metadata$peak)
      if (input[['log2_intensity']]){
        current.metadata$peak = log2(current.metadata$peak+1)
      } else {
        current.metadata$peak = current.metadata$peak
      }
      legend_title = ifelse(input[['log2_intensity']],paste0('log2 ',my_peak$m_z,'\n intensity'),paste0(my_peak$m_z,'\n intensity'))
      spatial.plot = ggplot2::ggplot(current.metadata,ggplot2::aes(x=x,y=y,color=peak,fill=peak))+geom_tile()+
        ggplot2::facet_wrap(~current.metadata$Sample, scales = 'free',ncol=floor(2*sqrt(length(input[['samplesToShow']]))))  +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(legend_title))+
        ggplot2::scale_fill_gradient(name=legend_title,low = "lightgrey", high = "brown")+
        ggplot2::scale_color_gradient(name=legend_title,low = "lightgrey", high = "brown")+
        ggplot2::theme(aspect.ratio = 1)
      return(list('spatial'=spatial.plot,'density'=density.plot))
    }) %>% bindEvent(input[["go_plot_peak"]])

    show_metadata <- reactive({
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
      current.metadata$metadata = current.metadata[,input[['metadataName']]]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      return(ggplot2::ggplot(current.metadata,ggplot2::aes(x=x,y=y,color=metadata,fill=metadata))+
               ggplot2::geom_tile()+
               ggplot2::facet_wrap(~current.metadata$Sample, ncol = floor(2*sqrt(length(input[['samplesToShow']]))), scales = 'free')  +
               ggplot2::theme_classic() +
               ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank(),
                     axis.line.x = ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_blank())+
               ggplot2::theme(aspect.ratio = 1))
    }) %>% bindEvent(input[["go_plot_metadata"]])


    output[['plotPeak']] <- renderPlot({
      show_peak()$spatial},height=600)
    output[['peakDensity']] <- renderPlot({
      show_peak()$density
    })
    output[['plotMetadata']] <- renderPlot({
      show_metadata()
    })

    # if (input[["pickShownSamples"]]){
    #   return(input[["samplesToShow"]])
    # } else {
      return(unique(bulk.metadata[,1]))
 #   }
  })
}
