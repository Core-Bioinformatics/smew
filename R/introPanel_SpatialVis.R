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
          dropMenu(
            circleButton(ns("info_peak"), icon = icon("info"),status = "success"),
            tags$div(
              tags$h3("Spatial peak visualisation"),
              tags$ul(
                tags$li("Select a peak and visualise its intensity across the samples selected at the top of this panel."),
                tags$li("The intensity can be capped at different percentiles using the sliders provided. If you select 0 and 100 as the limits then no capping is applied."),
                tags$li("Density plots showing the variation in intensity are also included to help select meaningful capping points."),
                tags$li("Log (base 2) transformations can also be applied to intensities (after capping is applied).")
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
          selectInput(ns("peakName"), "Peaks to include:", multiple = FALSE, choices = character(0)),
          checkboxInput(ns("log2_intensity"),label = 'Show log2 intensity',value = FALSE),
          sliderInput(ns("capRange"), "Cap scale on percentiles:",
                      min = 0, max = 100,
                      value = c(5,95)),
          checkboxInput(ns('splitDensity'),value = T,label = 'Split density plot by sample'),
          actionButton(ns("go_plot_peak"),'Show spatial visualisation'),
          div(style = "margin-top:10px"),
          dropMenu(
            circleButton(ns("downloadsMetadata"), icon = icon("download"),status = "success"),
            tags$div(
              tags$h3("Downloads"),
              fluidRow(
                column(5,offset=0,
                       tags$h4("Spatial distribution"),
                       textInput(ns('spatialFileName'),'File name for download', value ='spatial.png', placeholder = 'spatial.png'),
                       numericInput(ns('spatialWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('spatialHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadSpatial'), 'Download spatial figure')),
                column(5,offset=1,
                       tags$h4("Density plot"),
                       textInput(ns('densityFileName'),'File name for download', value ='density.png', placeholder = 'density.png'),
                       numericInput(ns('densityWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('densityHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadDensity'), 'Download density plot')),

                )),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
          ),
      mainPanel(
        plotOutput(ns('peakDensity')),
      fluidRow(column=10,plotOutput(ns('plotPeak'),click = ns('peak_click'))))),
      sidebarLayout(
        sidebarPanel(
          dropMenu(
            circleButton(ns("info_metadata"), icon = icon("info"),status = "success"),
            tags$div(
              tags$h3("Spatial peak visualisation"),
              tags$ul(
                tags$li("Select a metadata column and visualise it spatially across the samples selected at the top of this panel."),
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
            selectInput(ns("metadataName"), "Metadata to display:", multiple = FALSE, choices = colnames(full.metadata),selected=colnames(full.metadata)[length(colnames(full.metadata))]),
            actionButton(ns("go_plot_metadata"),'Show spatial visualisation'),
          div(style = "margin-top:10px"),
          dropMenu(
            circleButton(ns("downloads"), icon = icon("download"),status = "success"),
            tags$div(
              tags$h3("Downloads"),
              fluidRow(
                column(10,offset=0,
                       tags$h4("Metadata distribution"),
                       textInput(ns('metaFileName'),'File name for download', value ='metadata.png', placeholder = 'metadata.png'),
                       numericInput(ns('metaWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('metaHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadMeta'), 'Download metadata figure'))
              )),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          )),
        mainPanel(
          plotOutput(ns('plotMetadata'),click = ns('metadata_click'))),

      ))
  }else{
    NULL
  }
}

#' @rdname IntroSpatialVisPanel
#' @export
IntroSpatialVisPanelServer <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno){

  moduleServer(id, function(input, output, session){
    updateSelectizeInput(session, "peakName", choices = anno$display_name, server = TRUE, selected = anno$display_name[1])

    show_peak <- reactive({
      my_peak = anno[anno$display_name==input[['peakName']],]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
      current.intensity.matrix <- t(full.intensity.matrix)[,full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']]]
      caps = quantile(current.intensity.matrix[my_peak$m_z,],probs=input[['capRange']]/100)
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      current.metadata$peak = current.intensity.matrix[my_peak$m_z,]
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
        ggplot2::facet_wrap(~current.metadata$Sample, scales = 'free',nrow=max(1,floor(sqrt(length(input[['samplesToShow']])/2))))  +
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

    show_peak_zoom <- reactive({
      my_peak = anno[anno$display_name==input[['peakName']],]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] == input$peak_click$panelvar1,]
      current.intensity.matrix <- t(full.intensity.matrix)[,full.metadata[,colnames(bulk.metadata)[1]] == input$peak_click$panelvar1]
      caps = quantile(current.intensity.matrix[my_peak$m_z,],probs=input[['capRange']]/100)
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      current.metadata$peak = current.intensity.matrix[my_peak$m_z,]
      current.metadata$peak = pmin(caps[2],current.metadata$peak)
      current.metadata$peak = pmax(caps[1],current.metadata$peak)
      if (input[['log2_intensity']]){
        current.metadata$peak = log2(current.metadata$peak+1)
      } else {
        current.metadata$peak = current.metadata$peak
      }
      legend_title = ifelse(input[['log2_intensity']],paste0('log2 ',my_peak$m_z,'\n intensity'),paste0(my_peak$m_z,'\n intensity'))
      return(ggplot2::ggplot(current.metadata,ggplot2::aes(x=x,y=y,color=peak,fill=peak))+geom_tile()+
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
        ggplot2::theme(aspect.ratio = 1))
    }) %>% bindEvent(input[["go_plot_peak"]])


    show_metadata <- reactive({
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
      current.metadata$metadata = current.metadata[,input[['metadataName']]]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      return(ggplot2::ggplot(current.metadata,ggplot2::aes(x=x,y=y,color=metadata,fill=metadata))+
               ggplot2::geom_tile()+
               ggplot2::facet_wrap(~current.metadata$Sample, nrow = max(1,floor(sqrt(length(input[['samplesToShow']])/2))), scales = 'free')  +
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

    show_metadata_zoom <- reactive({
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] == input$metadata_click$panelvar1,]
      current.metadata$metadata = current.metadata[,input[['metadataName']]]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      return(ggplot2::ggplot(current.metadata,ggplot2::aes(x=x,y=y,color=metadata,fill=metadata))+
               ggplot2::geom_tile()+
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
      show_peak()$spatial})

    output[['plotPeakZoom']] <- renderPlot({
      show_peak_zoom()})

    output[['downloadSpatial']] <- downloadHandler(
      filename = function() { input[['spatialFileName']] },
      content = function(file) {
        ggsave(file, plot = show_peak()$spatial, dpi = 300,
               width=input[['spatialWidth']],height=input[['spatialHeight']])
      }
    )
    output[['peakDensity']] <- renderPlot({
      show_peak()$density
    })

    output[['downloadDensity']] <- downloadHandler(
      filename = function() { input[['densityFileName']] },
      content = function(file) {
        ggsave(file, plot = show_peak()$density, dpi = 300,
               width=input[['densityWidth']],height=input[['densityHeight']])
      }
    )
    output[['plotMetadata']] <- renderPlot({
      show_metadata()
    })

    output[['plotMetadataZoom']] <- renderPlot({
      show_metadata_zoom()
    })

    output[['downloadMeta']] <- downloadHandler(
      filename = function() { input[['metaFileName']] },
      content = function(file) {
        ggsave(file, plot = show_metadata(), dpi = 300,
               width=input[['metaWidth']],height=input[['metaHeight']])
      }
    )

    observeEvent(input$peak_click, {
      ns <- session$ns
      showModal(
        modalDialog(
          plotOutput(ns('plotPeakZoom')),
          easyClose = TRUE,
          footer = NULL
        )
      )
    })

    observeEvent(input$metadata_click, {
      ns <- session$ns
      showModal(
        modalDialog(
          plotOutput(ns('plotMetadataZoom')),
          easyClose = TRUE,
          footer = NULL
        )
      )
    })

  })
}
