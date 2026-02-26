##' Spatial Visualisation Panel UI
##'
##' Provides the UI for spatial visualisation of peaks and metadata in the SMEW app.
##'
##' @param id Shiny module id
##' @param bulk.metadata Data frame of bulk sample metadata
##' @param full.metadata Data frame of full sample metadata
##' @param show Logical; whether to show the panel (default TRUE)
##' @return A shiny tabPanel object for spatial visualisation
##' @export
IntroSpatialVisUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- shiny::NS(id)

  if(show){
    shiny::tabPanel(
      'Spatial visualisation',

      shiny::selectInput(
            inputId = ns("samplesToShow"),
            label = "Select samples to show",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[,1])[1],
            multiple = TRUE
          ),
      shiny::h2('Visualise individual peaks'),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shinyWidgets::dropMenu(
            shinyWidgets::circleButton(ns("info_peak"), icon = shiny::icon("info"),status = "success"),
            shiny::tags$div(
              shiny::tags$h3("Spatial peak visualisation"),
              shiny::tags$ul(
                shiny::tags$li("Select a peak and visualise its intensity across the samples selected at the top of this panel."),
                shiny::tags$li("The intensity can be capped at different percentiles using the sliders provided. If you select 0 and 100 as the limits then no capping is applied."),
                shiny::tags$li("Density plots showing the variation in intensity are also included to help select meaningful capping points."),
                shiny::tags$li("Log (base 2) transformations can also be applied to intensities (after capping is applied) to facilitate more meaningful visualisation.")
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
          shiny::selectInput(ns("peakName"), "Peaks to include:", multiple = FALSE, choices = character(0)),
          shiny::checkboxInput(ns("log2_intensity"),label = 'Show log2 intensity',value = FALSE),
          shiny::sliderInput(ns("capRange"), "Cap scale on percentiles:",
                      min = 0, max = 100,
                      value = c(5,95)),
          shiny::checkboxInput(ns('splitDensity'),value = T,label = 'Split density plot by sample'),
          shiny::actionButton(ns("go_plot_peak"),'Show spatial visualisation'),
          shiny::div(style = "margin-top:10px"),
          shinyWidgets::dropMenu(
            shinyWidgets::circleButton(ns("downloadsMetadata"), icon = shiny::icon("download"),status = "success"),
            shiny::tags$div(
              shiny::tags$h3("Downloads"),
              shiny::fluidRow(
                shiny::column(5,offset=0,
                       shiny::tags$h4("Spatial distribution"),
                       shiny::textInput(ns('spatialFileName'),'File name for download', value ='spatial.png', placeholder = 'spatial.png'),
                       shiny::numericInput(ns('spatialWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       shiny::numericInput(ns('spatialHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       shiny::downloadButton(ns('downloadSpatial'), 'Download spatial figure')),
                shiny::column(5,offset=1,
                       shiny::tags$h4("Density plot"),
                       shiny::textInput(ns('densityFileName'),'File name for download', value ='density.png', placeholder = 'density.png'),
                       shiny::numericInput(ns('densityWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       shiny::numericInput(ns('densityHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       shiny::downloadButton(ns('downloadDensity'), 'Download density plot')),
                )),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
          ),
      shiny::mainPanel(
        plotly::plotlyOutput(ns('peakDensity')),
        shiny::fluidRow(shiny::column(10,shiny::plotOutput(ns('plotPeak'),click = ns('peak_click'),height='600px')))
        )),
      shiny::h2('Visualise multiple peaks as colour channels'),
       shiny::sidebarLayout(
         shiny::sidebarPanel(
           shinyWidgets::dropMenu(
             shinyWidgets::circleButton(ns("info_multiple_peak"), icon = shiny::icon("info"),status = "info"),
             shiny::tags$div(
               shiny::tags$h3("Multiple spatial peak visualisation"),
               shiny::tags$ul(
                 shiny::tags$li("Select 3 peaks and visualise their intensity as RGB channels across the samples selected at the top of this panel."),
                 shiny::tags$li("Log (base 2) transformations can also be applied to intensities (after capping is applied).")
               )
             ),
             theme = "light-border",
             placement = "right",
             arrow = FALSE
           ),
           shiny::selectInput(ns("peakName1"), "Peak 1 (Red):", multiple = FALSE, choices = character(0)),
           shiny::selectInput(ns("peakName2"), "Peak 2 (Green):", multiple = FALSE, choices = character(0)),
           shiny::selectInput(ns("peakName3"), "Peak 3 (Blue):", multiple = FALSE, choices = character(0)),
           shiny::checkboxInput(ns("log2_intensity_multiple"),label = 'Show log2 intensity',value = FALSE),
           shiny::sliderInput(ns("capRange_multiple"), "Cap scale on percentiles:",
                       min = 0, max = 100,
                       value = c(5,95)),
           shiny::actionButton(ns("go_plot_multiple_peaks"),'Show spatial visualisation'),
           shiny::div(style = "margin-top:10px"),
           shinyWidgets::dropMenu(
             shinyWidgets::circleButton(ns("downloads_multiple"), icon = shiny::icon("download"),status = "success"),
             shiny::tags$div(
               shiny::tags$h3("Downloads"),
               shiny::tags$h4("Multiple peak distribution"),
               shiny::textInput(ns('multiPeakFileName'),'File name for download', value ='multiPeak.png', placeholder = 'multiPeak.png'),
               shiny::numericInput(ns('multiPeakWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
               shiny::numericInput(ns('multiPeakHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
               shiny::downloadButton(ns('downloadMultiPeak'), 'Download multiple peak figure')),
             theme = "light-border",
             placement = "right",
             arrow = FALSE
           )
         ),
         shiny::mainPanel(
           shiny::fluidRow(
             shiny::column(12, shiny::plotOutput(ns('plotMultiplePeaks'), click = ns('peak_multiple_click'), height='600px')),
           )
         )),
      shiny::h2('Visualise annotations'),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shinyWidgets::dropMenu(
            shinyWidgets::circleButton(ns("info_metadata"), icon = shiny::icon("info"),status = "success"),
            shiny::tags$div(
              shiny::tags$h3("Spatial peak visualisation"),
              shiny::tags$ul(
                shiny::tags$li("Select a metadata column and visualise it spatially across the samples selected at the top of this panel."),
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
            shiny::selectInput(ns("metadataName"), "Metadata to display:", multiple = FALSE, choices = colnames(full.metadata)[!(colnames(full.metadata)%in%c('x','y','Sample','x_tf','y_tf'))],selected=utils::tail(colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','Sample','x_tf','y_tf'))],1)),
            shiny::actionButton(ns("go_plot_metadata"),'Show spatial visualisation'),
          shiny::div(style = "margin-top:10px"),
          shinyWidgets::dropMenu(
            shinyWidgets::circleButton(ns("downloads"), icon = shiny::icon("download"),status = "success"),
            shiny::tags$div(
              shiny::tags$h3("Downloads"),
              shiny::fluidRow(
                shiny::column(10,offset=0,
                       shiny::tags$h4("Metadata distribution"),
                       shiny::textInput(ns('metaFileName'),'File name for download', value ='metadata.png', placeholder = 'metadata.png'),
                       shiny::numericInput(ns('metaWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       shiny::numericInput(ns('metaHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       shiny::downloadButton(ns('downloadMeta'), 'Download metadata figure'))
              )),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          )),
        shiny::mainPanel(
          shiny::plotOutput(ns('plotMetadata'),click = ns('metadata_click'),height='600px')),

      ))
  }else{
    NULL
  }
}

##' Spatial Visualisation Panel Server
##'
##' Provides the server logic for spatial visualisation of peaks and metadata in the SMEW app.
##'
##' @param id Shiny module id
##' @param bulk.metadata Data frame of bulk sample metadata
##' @param full.metadata Data frame of full sample metadata
##' @param full.intensity.matrix Matrix of intensities (features x samples)
##' @param anno Data frame of peak annotations
##' @return None; called for side effects in Shiny module
##' @export
IntroSpatialVisServer <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno){
  shiny::moduleServer(id, function(input, output, session){
    # remove constant peaks from list
    shiny::updateSelectizeInput(session, "peakName", choices = anno$display_name, server = TRUE, selected = anno$display_name[1])
    shiny::updateSelectizeInput(session, "peakName1", choices = anno$display_name, server = TRUE, selected = anno$display_name[1])
    shiny::updateSelectizeInput(session, "peakName2", choices = anno$display_name, server = TRUE, selected = anno$display_name[2])
    shiny::updateSelectizeInput(session, "peakName3", choices = anno$display_name, server = TRUE, selected = anno$display_name[3])

    show_metadata <- shiny::reactive({
      shiny::withProgress(message = 'Rendering metadata spatial plot...', value = 0, {
        current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
        current.metadata$metadata = current.metadata[,input[['metadataName']]]
        current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
        shiny::incProgress(0.7)
        ggplot2::ggplot(current.metadata,ggplot2::aes(x=.data$x_tf,y=.data$y_tf,color=.data$metadata,fill=.data$metadata))+ 
               ggplot2::geom_tile()+
               ggplot2::facet_wrap(~.data$Sample, nrow = max(1,floor(sqrt(length(input[['samplesToShow']]))/1.5)), scales = 'free')  +
               ggplot2::theme_classic() +
               ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                              axis.text.x=ggplot2::element_blank(),
                              axis.ticks.x=ggplot2::element_blank(),
                              axis.line.x = ggplot2::element_blank(),
                              axis.title.y=ggplot2::element_blank(),
                              axis.text.y=ggplot2::element_blank(),
                              axis.ticks.y=ggplot2::element_blank(),
                              axis.line.y = ggplot2::element_blank())+
               ggplot2::theme(aspect.ratio = 1)
      })
    }) |> shiny::bindEvent(input[["go_plot_metadata"]])

    show_metadata_zoom <- shiny::reactive({
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] == input$metadata_click$panelvar1,]
      current.metadata$metadata = current.metadata[,input[['metadataName']]]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      return(ggplot2::ggplot(current.metadata,ggplot2::aes(x=.data$x_tf,y=.data$y_tf,color=.data$metadata,fill=.data$metadata))+ 
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
    })

    show_peak <- shiny::reactive({
      shiny::withProgress(message = 'Rendering peak spatial plot...', value = 0, {
        my_peak = anno[anno$display_name == input[['peakName']],]
        current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
        current.intensity.matrix <- t(full.intensity.matrix)[,full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']]]
        caps = stats::quantile(current.intensity.matrix[my_peak$m_z,], probs = input[['capRange']]/100)
        current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
        current.metadata$peak = current.intensity.matrix[my_peak$m_z,]
        shiny::incProgress(0.4)
        if (input[['splitDensity']]){
          density.plot = ggplot2::ggplot(current.metadata, ggplot2::aes(x = .data$peak, color = .data$Sample)) +
            ggplot2::geom_density() + ggplot2::theme_classic() +
            ggplot2::geom_vline(xintercept = caps[1]) + ggplot2::geom_vline(xintercept = caps[2])
        } else {
          density.plot = ggplot2::ggplot(current.metadata, ggplot2::aes(x = .data$peak)) +
            ggplot2::geom_density() + ggplot2::theme_classic() +
            ggplot2::geom_vline(xintercept = caps[1]) + ggplot2::geom_vline(xintercept = caps[2])
        }
        current.metadata$peak = pmin(caps[2], current.metadata$peak)
        current.metadata$peak = pmax(caps[1], current.metadata$peak)
        if (input[['log2_intensity']]){
          current.metadata$peak = log2(current.metadata$peak + 1)
        }
        upper.lim = max(current.metadata$peak)
        legend_title = ifelse(input[['log2_intensity']],
                              paste0('log2 ', my_peak$m_z, '\n intensity'),
                              paste0(my_peak$m_z, '\n intensity'))
        spatial.plot = ggplot2::ggplot(current.metadata, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data$peak, fill = .data$peak)) +
          ggplot2::geom_tile() +
          ggplot2::facet_wrap(~ .data$Sample, scales = 'free',
                              nrow = max(1, floor(sqrt(length(input[['samplesToShow']])/1.5)))) +
          ggplot2::theme_classic() +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.line.x  = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank(),
                         axis.text.y  = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank(),
                         axis.line.y  = ggplot2::element_blank(),
                         legend.title = ggplot2::element_text(legend_title),
                         strip.text   = ggplot2::element_text(size = 10),
                         aspect.ratio = 1) +
          ggplot2::scale_fill_gradient(name = legend_title, low = "lightgrey", high = "brown", limits = c(0, upper.lim)) +
          ggplot2::scale_color_gradient(name = legend_title, low = "lightgrey", high = "brown", limits = c(0, upper.lim))
        shiny::incProgress(0.9)
        return(list('spatial' = spatial.plot,
                    'density' = density.plot,
                    'upper_lim' = upper.lim,
                    'full_df' = current.metadata))
      })
    }) |> shiny::bindEvent(input[["go_plot_peak"]])

    show_peak_zoom <- shiny::reactive({
      sp <- show_peak()
      shiny::req(sp$full_df)
      upper.lim <- sp$upper_lim
      my_peak <- anno[anno$display_name == input[['peakName']],]
      selected_sample <- input$peak_click$panelvar1
      df_zoom <- subset(sp$full_df, .data$Sample == selected_sample)

      legend_title <- ifelse(input[['log2_intensity']],
                             paste0('log2 ', my_peak$m_z, '\n intensity'),
                             paste0(my_peak$m_z, '\n intensity'))

      ggplot2::ggplot(df_zoom, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data$peak, fill = .data$peak)) +
        ggplot2::geom_tile() +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line.x  = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y  = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.line.y  = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(legend_title),
                       aspect.ratio = 1) +
        ggplot2::scale_fill_gradient(name = legend_title, low = "lightgrey", high = "brown", limits = c(0, upper.lim)) +
        ggplot2::scale_color_gradient(name = legend_title, low = "lightgrey", high = "brown", limits = c(0, upper.lim))
    })

    show_multiple_peaks <- shiny::reactive({
      shiny::withProgress(message = 'Rendering multi-peak spatial plot...', value = 0, {
        red_peak = anno[anno$display_name==input[['peakName1']],]
        green_peak = anno[anno$display_name==input[['peakName2']],]
        blue_peak = anno[anno$display_name==input[['peakName3']],]
        current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']],]
        current.intensity.matrix <- t(full.intensity.matrix)[,full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToShow']]]
        current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
        shiny::incProgress(0.3)
        red.caps = stats::quantile(current.intensity.matrix[red_peak$m_z,],probs=input[['capRange_multiple']]/100)
        current.metadata$red = current.intensity.matrix[red_peak$m_z,]
        current.metadata$red = pmin(red.caps[2],current.metadata$red)
        current.metadata$red = pmax(red.caps[1],current.metadata$red)
        green.caps = stats::quantile(current.intensity.matrix[green_peak$m_z,],probs=input[['capRange_multiple']]/100)
        current.metadata$green = current.intensity.matrix[green_peak$m_z,]
        current.metadata$green = pmin(green.caps[2],current.metadata$green)
        current.metadata$green = pmax(green.caps[1],current.metadata$green)
        blue.caps = stats::quantile(current.intensity.matrix[blue_peak$m_z,],probs=input[['capRange_multiple']]/100)
        current.metadata$blue = current.intensity.matrix[blue_peak$m_z,]
        current.metadata$blue = pmin(blue.caps[2],current.metadata$blue)
        current.metadata$blue = pmax(blue.caps[1],current.metadata$blue)
        if (input[['log2_intensity_multiple']]){
          current.metadata$red= log2(current.metadata$red+1)
          current.metadata$green= log2(current.metadata$green+1)
          current.metadata$blue= log2(current.metadata$blue+1)
        }
        caps = c(red.caps,green.caps,blue.caps)
        minima = c(min(current.metadata$red),min(current.metadata$green),min(current.metadata$blue))
        maxima = c(max(current.metadata$red),max(current.metadata$green),max(current.metadata$blue))
        current.metadata$red = (current.metadata$red-min(current.metadata$red))/(max(current.metadata$red)-min(current.metadata$red))
        current.metadata$green = (current.metadata$green-min(current.metadata$green))/(max(current.metadata$green)-min(current.metadata$green))
        current.metadata$blue = (current.metadata$blue-min(current.metadata$blue))/(max(current.metadata$blue)-min(current.metadata$blue))
        upper.lims = c(max(current.metadata$red),max(current.metadata$green),max(current.metadata$blue))
        current.metadata$rgb_hex <- grDevices::rgb(
          red = pmin(pmax(current.metadata$red, 0), 1),
          green = pmin(pmax(current.metadata$green, 0), 1),
          blue = pmin(pmax(current.metadata$blue, 0), 1)
        )
        shiny::incProgress(0.9)
        p <- ggplot2::ggplot(current.metadata, ggplot2::aes(x = .data$x_tf, y = .data$y_tf)) +
          ggplot2::geom_tile(ggplot2::aes(fill = .data$rgb_hex)) +
          ggplot2::scale_fill_identity() +
          ggplot2::facet_wrap(~ .data$Sample, scales = 'free',
                              nrow = max(1, floor(sqrt(length(input[['samplesToShow']]) / 1.5)))) +
          ggplot2::theme_classic() +
          ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 10),
            aspect.ratio = 1
          )
        list('plot' = p,'minima' = minima,'maxima' = maxima,'upper_lim' = upper.lims,'caps' = caps,'full_df' = current.metadata)
      })
    }) |> shiny::bindEvent(input[["go_plot_multiple_peaks"]])

    # Render main multi-peak RGB plot and compact legend
    output[['plotMultiplePeaks']] <- shiny::renderPlot({
      show_multiple_peaks()$plot
    }) 

    output[['multiLegend']] <- shiny::renderPlot({
      shiny::req(input[['peakName1']], input[['peakName2']], input[['peakName3']])
      r_title <- stringr::word(input[['peakName1']], sep = '_', start = 1, end = 2)
      g_title <- stringr::word(input[['peakName2']], sep = '_', start = 1, end = 2)
      b_title <- stringr::word(input[['peakName3']], sep = '_', start = 1, end = 2)

      vals <- seq(0, 1, length.out = 100)
      legend_df <- rbind(
        data.frame("channel" = r_title, "val" = vals, "fill_hex" = grDevices::rgb(vals, 0, 0)),
        data.frame("channel" = g_title, "val" = vals, "fill_hex" = grDevices::rgb(0, vals, 0)),
        data.frame("channel" = b_title, "val" = vals, "fill_hex" = grDevices::rgb(0, 0, vals))
      )
      legend_df$channel <- factor(legend_df$channel, levels = c(r_title, g_title, b_title))

      ggplot2::ggplot(legend_df, ggplot2::aes(x = .data$val, y = .data$channel, fill = .data$fill_hex)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_identity() +
        ggplot2::theme_classic() +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 9))
    }) |> shiny::bindEvent(input[["go_plot_multiple_peaks"]])

    

    show_multiple_peaks_zoom <- shiny::reactive({
      multiple_peaks <- show_multiple_peaks()
      shiny::req(multiple_peaks$full_df)
      sample_key <- input$peak_multiple_click$panelvar1
      df_zoom <- subset(multiple_peaks$full_df, .data$Sample == sample_key | .data$Group == sample_key)

      p <-
        ggplot2::ggplot(df_zoom, ggplot2::aes(x = .data$x_tf, y = .data$y_tf)) +
        ggplot2::geom_tile(ggplot2::aes(fill = .data$rgb_hex)) +
        ggplot2::scale_fill_identity() +
        ggplot2::theme_classic() +
        ggplot2::theme(
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank(),
          aspect.ratio = 1
        )
      return(p)
    })
    output[['plotMultiplePeakZoom']] <- shiny::renderPlot({
      show_multiple_peaks_zoom()})

    output[['downloadSpatial']] <- create_download_plot_handler(
      plot_func = function() show_peak()$spatial,
      filename_func = function() input[['spatialFileName']],
      width_func = function() input[['spatialWidth']],
      height_func = function() input[['spatialHeight']],
      units = 'in',
      dpi = 300
    )

    output[['downloadMultiPeak']] <- create_download_plot_handler(
      plot_func = function() show_multiple_peaks()$plot,
      filename_func = function() input[['multiPeakFileName']],
      width_func = function() input[['multiPeakWidth']],
      height_func = function() input[['multiPeakHeight']],
      units = 'in',
      dpi = 300
    )

    output[['peakDensity']] <- plotly::renderPlotly({
      plotly::ggplotly(show_peak()$density)
    })

    output[['downloadDensity']] <- create_download_plot_handler(
      plot_func = function() show_peak()$density,
      filename_func = function() input[['densityFileName']],
      width_func = function() input[['densityWidth']],
      height_func = function() input[['densityHeight']],
      units = 'in',
      dpi = 300
    )

    output[['plotPeak']] <- shiny::renderPlot({
      show_peak()$spatial}
      )

    output[['plotPeakZoom']] <- shiny::renderPlot({
      show_peak_zoom()})

    output[['plotMetadata']] <- shiny::renderPlot({
      show_metadata()
    }#, height = function() {plot_height()}
    )

    output[['plotMetadataZoom']] <- shiny::renderPlot({
      show_metadata_zoom()
    })

    output[['downloadMeta']] <- create_download_plot_handler(
      plot_func = function() show_metadata(),
      filename_func = function() input[['metaFileName']],
      width_func = function() input[['metaWidth']],
      height_func = function() input[['metaHeight']],
      units = 'in',
      dpi = 300
    )

    shiny::observeEvent(input$peak_click, {
      ns <- session$ns
      shiny::showModal(
        shiny::modalDialog(
          shiny::plotOutput(ns('plotPeakZoom')),
          easyClose = TRUE,
          footer = NULL
        )
      )
    })

    shiny::observeEvent(input$peak_multiple_click, {
      ns <- session$ns
      shiny::showModal(
        shiny::modalDialog(
          shiny::plotOutput(ns('plotMultiplePeakZoom')),
          easyClose = TRUE,
          footer = NULL
        )
      )
    })


    shiny::observeEvent(input$metadata_click, {
      ns <- session$ns
      shiny::showModal(
        shiny::modalDialog(
          shiny::plotOutput(ns('plotMetadataZoom')),
          easyClose = TRUE,
          footer = NULL
        )
      )
    })

  })
}
