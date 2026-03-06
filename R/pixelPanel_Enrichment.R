
#' Performs and visualises spatial pathway enrichment analysis
#'
#' @description UI and server logic for spatial pathway enrichment analysis at the pixel level in the SMEW app. Visualizes spatially-resolved pathway enrichment p-values, distributions, and significance across samples.
#'
#' @details
#' \itemize{
#'  \item{Visualize pre-computed spatial pathway enrichment p-values across samples.}
#'  \item{Show adjusted p-values or -log10(p-values) spatially.}
#'  \item{Density and bar plots for p-value distributions and proportion of significant pixels.}
#'  \item{Download handlers for spatial and bar plots.}
#'  \item{Interactive controls for pathway selection, capping, and significance threshold.}
#' }
#' @param id Shiny module id
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param full.metadata Data frame of full sample metadata
#' @param pixel.enrichment List or data structure with pixel-level enrichment results
#' @param show Logical; whether to show the panel (default TRUE, UI only)
#' @param full.intensity.matrix Matrix of intensities (features x samples, server only)
#' @param anno Data frame of peak annotations (server only)
#'
#' @return UI: A shiny::tabPanel object for the enrichment panel. Server: None; called for side effects in Shiny module.
#' @export
#' @name PixelPanel_EnrichmentTab
#' @rdname PixelPanel_EnrichmentTab
PixelPanel_EnrichmentTabUI <- function(id, bulk.metadata, full.metadata, pixel.enrichment, show = TRUE) {
  ns <- shiny::NS(id)
  if (show) {
    shiny::tabPanel(
      'Spatial Enrichment',
      bslib::accordion(
        bslib::accordion_panel(
          title = "Information",
          icon = bsicons::bs_icon("info-circle"),
          shiny::tags$ul(
            shiny::tags$li("Select a pathway and visualise its pre-computed spatial p-values across the samples selected at the top of this panel."),
            shiny::tags$li("The adjusted p-value can be shown directly or the -log10 p-value"),
            shiny::tags$li("Density plots and barplots showing the distribution of p-values and proportion of significant pixels are also shown"),
            shiny::tags$li("Caps can also be applied to the colour scale to aid visualisation.")
          )
        ),
        bslib::accordion_panel(
          title = "Sample selection",
          icon = bsicons::bs_icon("gear"),
          shiny::selectInput(
            inputId = ns("samplesToShow"),
            label = "Select samples to show",
            choices = c(),
            multiple = TRUE
          )
        ),
        id = ns("acc"),
        open = "Sample selection"
      ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(
            style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_peak"), icon = shiny::icon("info"), status = "info"),
              shiny::tags$div(
                shiny::tags$h3("Spatial enrichment visualisation"),
                shiny::tags$ul(
                  shiny::tags$li("Select a pathway and visualise its p-value across the samples selected at the top of this panel."),
                  shiny::tags$li("The adjusted p-value can be shown directly or the -log10 p-value"),
                  shiny::tags$li("Density plots and barplots showing the distribution of p-values and proportion of significant pixels are also shown"),
                  shiny::tags$li("Caps can also be applied to the colour scale to aid visualisation.")
                )
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloadsMetadata"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Downloads"),
                shiny::fluidRow(
                  shiny::column(5, offset = 0,
                    shiny::tags$h4("Spatial distribution"),
                    shiny::textInput(ns('spatialFileName'), 'File name for download', value = 'spatial.png', placeholder = 'spatial.png'),
                    shiny::numericInput(ns('spatialWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                    shiny::numericInput(ns('spatialHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                    shiny::downloadButton(ns('downloadSpatial'), 'Download spatial figure')
                  ),
                  shiny::column(5, offset = 1,
                    shiny::tags$h4("Proportion of significant pixels plot"),
                    shiny::textInput(ns('barFileName'), 'File name for download', value = 'bar.png', placeholder = 'bar.png'),
                    shiny::numericInput(ns('barWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                    shiny::numericInput(ns('barHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                    shiny::downloadButton(ns('downloadBar'), 'Download bar plot')
                  )
                )
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::actionButton(ns("fix_samples"), 'Fix samples and sort pathways by significance'),
          shiny::selectInput(ns("pathwayName"), "Pathway to show:", multiple = FALSE, choices = c()),
          shiny::selectInput(ns("pvalShown"), "How to show p-value:", choices = c('Adjusted p-value', '-log10(adjusted p-value)')),
          shiny::numericInput(ns('sigThreshold'), 'Significance threshold', value = 0.05, min = 0.001, max = 1, step = 0.001),
          shiny::conditionalPanel(
            id = ns('pvalCap_log'),
            ns = ns,
            condition = "input.pvalShown == '-log10(adjusted p-value)'",
            shiny::numericInput(ns('logCap'), 'Upper cap on log scale', value = 10, min = 1, max = 1000, step = 1)
          ),
          shiny::radioButtons(ns('barplot_metadata'), label = "Color bar plot by",
            choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[ncol(bulk.metadata)])
        ),
        shiny::mainPanel(
          shiny::tags$h3('Spatial enrichment across samples'),
          shiny::fluidRow(shiny::column(10, shiny::plotOutput(ns('plotPeak'), click = ns('peak_click'), height = "600px"))),
          shiny::tags$h3('Proportion of significant pixels'),
          shiny::plotOutput(ns('plotBar')),
        )
      )
    )
  } else {
    NULL
  }
}

#' @rdname PixelPanel_EnrichmentTab
#' @param full.intensity.matrix Matrix of intensities (features x samples)
#' @param anno Data frame of peak annotations
#' @return None; called for side effects in Shiny module
#' @export
PixelPanel_EnrichmentTabServer <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno, pixel.enrichment){

  shiny::moduleServer(id, function(input, output, session){

    shiny::updateSelectInput(
      session,
      "samplesToShow",
      choices = names(pixel.enrichment),
      )
    ordered_pathway_list <- shiny::reactive({
      shiny::req(input[['samplesToShow']])
      shiny::withProgress(message = 'Sorting pathways by significance...', value = 0, {
        shiny::incProgress(0.3)
        pathway_table = pixel.enrichment[input[['samplesToShow']]]
        pathway_table = dplyr::bind_rows(pathway_table, .id = "sample")
        pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
        pathway_table = pathway_table |> dplyr::group_by(.data$pathway) |> dplyr::summarise(non_na_count = sum(!is.na(.data$pathwaySig)), .groups = 'drop') |>
          dplyr::arrange(dplyr::desc(.data$non_na_count), .by_group = FALSE) |> dplyr::filter(.data$non_na_count!=0) |>
          dplyr::select(.data$pathway)
        shiny::incProgress(0.7)
        return(pathway_table$pathway)
      })
    })  |>
      shiny::bindEvent(input[['fix_samples']])

    shiny::observe({
    shiny::updateSelectizeInput(session, "pathwayName", choices = ordered_pathway_list(), server = TRUE, selected = ordered_pathway_list()[1])
    })

    selected_pathway_table <- shiny::reactive({
      shiny::req(input[['samplesToShow']])
      shiny::withProgress(message = 'Preparing pathway table...', value = 0, {
        shiny::incProgress(0.4)
        pathway_table = pixel.enrichment[input[['samplesToShow']]]
        pathway_table = dplyr::bind_rows(pathway_table, .id = "sample")
        shiny::incProgress(0.6)
        return(pathway_table)
      })
    }) |>
      shiny::bindEvent(input[['fix_samples']])
    selected.samples.metadata <- shiny::reactive({
      return(full.metadata |> dplyr::filter(.data$Sample %in% input[['samplesToShow']]))
    }) |>
      shiny::bindEvent(input[['fix_samples']])
    show_peak <- shiny::reactive({
      shiny::req(input[['pathwayName']])
      pathway_table = selected_pathway_table() |> dplyr::filter(.data$pathway==input[['pathwayName']])
      pathway_table = merge(selected.samples.metadata(),pathway_table,all.x=T)
      pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
      if (input[['pvalShown']]=='-log10(adjusted p-value)'){
        pathway_table$pathwaySig = -log10(pathway_table$pathwaySig)
        pathway_table$pathwaySig = ifelse(pathway_table$pathwaySig>input[['logCap']],input[['logCap']],pathway_table$pathwaySig)
      }
      legend_title=input[['pvalShown']]
      spatial.plot = ggplot2::ggplot(pathway_table,ggplot2::aes(x=.data$x_tf,y=.data$y_tf,color=.data$pathwaySig,fill=.data$pathwaySig)) +
      ggplot2::geom_tile()+
        ggplot2::facet_wrap(~.data$Sample,nrow=max(1,floor(sqrt(length(input[['samplesToShow']])/1.5))))  +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(legend_title),
                       aspect.ratio = 1)
      if (input[['pvalShown']]=='-log10(adjusted p-value)'){
        spatial.plot <- spatial.plot +
          ggplot2::scale_fill_gradient(name=legend_title,low = "lightgrey", high = "brown",na.value = 'lightgrey')+
          ggplot2::scale_color_gradient(name=legend_title,low = "lightgrey", high = "brown",na.value = 'lightgrey')
      } else {
        spatial.plot <- spatial.plot +
        ggplot2::scale_fill_gradient(name=legend_title,high = "lightgrey", low = "brown",na.value = 'lightgrey')+
        ggplot2::scale_color_gradient(name=legend_title,high = "lightgrey", low = "brown",na.value = 'lightgrey')}

         density.plot = ggplot2::ggplot(pathway_table,ggplot2::aes(x=.data$FDR,color=.data$Sample))+ggplot2::geom_density()+ggplot2::theme_classic()
      return(list('spatial'=spatial.plot,'density'=density.plot))
    }) |>
      shiny::bindCache(input[['pathwayName']],input[['pvalShown']],input[['logCap']],input[['sigThreshold']],input[['fix_samples']])

    prop_significant <- shiny::reactive({
      shiny::req(input[['pathwayName']])
      if (!is.null(input[['pathwayName']])){
      pathway_table = selected_pathway_table() |> dplyr::filter(.data$pathway==input[['pathwayName']])
      pathway_table = merge(selected.samples.metadata(),pathway_table,all.x=T)
      pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
      pathway_table = pathway_table |> dplyr::group_by(dplyr::across(dplyr::all_of(c('Sample',input[['barplot_metadata']])))) |> dplyr::summarise("non_na_count" = sum(!is.na(.data$pathwaySig)),count = dplyr::n(), .groups = 'drop') |>
        dplyr::arrange(dplyr::desc(.data$non_na_count), .by_group = FALSE) |> dplyr::filter(.data$non_na_count!=0)
      ggplot2::ggplot(pathway_table,ggplot2::aes(x=.data$Sample,y=.data$non_na_count/.data$count,fill=.data[[input[['barplot_metadata']]]]))+
        ggplot2::geom_bar(stat='identity')+
        ggplot2::ylab('Proportion of significant pixels')+
        ggplot2::xlab('Sample')+
        ggplot2::theme_classic()+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))+
        ggplot2::scale_fill_discrete(name=input[['barplot_metadata']])
      }
    })

    show_peak_zoom <- shiny::reactive({
      shiny::req(input[['pathwayName']], input$peak_click$panelvar1)
      pathway_table = pixel.enrichment[[input$peak_click$panelvar1]] |> dplyr::filter(.data$pathway==input[['pathwayName']])
      pathway_table = merge(full.metadata |> dplyr::filter(.data$Sample == input$peak_click$panelvar1),pathway_table,all.x=T)

      pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
      if (input[['pvalShown']]=='-log10(adjusted p-value)'){
        pathway_table$pathwaySig = -log10(pathway_table$pathwaySig)
        pathway_table$pathwaySig = ifelse(pathway_table$pathwaySig>input[['logCap']],input[['logCap']],pathway_table$pathwaySig)
      }
      legend_title=input[['pvalShown']]
      spatial.plot = ggplot2::ggplot(pathway_table,ggplot2::aes(x=.data$x,y=.data$y,color=.data$pathwaySig,fill=.data$pathwaySig))+ggplot2::geom_tile()+
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(legend_title),
                       aspect.ratio = 1)
      if (input[['pvalShown']]=='-log10(adjusted p-value)'){
        spatial.plot <- spatial.plot +
          ggplot2::scale_fill_gradient(name=legend_title,low = "lightgrey", high = "brown",na.value = 'lightgrey')+
          ggplot2::scale_color_gradient(name=legend_title,low = "lightgrey", high = "brown",na.value = 'lightgrey')
      } else {
        spatial.plot <- spatial.plot +
          ggplot2::scale_fill_gradient(name=legend_title,high = "lightgrey", low = "brown",na.value = 'lightgrey')+
          ggplot2::scale_color_gradient(name=legend_title,high = "lightgrey", low = "brown",na.value = 'lightgrey')}
     return(spatial.plot)
    })


    output[['plotBar']] <- shiny::renderPlot({
      prop_significant()})

    output[['plotPeak']] <- shiny::renderPlot({
      show_peak()$spatial})


    output[['plotPeakZoom']] <- shiny::renderPlot({
      show_peak_zoom()})

    output[['peakDensity']] <- shiny::renderPlot({
      show_peak()$density
    })

    output[['downloadSpatial']] <- utils_create_download_plot_handler(
      plot_func = function() show_peak()$spatial,
      filename_func = function() input[['spatialFileName']],
      width_func = function() input[['spatialWidth']],
      height_func = function() input[['spatialHeight']]
    )

    output[['downloadBar']] <- utils_create_download_plot_handler(
      plot_func = prop_significant,
      filename_func = function() input[['barFileName']],
      width_func = function() input[['barWidth']],
      height_func = function() input[['barHeight']]
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


  })
}
