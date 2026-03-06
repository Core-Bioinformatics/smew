#' Allows users to manually select regions of interest, optionally using histology images
#'
#' @description UI and server logic for the Manual Region Selection and Histology Overlay panel, enabling users to overlay metabolic data with histology images, draw custom regions of interest, and export region/cluster assignments in spatial omics data. Supports sample selection, overlay customisation, and download of results.
#'
#' @details
#' \itemize{
#'   \item{Overlay metabolic information with histology images (if available) and manually draw custom regions of interest.}
#'   \item{Customise cluster/region names.}
#'   \item{Select samples to display and interactively define regions.}
#'   \item{Export region or cluster assignments for downstream analysis.}
#'   \item{All results and overlays are available for download.}
#' }
#'
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame of bulk sample metadata (UI)
#' @param full.metadata Data frame of full sample metadata (UI, server)
#' @param full.intensity.matrix Matrix of intensities (features x samples, UI, server)
#' @param show Logical; whether to show the panel (default TRUE, UI)
#' @param anno Data frame of peak annotations (server)
#' @param shared_data Reactive or shared data object (server)
#' @name RegionPanel_HistologyTab
#' @rdname RegionPanel_HistologyTab
#' @return UI: A shiny tabPanel object for the manual region selection tab. Server: None (side effects in Shiny module).
#' @export
RegionPanel_HistologyTabUI <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, show = TRUE){
  ns <- shiny::NS(id)
  if(show){
    shiny::tabPanel('Manual Region Selection',
                    shiny::br(),
                    bslib::accordion(
                      bslib::accordion_panel(
                      title = "Information",
                      icon = bsicons::bs_icon("info-circle"),
                      shiny::tags$ul(
                        shiny::tags$li("This tab will overlay metabolic information with histology images (if available) and users can draw out custom regions of interest."),
                        shiny::tags$li("The appearance of overlay images can be customised, as can cluster names."),
                      )
                    ),
                    bslib::accordion_panel(
                      title = "Sample selection",
                      icon = bsicons::bs_icon("gear"),
                      shiny::selectInput(
                        inputId = ns("samplesToDisplay"),
                        label = "Select samples to display",
                        choices = unique(bulk.metadata[,1]),
                        selected = unique(bulk.metadata[,1]),
                        multiple = FALSE, width = '100%'
                      )
                    ),
                    bslib::accordion_panel(
                      title = "Downloads",
                      icon = bsicons::bs_icon("download"),
                      shiny::tags$div(
                        shiny::tags$h4("Download output plots"),
                        shiny::fluidRow(
                          shiny::column(6,
                                        shiny::tags$strong("Overlay plot"),
                                        shiny::textInput(ns('overlayFileName'), 'File name', value ='Overlay.png'),
                                        shiny::numericInput(ns('overlayWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                                        shiny::numericInput(ns('overlayHeight'),value = 8,label = 'Height (in)',min = 1,max = 50,step = 1),
                                        shiny::downloadButton(ns('downloadOverlay'), 'Download overlay')
                          ),
                          shiny::column(6,
                                        shiny::tags$strong("Selected plot"),
                                        shiny::textInput(ns('selectedFileName'), 'File name', value ='Selected.png'),
                                        shiny::numericInput(ns('selectedWidth'),value = 6,label = 'Width (in)',min = 1,max = 50,step = 1),
                                        shiny::numericInput(ns('selectedHeight'),value = 4,label = 'Height (in)',min = 1,max = 50,step = 1),
                                        shiny::downloadButton(ns('downloadSelected'), 'Download selected')
                          )
                        ),
                        shiny::fluidRow(
                          shiny::column(6,
                                        shiny::tags$strong("Overview plot"),
                                        shiny::textInput(ns('overviewFileName'), 'File name', value ='Overview.png'),
                                        shiny::numericInput(ns('overviewWidth'),value = 6,label = 'Width (in)',min = 1,max = 50,step = 1),
                                        shiny::numericInput(ns('overviewHeight'),value = 4,label = 'Height (in)',min = 1,max = 50,step = 1),
                                        shiny::downloadButton(ns('downloadOverview'), 'Download overview')
                          )
                        )
                      )
                    ),
                    id = ns("acc"),
                    open = "Sample selection"
                  ),
                  shiny::column(12,
                                shiny::fluidRow(
                                  shiny::selectizeInput(ns("feature"), "Select feature to color", choices = c()),
                                  shiny::textInput(ns("region_name"), "Region name", placeholder = "e.g. Tumor"),
                                  shiny::actionButton(ns("add_region"), "Add selection"),
                                  shiny::actionButton(ns("clear_regions"), "Clear regions"),
                                  shiny::sliderInput(ns("point_size"), "Spot size", min = 0.2, max = 3, value = 1, step = 0.1),
                                  shiny::sliderInput(ns("alpha"), "Transparency", min = 0.1, max = 1, value = 0.7, step = 0.1),
                                  plotly::plotlyOutput(ns("overlay"), height = "800px")
                                )
                  ),
                  shiny::fluidRow(
                    shiny::column(6,
                                  shiny::plotOutput(ns("selected"), height = "400px")
                    ),
                    shiny::column(6,
                                  shiny::plotOutput(ns("overview"), height = "400px")
                    )
                  ),
                  shiny::textInput(ns('clusterName'),label = 'Name for clusters'),
                  shiny::actionButton(ns('lock_cluster'),'Add clusters to object')
  )
  } else {
    NULL
  }
}

#' @rdname RegionPanel_HistologyTab
#' @export
RegionPanel_HistologyTabServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, anno, shared_data) {
  shiny::moduleServer(id, function(input, output, session){
    
    coords_df <- shiny::reactiveVal(NULL)
    all_coords_df <- shiny::reactiveVal(
      cbind(full.metadata, 'Region' = NA)
    )
    
    # Replace with your actual objects
    tissue_img <- shiny::reactive({
      if (file.exists(paste0('images/',input$samplesToDisplay,'/histology_aligned.jpg'))){
        img <- jpeg::readJPEG(paste0('images/',input$samplesToDisplay,'/histology_aligned.jpg'))
        return(img)
      }
    })
    
    static_coords = shiny::reactive({
      msi_coords = full.metadata[full.metadata$Sample == input$samplesToDisplay,]
      msi_coords['x'] = (msi_coords['x']-min(msi_coords['x']))/2
      msi_coords['y'] = (msi_coords['y']-min(msi_coords['y']))/2
      return(msi_coords)
    })
    
    shiny::observeEvent(input$samplesToDisplay, {
      msi_coords = static_coords()
      msi_coords$Region = NA
      coords_df(msi_coords)
    })
    
    coords_intensity <- shiny::reactive({
      msi_int = as.data.frame(full.intensity.matrix[full.metadata$Sample == input$samplesToDisplay,])
      return(msi_int)
    })
    
    # Feature selector
    shiny::observe({
      feats <- setdiff(colnames(full.intensity.matrix), c("pixel_id"))
      shiny::updateSelectizeInput(session = session, "feature", choices = c("None",feats), server = TRUE)
      shiny::updateSelectizeInput(session = session, "featureBox", choices = feats, server = TRUE)
    })
    
    overlay_plot <- function(){
      df <- static_coords()
      shiny::req(df)
      shiny::req(input$feature)
      if (file.exists(paste0('images/', input$samplesToDisplay, '/histology_aligned.jpg'))) {
        img <- tissue_img()
        if (input$feature == "None") {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = -.data$y)) +
            ggplot2::geom_point(size = input$point_size, alpha = 0) +
            ggplot2::annotation_raster(img, 0, ncol(img), -nrow(img), 0) +
            ggplot2::theme_void() + ggplot2::coord_fixed()
        } else {
          df[, input$feature] = coords_intensity()[, input$feature]
          p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = -.data$y,
                                                color = .data[[input$feature]])) +
            ggplot2::annotation_raster(img, 0, ncol(img), -nrow(img), 0) +
            ggplot2::geom_point(size = input$point_size, alpha = input$alpha) +
            ggplot2::theme_void() + ggplot2::coord_fixed() +
            ggplot2::scale_color_viridis_c()
        }
      } else {
        if (input$feature == "None") {
          p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = -.data$y)) +
            ggplot2::geom_point(size = input$point_size, alpha = 0) +
            ggplot2::theme_void() + ggplot2::coord_fixed()
        } else {
          df[,input$feature] = coords_intensity()[,input$feature]
          p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = -.data$y,
                                                color = .data[[input$feature]])) +
            ggplot2::geom_point(size = input$point_size, alpha = input$alpha) +
            ggplot2::theme_void() + ggplot2::coord_fixed() +
            ggplot2::scale_color_viridis_c()
        }    }
      p
    }
    # Overlay with lasso
    output$overlay <- plotly::renderPlotly({
      plotly::ggplotly(overlay_plot()) |> plotly::layout(dragmode = "lasso")
    })
    output[['downloadOverlay']] <- utils_create_download_plot_handler(
      plot_func = overlay_plot,
      filename_func = function() input[['overlayFileName']],
      width_func = function() input[['overlayWidth']],
      height_func = function() input[['overlayHeight']]
    )
    
    # For some reason this updating now doesn't work, need to find old verion
    # Add region when button clicked
    shiny::observeEvent(input$add_region, {
      sel <- plotly::event_data("plotly_selected")
      shiny::req(sel)
      
      df <- coords_df()
      all_clusters = all_coords_df()
      selected_rows <- sel$pointNumber + 1
      
      region_name <- trimws(input$region_name)
      if (region_name == "") {
        region_name <- paste0("Region_", length(unique(stats::na.omit(all_clusters$Region))) + 1)
      }
      df$Region[selected_rows] <- region_name
      coords_df(df)
      
      all_df = all_coords_df()
      idx <- match(all_df$pixel_id, df$pixel_id)
      to_update <- !is.na(idx)
      
      # Fill in cluster values where available
      all_df[['Region']][to_update] <- df$Region[idx[to_update]]
      all_coords_df(all_df)
    })
    
    # Clear regions
    shiny::observeEvent(input$clear_regions, {
      df <- coords_df()
      df$Region <- NA
      coords_df(df)
      
      all_df = all_coords_df()
      idx <- match(all_df$pixel_id, df$pixel_id)
      to_update <- !is.na(idx)
      
      # Fill in cluster values where available
      all_df[['Region']][to_update] <- NA
      all_coords_df(all_df)
    })
    
    selected_plot <- function(){
      df <- coords_df()
      if (file.exists(paste0('images/', input$samplesToDisplay, '/histology_aligned.jpg'))) {
        img <- tissue_img()
        ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = -.data$y)) +
          ggplot2::annotation_raster(img, 0, ncol(img), -nrow(img), 0) +
          ggrastr::geom_point_rast(ggplot2::aes(color = .data$Region),
                                   size = input$point_size,
                                   alpha = input$alpha) +
          ggplot2::scale_color_manual(values = stats::setNames(RColorBrewer::brewer.pal(8, "Set2"),
                                                               unique(stats::na.omit(df$Region))),
                                      na.translate = FALSE) +
          ggplot2::theme_void() + ggplot2::coord_fixed()
        
      } else {
        ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = -.data$y)) +
          ggrastr::geom_point_rast(ggplot2::aes(color = .data$Region),
                                   size = input$point_size,
                                   alpha = input$alpha) +
          ggplot2::scale_color_manual(values = stats::setNames(RColorBrewer::brewer.pal(8,"Set2"),
                                                               unique(stats::na.omit(df$Region))),
                                      na.translate = FALSE) +
          ggplot2::theme_void() + ggplot2::coord_fixed()
      }
    }
    # Show updated selection
    output$selected <- shiny::renderPlot({ selected_plot() })
    output[['downloadSelected']] <- utils_create_download_plot_handler(
      plot_func = selected_plot,
      filename_func = function() input[['selectedFileName']],
      width_func = function() input[['selectedWidth']],
      height_func = function() input[['selectedHeight']]
    )
    
    
    overview_plot <- function(){
      df <- all_coords_df()
      if (all(is.na(df$Region))){
        ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = - .data$y_tf)) + 
          ggplot2::geom_tile(size = input$point_size,fill='gray') +
          ggplot2::theme_void() + 
          ggplot2::facet_wrap(ggplot2::vars(.data$Sample),
            nrow = max(1, floor(sqrt(length(unique(df$Sample)) / 1.5)))) + 
          ggplot2::coord_fixed()
      } else {
        ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = - .data$y_tf, fill = .data$Region)) + 
          ggplot2::geom_tile(size = input$point_size) +
          ggplot2::theme_void() + 
          ggplot2::facet_wrap(ggplot2::vars(.data$Sample), nrow = max(1, floor(sqrt(length(unique(df$Sample)) / 1.5)))) + 
          ggplot2::scale_fill_discrete(na.value = "gray") +
          ggplot2::coord_fixed()
      }
    }
    output$overview <- shiny::renderPlot({ overview_plot() })
    output[['downloadOverview']] <- utils_create_download_plot_handler(
      plot_func = overview_plot,
      filename_func = function() input[['overviewFileName']],
      width_func = function() input[['overviewWidth']],
      height_func = function() input[['overviewHeight']]
    )
    
    shiny::observeEvent(input$lock_cluster, {
      meta <- shared_data$updated.metadata
      clusterdf = all_coords_df()
      # meta[input$clusterName] <- clusterdf$cluster[
      #   match(meta$pixel_id, clusterdf$pixel_id)
      # ]
      if (!input$clusterName %in% colnames(meta)) {
        meta[[input$clusterName]] <- NA
      }
      
      # Update only rows that appear in metadata_sub
      idx <- match(meta$pixel_id, clusterdf$pixel_id)
      to_update <- !is.na(idx)
      
      # Fill in cluster values where available
      meta[[input$clusterName]][to_update] <- clusterdf$Region[idx[to_update]]
      
      shared_data$updated.metadata <- shiny::isolate(meta)   # update shared dataset
    })
    
  }
  )
}
