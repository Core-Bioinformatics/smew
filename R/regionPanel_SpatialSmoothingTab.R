#' Spatial smoothing helpers
#' @description Functions for smoothing cluster labels across spatial neighborhoods (used by the Cluster Smoothing tab).
#' @keywords internal

smoothed.cluster <- function(metadata.sub,cluster.col){
  sample = base::unique(metadata.sub$Sample)
  metadata.sub[,cluster.col] = base::as.character(metadata.sub[,cluster.col])
  base::rownames(metadata.sub)=metadata.sub$pixel_id
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('metadata.sub','select.label'),envir=base::environment())
    smoothed.clusters = pbapply::pbsapply(base::rownames(metadata.sub),FUN=function(x)select.label(x,metadata.sub,cluster.col),simplify = TRUE, cl=cl)
  parallel::stopCluster(cl)
  base::names(smoothed.clusters)=metadata.sub$pixel_id
  return(smoothed.clusters)
}


select.label = function(i, current.metadata, cluster.col){
  current.col = current.metadata[i,'x_tf']
  current.row = current.metadata[i,'y_tf']
    selected.cluster = dplyr::select(dplyr::filter(current.metadata, .data$x_tf == current.col & .data$y_tf == current.row), dplyr::one_of(cluster.col))
  neighbour.clusters = dplyr::select(
    dplyr::filter(
      dplyr::filter(
        dplyr::filter(current.metadata, .data$x_tf %in% base::seq(current.col-1,current.col+1)),
        .data$y_tf %in% base::seq(current.row-1,current.row+1)),
      .data$x_tf!=current.col | .data$y_tf!=current.row),
    dplyr::one_of(cluster.col))
  if (base::nrow(neighbour.clusters)<5){
    return(selected.cluster[,cluster.col])
  } else {
    top.cluster = base::sort(base::table(neighbour.clusters[,cluster.col]),decreasing = TRUE)[1]
    if (top.cluster > (0.5 * base::nrow(neighbour.clusters))){
      return(base::names(top.cluster))
    } else {
      return(selected.cluster[,cluster.col])
    }
  }
}

# Jaccard similarity function
jaccard_index <- function(x, y) {
  base::mean(x == y, na.rm = TRUE)
}


#' RegionSmoothingUI
#'
#' UI for the region smoothing tab.
#'
#' @param id Shiny module id
#' @return A shiny tabPanel object for the region smoothing tab
#' @export
RegionSmoothingUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel(
    'Region Smoothing',
    bslib::accordion(
      bslib::accordion_panel(
        title = "Information",
        icon = bsicons::bs_icon("arrow-right-circle"),
        shiny::tags$ul(
          shiny::tags$li("Apply spatial smoothing to clusters and regions to reduce noise and enhance spatial patterns."),
          shiny::tags$li("Smoothed clusters can then be used by the following tabs.")
        )
      ),
      bslib::accordion_panel(
        title = "Downloads",
        icon = bsicons::bs_icon("download"),
        shiny::fluidRow(
          shiny::column(6,
            shiny::tags$h5("Download non-smoothed plot"),
            shiny::textInput(ns('originalFileName'), 'File name for original plot', value ='OriginalClusters.png'),
            shiny::numericInput(ns('originalWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
            shiny::numericInput(ns('originalHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
            shiny::downloadButton(ns('downloadOriginal'), 'Download original')
          ),
          shiny::column(6,
            shiny::tags$h5("Download smoothed plot"),
            shiny::textInput(ns('smoothedFileName'), 'File name for smoothed plot', value ='SmoothedClusters.png'),
            shiny::numericInput(ns('smoothedWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
            shiny::numericInput(ns('smoothedHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
            shiny::downloadButton(ns('downloadSmoothed'), 'Download smoothed')
          )
        )
      ),

      open = NULL
    ),
    bslib::card(
      class = "mb-3",
      style = "width:100%;",
      shiny::fluidRow(
        shiny::column(4,
          shiny::actionButton(ns('update_options'),'Update clustering options')
        ),
        shiny::column(4,
          shiny::selectInput(ns("col1"), "Select cluster column:", choices = c())
        ),
        shiny::column(4,
          shiny::actionButton(ns('run_smoothing'), label ='Run spatial smoothing')
        )
      )
    ),
    shiny::tags$h4("Non-smoothed regions"),
    shiny::plotOutput(ns("originalPlot"),height = '600px'),
    shiny::tags$h4("Smoothed regions"),
    shiny::plotOutput(ns("smoothedPlot"),height = '600px')
  )
}

#' RegionSmoothingServer
#'
#' Server logic for the region smoothing tab.
#'
#' @param id Shiny module id
#' @param shared_data Reactive or shared data object
#' @param full.metadata Data frame of full sample metadata
#' @return None; called for side effects in Shiny module
#' @export
RegionSmoothingServer <- function(id, shared_data, full.metadata) {
  shiny::moduleServer(id, function(input, output, session) {
    
    shiny::observe({
      categorical_cols <- colnames(full.metadata)[sapply(full.metadata, function(col) {
        n_unique <- length(unique(stats::na.omit(col)))
        n_unique >= 1 && n_unique < 20  # exclude constants & numeric-like columns
      })]
      
      shiny::updateSelectInput(session,'col1',choices = categorical_cols)
    })
    
    shiny::observeEvent(input$update_options, {
      categorical_cols <- colnames(shared_data$updated.metadata)[sapply(shared_data$updated.metadata, function(col) {
        n_unique <- length(unique(stats::na.omit(col)))
        n_unique >= 1 && n_unique < 20  # exclude constants & numeric-like columns
      })]
      
      shiny::updateSelectInput(session,'col1',choices = categorical_cols)
    })
    
    output$originalPlot <- shiny::renderPlot({
      df <- shared_data$updated.metadata
      shown_samples = unique(df[!is.na(df[,input$col1]),]$Sample)
      df = df[df$Sample %in% shown_samples,]
      if (is.numeric(df[,input$col1])){
        df[,input$col1] = factor(df[,input$col1])
      }
      ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, fill = .data[[input$col1]], color = .data[[input$col1]])) +
        ggplot2::geom_tile() +
        ggplot2::theme_classic() +
        ggplot2::facet_wrap(~.data$Sample,
                   nrow = max(1,floor(sqrt(length(unique(df$Sample))/1.5))),
                   scales='free')+
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
    output[['downloadOriginal']] <- create_download_plot_handler(
      plot_func = function(){
        df <- shared_data$updated.metadata
        shown_samples = unique(df[!is.na(df[,input$col1]),]$Sample)
        df = df[df$Sample %in% shown_samples,]
        if (is.numeric(df[,input$col1])){
          df[,input$col1] = factor(df[,input$col1])
        }
        ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, fill = .data[[input$col1]], color = .data[[input$col1]])) +
          ggplot2::geom_tile() +
          ggplot2::theme_classic() +
          ggplot2::facet_wrap(~.data$Sample,
                     nrow = max(1,floor(sqrt(length(unique(df$Sample))/1.5))),
                     scales='free')+
          ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                         axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank(),
                         axis.line.x = ggplot2::element_blank(),
                         axis.title.y=ggplot2::element_blank(),
                         axis.text.y=ggplot2::element_blank(),
                         axis.ticks.y=ggplot2::element_blank(),
                         axis.line.y = ggplot2::element_blank())+
          ggplot2::theme(aspect.ratio = 1)
      },
      filename_func = function() input[['originalFileName']],
      width_func = function() input[['originalWidth']],
      height_func = function() input[['originalHeight']]
    )
    
    smoothed_clusters <- shiny::reactive({
      shiny::withProgress(message = 'Running spatial smoothing...', value = 0, {
        current.metadata = shared_data$updated.metadata
        current.metadata = current.metadata[!is.na(current.metadata[,input$col1]),]
        smoothed.clusters = c()
        samples <- unique(current.metadata$Sample)
        for (i in seq_along(samples)){
          my.sample <- samples[i]
          smoothed.clusters = c(smoothed.clusters, smoothed.cluster(current.metadata[current.metadata$Sample==my.sample,], input$col1))
          shiny::incProgress(i/length(samples))
        }
        current.metadata[,paste0(input$col1,'_smoothed')] = smoothed.clusters
        meta <- shared_data$updated.metadata
        if (!paste0(input$col1,'_smoothed') %in% colnames(meta)) {
          meta[[paste0(input$col1,'_smoothed')]] <- NA
        }
        # Update only rows that appear in metadata_sub
        idx <- match(meta$pixel_id, current.metadata$pixel_id)
        to_update <- !is.na(idx)
        # Fill in cluster values where available
        meta[[paste0(input$col1,'_smoothed')]][to_update] <- current.metadata[,paste0(input$col1,'_smoothed')][idx[to_update]]
        shared_data$updated.metadata <- shiny::isolate(meta)   # update shared dataset
        return(current.metadata)
      })
    }) |> shiny::bindEvent(input$run_smoothing)
    
    output$smoothedPlot <- shiny::renderPlot({
      
      df <- smoothed_clusters()
      shown_samples = unique(df[!is.na(df[,input$col1]),]$Sample)
      df = df[df$Sample %in% shown_samples,]
      
      ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, fill = .data[[paste0(input$col1,'_smoothed')]], color = .data[[paste0(input$col1,'_smoothed')]])) +
        ggplot2::geom_tile() +
        ggplot2::theme_classic() +
        ggplot2::facet_wrap(~.data$Sample,
                   nrow = max(1,floor(sqrt(length(unique(df$Sample))/1.5))),
                   scales='free')+
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
    output[['downloadSmoothed']] <- create_download_plot_handler(
      plot_func = function(){
        df <- smoothed_clusters()
        shown_samples = unique(df[!is.na(df[,input$col1]),]$Sample)
        df = df[df$Sample %in% shown_samples,]
        ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, fill = .data[[paste0(input$col1,'_smoothed')]], color = .data[[paste0(input$col1,'_smoothed')]])) +
          ggplot2::geom_tile() +
          ggplot2::theme_classic() +
          ggplot2::facet_wrap(~.data$Sample,
                     nrow = max(1,floor(sqrt(length(unique(df$Sample))/1.5))),
                     scales='free')+
          ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                         axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank(),
                         axis.line.x = ggplot2::element_blank(),
                         axis.title.y=ggplot2::element_blank(),
                         axis.text.y=ggplot2::element_blank(),
                         axis.ticks.y=ggplot2::element_blank(),
                         axis.line.y = ggplot2::element_blank())+
          ggplot2::theme(aspect.ratio = 1)
      },
      filename_func = function() input[['smoothedFileName']],
      width_func = function() input[['smoothedWidth']],
      height_func = function() input[['smoothedHeight']]
    )
    
    
  })
}