#' Analyses and visualises relationships between molecular intensity and distance from a region of interest
#'
#' @description UI and server logic for the ROI Distance Correlation panel, enabling users to correlate molecular intensity with distance from a selected region of interest (ROI) in spatial omics data. Supports region selection, correlation analysis, visualisation of results (volcano plots, spatial maps, tables), and download of results.
#'
#' @details
#' \itemize{
#'   \item{Correlate molecular intensity with distance from a selected region of interest (ROI) for each sample.}
#'   \item{Choose a region column and label to define the ROI; distance from each pixel to the ROI is calculated.}
#'   \item{Select the correlation method (Pearson or Spearman) to assess the relationship between distance and feature intensity.}
#'   \item{View volcano plots of correlation results, spatial maps of distances, and tables of feature correlations.}
#'   \item{Visualise smoothed intensity profiles as a function of distance from the ROI for selected features, optionally split by metadata.}
#'   \item{All results and plots can be downloaded for further analysis or reporting.}
#' }
#'
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame of bulk sample metadata (UI)
#' @param show Logical; whether to show the panel (default TRUE, UI)
#' @param full.intensity.matrix Matrix of intensities (features x samples, server)
#' @param full.metadata Data frame of full sample metadata (server)
#' @param shared_data Reactive or shared data object (server)
#' @param anno Data frame of peak annotations (server)
#' @return UI: A shiny tabPanel object for the radial distance correlation tab. Server: None (side effects in Shiny module).
#' @name RegionPanel_RadialDistanceTab
#' @rdname RegionPanel_RadialDistanceTab
#' @export
RegionPanel_RadialDistanceTabUI <- function(id, bulk.metadata, show = TRUE) {
  ns <- shiny::NS(id)
  if (show){
  shiny::tabPanel(
    "Radial Distance",
    bslib::accordion(
      bslib::accordion_panel(
        title = "Information",
        icon = bsicons::bs_icon("arrow-right-circle"),
        shiny::tags$ul(
          shiny::tags$li("This tab allows you to correlate molecular intensity with distance from a selected region of interest (ROI)."),
          shiny::tags$li("Choose a region column and label to define the ROI. The distance from each pixel to the ROI will be calculated for each sample."),
          shiny::tags$li("You can select the correlation method (Pearson or Spearman) to assess the relationship between distance and feature intensity."),
          shiny::tags$li("After running the correlation, you can view a volcano plot of correlation results, a spatial map of distances, and a table of feature correlations."),
          shiny::tags$li("The 'Distance vs intensity trends' tab allows you to visualize smoothed intensity profiles as a function of distance from the ROI for selected features, optionally split by metadata."),
          shiny::tags$li("All results and plots can be downloaded for further analysis or reporting.")
        )
      ),
      bslib::accordion_panel(
        title = "Sample selection",
        icon = bsicons::bs_icon("gear"),
        shiny::selectInput(ns("selected_samples"), "Select samples for analysis:",
          choices = unique(bulk.metadata$Sample),
          selected = unique(bulk.metadata$Sample),
          multiple = TRUE)
      ),
      bslib::accordion_panel(
        title = "Downloads",
        icon = bsicons::bs_icon("download"),
        shiny::tags$div(
          shiny::tags$h4("Download output plots and tables"),
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$strong("Spatial distance plot"),
              shiny::textInput(ns('spatialCorrFileName'), 'File name', value ='SpatialDistance.png'),
              shiny::numericInput(ns('spatialCorrWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('spatialCorrHeight'),value = 4,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadSpatialCorr'), 'Download spatial distance')
            ),
            shiny::column(6,
              shiny::tags$strong("Volcano plot"),
              shiny::textInput(ns('volcanoCorrFileName'), 'File name', value ='CorrelationVolcano.png'),
              shiny::numericInput(ns('volcanoCorrWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('volcanoCorrHeight'),value = 4,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadVolcanoCorr'), 'Download volcano')
            )
          ),
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$strong("Distance vs intensity trends"),
              shiny::textInput(ns('smoothFileName'), 'File name', value ='DistanceIntensity.png'),
              shiny::numericInput(ns('smoothWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('smoothHeight'),value = 5,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadSmooth'), 'Download distance trends')
            ),
            shiny::column(6,
              shiny::tags$strong("Correlation results table"),
              shiny::br(),
              shiny::downloadButton(ns("download_corr"), "Download correlation results")
            )
          )
        )
      ),
      id = ns("acc"),
      open = FALSE
    ),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput(ns("roi_column"), "ROI column:", choices = NULL),
        shiny::selectInput(ns("roi_label"), "ROI label (region of interest):", choices = NULL),
        shiny::selectInput(ns("corr_method"), "Correlation method:", choices = c("pearson", "spearman")),
        shiny::actionButton(ns("run_corr"), "Run correlation"),
        shiny::hr(),
        shiny::uiOutput(ns("feature_select_ui")),
        shiny::checkboxInput(ns('split_smoothed'),'Should smoothed intensity vs distance plot be split by sample-wide metadata?',value = F),
        shiny::conditionalPanel(condition = "input.split_smoothed",
          shiny::column(6,shiny::selectInput(ns('split_var'),label = 'Metadata for split',choices = base::colnames(bulk.metadata),selected = base::colnames(bulk.metadata)[1])),
          ns=ns)
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel("Correlation summary", 
            shiny::plotOutput(ns("spatial_plot"), height = '600px'),
            plotly::plotlyOutput(ns("volcano_plot"), height = 400),
            DT::dataTableOutput(ns("corr_table"))
          ),
          shiny::tabPanel("Distance vs intensity trends", 
            shiny::plotOutput(ns("smooth_plot"), height = 500)
          )
        )
      )
    )
  )
  }
}

#' @rdname RegionPanel_RadialDistanceTab
#' @export
RegionPanel_RadialDistanceTabServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, shared_data, anno) {
  shiny::moduleServer(id, function(input, output, session) {
    shiny::observe({
      categorical_cols <- colnames(shared_data$updated.metadata)[sapply(shared_data$updated.metadata, function(col) {
        n_unique <- length(unique(stats::na.omit(col)))
        n_unique >= 1 && n_unique < 20  # exclude constants & numeric-like columns
      })]
      sub_sample_cols = categorical_cols[!(categorical_cols %in% colnames(bulk.metadata))]
      shiny::updateSelectInput(session, "roi_column", choices = sub_sample_cols)
    })
    
    shiny::observe({
      shiny::req(input$roi_column)
      shiny::updateSelectInput(session, "roi_label", 
                        choices = unique(shared_data$updated.metadata[[input$roi_column]]))
    })
    
    # compute once when button is pressed
    corr_results <- shiny::eventReactive(input$run_corr, {
      # Robust input validation
      if (is.null(input$roi_column) || is.null(input$roi_label) || input$roi_column == '' || input$roi_label == '') {
        return(NULL)
      }
      shiny::withProgress(message = 'Calculating feature vs distance correlations...', value = 0, {
        tryCatch({
          meta <- shared_data$updated.metadata
          expr <- full.intensity.matrix
          # Filter by selected samples
          selected_samples <- input$selected_samples
          expr <- expr[meta$Sample %in% selected_samples, , drop = FALSE]
          meta <- meta[meta$Sample %in% selected_samples, , drop = FALSE]
          shiny::incProgress(0.2)
          Sys.sleep(0.1)
          # distance from ROI
          roi_mask <- meta[[input$roi_column]] == input$roi_label
          dist_vec <- rep(NA, nrow(meta))
          # compute per sample if needed
          for (s in unique(meta$Sample)) {
            idx <- meta$Sample == s
            coords <- meta[idx, c("x", "y")]
            roi_coords <- coords[roi_mask[idx], ]
            if (nrow(roi_coords) > 0) {
              dists <- apply(coords, 1, function(pt)
                min(sqrt((pt[1] - roi_coords[,1])^2 + (pt[2] - roi_coords[,2])^2))
              )
              dist_vec[idx] <- dists
            }
          }
          shiny::incProgress(0.5)
          Sys.sleep(0.1)
          # correlate each feature with distance
          cors <- apply(expr, 2, function(f)
            suppressWarnings(stats::cor.test(f, dist_vec, method = input$corr_method))
          )
          shiny::incProgress(0.8)
          Sys.sleep(0.1)
          output_table = tibble::tibble(
            "feature" = colnames(expr),
            "correlation" = sapply(cors, function(x) x$estimate),
            "p.value" = sapply(cors, function(x) x$p.value)
          ) |>
            dplyr::mutate(adj.p = stats::p.adjust(.data$p.value, "BH"))
          # Merge annotation info (by m_z or feature name)
          anno_col = NULL
          if (!is.null(anno)) {
            # Try to merge by m_z if present, else by feature
            merge_col = intersect(c("m_z", "feature"), colnames(anno))
            if (length(merge_col) > 0) {
              output_table = merge(output_table, anno, by.x = "feature", by.y = merge_col[1], all.x = TRUE, sort = FALSE)
            }
          }
          shiny::incProgress(1)
          Sys.sleep(0.1)
          return(list('table'=output_table, 'dist_vec'=dist_vec, 'meta'=meta, 'expr'=expr))
        }, error = function(e) {
          return(NULL)
        })
      })
    })
    
    output$feature_select_ui <- shiny::renderUI({
      shiny::req(corr_results())
      ns <- session$ns
      shiny::selectizeInput(ns("selected_features"), "Select features to plot:",
                     choices = corr_results()$table[order(corr_results()$table$correlation),]$feature, selected = utils::head(corr_results()$table[order(corr_results()$table$correlation),]$feature), multiple = TRUE)
    })
    
    
    output$corr_table <- DT::renderDataTable({
      shiny::req(corr_results())
      corr_tab <- corr_results()$table[order(-corr_results()$table$correlation),]
      # If annotation columns exist, put them first
      anno_cols <- intersect(c("display_name", "name", "formula", "adduct", "Annotation"), colnames(corr_tab))
      if (length(anno_cols) > 0) {
        corr_tab <- corr_tab[, c(anno_cols, setdiff(names(corr_tab), anno_cols))]
      }
      DT::datatable(corr_tab, options = list(pageLength = 10), selection = 'none', rownames = FALSE)
    })
    
    volcanoPlot <- shiny::reactive({
      shiny::req(corr_results())
      shiny::withProgress(message = 'Rendering volcano plot...', value = 0, {
        shiny::incProgress(0.5)
        volcano_df <- corr_results()$table
        # Compose annotation string for tooltip if available
        anno_cols <- intersect(c("display_name", "name", "formula", "adduct"), colnames(volcano_df))
        anno_text <- if (length(anno_cols) > 0) {
          apply(volcano_df[, anno_cols, drop = FALSE], 1, function(row) {
            # Wrap display_name at 30 chars for readability
            txts <- lapply(seq_along(row), function(i) {
              val <- as.character(row[i])
              if (anno_cols[i] == "display_name" && nchar(val) > 30) {
                # Insert <br> every 30 chars
                paste(anno_cols[i], paste(strwrap(val, 30), collapse = "<br>"), sep = ": ")
              } else {
                paste(anno_cols[i], val, sep = ": ")
              }
            })
            paste(txts, collapse = "<br>")
          })
        } else {
          rep("", nrow(volcano_df))
        }
        volcano_df$log_adj_p = pmin(-log10(volcano_df$adj.p), 300)  # cap for better visualization
        p <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = .data$correlation, y = .data$log_adj_p,
                                   text = paste0(
                                     'Feature: ', .data$feature, '<br>',
                                     'Correlation: ', signif(.data$correlation, 3), '<br>',
                                     'p-value: ', signif(.data$p.value, 3), '<br>',
                                     'adj.p: ', signif(.data$adj.p, 3),
                                     ifelse(anno_text != "", paste0('<br>', anno_text), "")
                                   ))) +
          ggplot2::geom_point() +
          ggplot2::theme_minimal() +
          ggplot2::ylab("-log10(adj.p)")
        shiny::incProgress(1)
        return(p)
    })
    })

    output$volcano_plot <- plotly::renderPlotly({
      p_plotly = plotly::ggplotly(volcanoPlot(), tooltip = "text") |> plotly::layout(hoverlabel = list(bgcolor = "white", font = list(size = 12)), width = NULL)
      })

    output[['downloadVolcanoCorr']] <- utils_create_download_plot_handler(
      plot_func = volcanoPlot,
      filename_func = function() input[['volcanoCorrFileName']],
      width_func = function() input[['volcanoCorrWidth']],
      height_func = function() input[['volcanoCorrHeight']]
    )
    
    spatial_plot <- shiny::reactive({
      shiny::req(corr_results())
      shiny::withProgress(message = 'Rendering spatial distance plot...', value = 0, {
        shiny::incProgress(0.5)
        overview = corr_results()$meta
        overview$dist_vec = corr_results()$dist_vec
        # Only keep samples with at least one non-NA distance value
        valid_samples <- unique(overview$Sample[!is.na(overview$dist_vec)])
        overview <- overview[overview$Sample %in% valid_samples, , drop = FALSE]
        p <- ggplot2::ggplot(overview, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data$dist_vec, fill = .data$dist_vec)) +
          ggplot2::geom_tile() +
          ggplot2::facet_wrap(~.data$Sample, scales = 'free', nrow = max(1,floor(sqrt(length(unique(overview$Sample)))/1.5))) +
          ggplot2::theme_classic() +
          ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                        axis.text.x=ggplot2::element_blank(),
                        axis.ticks.x=ggplot2::element_blank(),
                        axis.line.x = ggplot2::element_blank(),
                        axis.title.y=ggplot2::element_blank(),
                        axis.text.y=ggplot2::element_blank(),
                        axis.ticks.y=ggplot2::element_blank(),
                        axis.line.y = ggplot2::element_blank()) +
          ggplot2::theme(aspect.ratio = 1)
        shiny::incProgress(1)
      })
      return(p)
    })
    output$spatial_plot <- shiny::renderPlot({
      spatial_plot()
    })
    output[['downloadSpatialCorr']] <- utils_create_download_plot_handler(
      plot_func = spatial_plot,
      filename_func = function() input[['spatialCorrFileName']],
      width_func = function() input[['spatialCorrWidth']],
      height_func = function() input[['spatialCorrHeight']]
    )
    
    # Smoothed intensity plot

    smooth_line <- shiny::reactive({
      shiny::req(input$selected_features, corr_results())
      shiny::withProgress(message = 'Rendering smoothed intensity plot...', value = 0, {
        meta <- corr_results()$meta
        expr <- corr_results()$expr
        dist_vec <- corr_results()$dist_vec
        shiny::incProgress(0.2)
        plot_df <- tibble::tibble(
          Distance = dist_vec,
          Sample = meta$Sample
        )
        # gather intensity for selected features
        for (feat in input$selected_features) {
          plot_df[[feat]] <- expr[, feat]
        }
        plot_df[,input$selected_features]=scale(plot_df[,input$selected_features])
        if (input$split_smoothed){
          plot_df$split = meta[,input$split_var]
        }
        shiny::incProgress(0.5)
        plot_df <- tidyr::pivot_longer(plot_df, 
                                       cols = dplyr::all_of(input$selected_features),
                                       names_to = "Feature", values_to = "Intensity")
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Distance, y = .data$Intensity, color = .data$Feature)) +
          ggplot2::geom_smooth(method = 'gam', se = FALSE, linewidth = 1.2) +
          ggplot2::theme_classic() +
          ggplot2::labs(title = "Smoothed Intensity vs Distance from ROI",
               x = "Distance", y = "Intensity") +
          ggplot2::theme(legend.position = "bottom")
        if (input$split_smoothed){
          p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$split))
        }
        shiny::incProgress(1)
      })
    return(p)
    })
    
    output$smooth_plot <- shiny::renderPlot({
      smooth_line()
    })   
    output[['downloadSmooth']] <- utils_create_download_plot_handler(
      plot_func = smooth_line,
      filename_func = function() input[['smoothFileName']],
      width_func = function() input[['smoothWidth']],
      height_func = function() input[['smoothHeight']]
    )
    
    output$download_corr <- shiny::downloadHandler(
      filename = function() paste0("feature_distance_correlations.csv"),
      content = function(file) {
        utils::write.csv(corr_results()$table, file, row.names = FALSE)
      }
    )
  })
}
