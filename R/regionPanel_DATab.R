#' Performs and visualises region-based differential analysis
#'
#' @description UI and server logic for the Region-Based Differential Analysis (DA) panel, enabling statistical comparison of user-defined regions (e.g., clusters, spatial regions, or metadata-defined groups) in spatial omics data. Supports pseudobulking, group selection, statistical testing, thresholding, and visualisation of results (tables, boxplots, volcano plots). Users can download results and plots for downstream analysis.
#'
#' @details
#' \itemize{
#'   \item{Allows users to perform differential analysis between two groups of pixels or regions, such as clusters or metadata-defined groups.}
#'   \item{Select a grouping variable (e.g., experimental condition, cluster), and specify the two groups to compare.}
#'   \item{Choose a statistical test (t-test or Wilcoxon rank sum), adjust for multiple testing, and set thresholds for significance.}
#'   \item{Results are displayed as interactive volcano plots, tables, and boxplots.}
#'   \item{Download significant features and plots for downstream analysis or reporting.}
#'   \item{Supports pseudobulking of regions for robust group-level comparison.}
#' }
#'
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame of bulk sample metadata (UI)
#' @param full.metadata Data frame of full sample metadata (UI, server)
#' @param show Logical; whether to show the panel (default TRUE, UI)
#' @param full.intensity.matrix Matrix of intensities (features x samples, server)
#' @param shared_data Reactive or shared data object (server)
#' @param anno Data frame of peak annotations (server)
#' @return UI: A shiny tabPanel object for the region-based DA tab. Server: None; called for side effects in Shiny module.
#' @export
#' @name RegionPanel_DATab
#' @rdname RegionPanel_DATab
RegionPanel_DATabUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- shiny::NS(id)
  
  if(show){
    shiny::tabPanel(
      'Region Differential Analysis',
       bslib::accordion(
          bslib::accordion_panel(
            title = "Information",
            icon = bsicons::bs_icon("arrow-right-circle"),
            shiny::tags$ul(
              shiny::tags$li("Perform differential analysis between two groups of pixels (e.g., clusters, regions, or metadata-defined groups)."),
              shiny::tags$li("Select a grouping variable, such as experimental conditions, and specify the two groups to compare."),
              shiny::tags$li("Choose a statistical test, adjust for multiple testing, and view results as volcano plots and tables."),
              shiny::tags$li("Download significant features and plots for downstream analysis or reporting.")
            )
          ),
          bslib::accordion_panel(
            title = "Downloads",
            icon = bsicons::bs_icon("download"),
            shiny::tags$div(
              shiny::tags$h4("Download results table"),
              shiny::fluidRow(
                shiny::column(6,
                  shiny::textInput(ns('fileName'),'File name for download', value ='DIAset.csv', placeholder = 'DIAset.csv'),
                  shiny::downloadButton(ns('download'), 'Download Table')
                ),
                shiny::column(6,
                  shiny::tags$strong("Boxplot for selected peak"),
                  shiny::textInput(ns('boxplotFileName'), 'File name', value ='PeakBoxplot.png'),
                  shiny::numericInput(ns('boxplotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                  shiny::numericInput(ns('boxplotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                  shiny::downloadButton(ns('downloadBoxplot'), 'Download Boxplot')
                )
              )
            )
          ),
          id = ns("acc"),
          open = NULL
        ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput(
            inputId = ns("regionToGroupOn"),
            label = "Select region to group on",
            choices = c(colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','x_tf','y_tf','Sample',colnames(bulk.metadata)))]),
            selected = c(colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','x_tf','y_tf','Sample',colnames(bulk.metadata)))])[1]
          ),
          shiny::numericInput(inputId = ns('minimumPixelsPerSample'),
            label = 'Minimum number of pixels per region per sample to be considered',
            value = 10,
            min = 1,
            max = 1000),
          shiny::selectInput(ns('condition'), 'Metadata column to use:', c('AllSamples',colnames(bulk.metadata)[-1]),
            selected = 'AllSamples'),
          shiny::actionButton(ns('pseudoBulk'), label = 'Start Pseudobulking'),
          shiny::selectInput(ns('variable1'), 'Condition 1:', c(),multiple = T),
          shiny::selectInput(ns('variable2'), 'Condition 2:', c(),multiple = T),
          shiny::selectInput(ns('pipeline'), 'DE pipeline:', c("t-test", "Wilcox rank sum")),
          shiny::sliderInput(ns('pvalThreshold'), label = 'Adjusted p-value threshold',
            min = 0, value = 0.05, max = 1, step = 0.005),
          shiny::sliderInput(ns('lfcThreshold'), label = 'log2 fold change threshold',
            min = 0, value = 1, max = 5, step = 0.1),
          shiny::actionButton(ns('goDE'), label = 'Start differential intensity analysis')
        ),
        shiny::mainPanel(
          shinyjs::useShinyjs(),
          shiny::textOutput(ns('print_pseudobulk')),
          DT::DTOutput(ns('data')),
          shiny::div("Click a peak row to show its boxplot across compared conditions."),
          shiny::plotOutput(ns('peakBoxplot'))
        )
      )
    )
  }else{
    NULL
  }
}


#' @rdname RegionPanel_DATab
#' @export
RegionPanel_DATabServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, shared_data, anno){

  # Validate non-reactive constants
  stopifnot({
    !shiny::is.reactive(anno)
  })
  
  shiny::moduleServer(id, function(input, output, session){
    shiny::observe({
      shiny::req(!is.null(shared_data$updated.metadata))
      choices <- colnames(shared_data$updated.metadata)[!(colnames(shared_data$updated.metadata) %in% c('pixel_id','x','y','x_tf','y_tf','Sample', colnames(bulk.metadata)))]
      shiny::validate(shiny::need(length(choices) > 0, 'No region columns available to group on'))
      shiny::updateSelectInput(session, 'regionToGroupOn', choices = choices, selected = choices[1])
    })
    
    pseudobulk_samples <- shiny::reactive({
      shiny::req(input[['regionToGroupOn']], input[['minimumPixelsPerSample']])
      shiny::validate(shiny::need(input[['minimumPixelsPerSample']] >= 1, 'Minimum pixels must be >= 1'))
      pseudobulked <- NULL
      shiny::withProgress(message = 'Pseudobulking regions...', value = 0, {
        shiny::incProgress(0.2, detail = 'Preparing inputs')
        pseudobulked <- region_utils_create_bulk_exp_regions(as.data.frame(full.intensity.matrix),
                                                shared_data$updated.metadata,
                                                bulk.metadata,
                                                colnames(bulk.metadata)[1],
                                                shared_data$updated.metadata[,input[['regionToGroupOn']]],
                                                minimum.pixels = input[['minimumPixelsPerSample']])
        shiny::incProgress(0.6, detail = 'Aggregating per-sample means')
        shiny::incProgress(0.2, detail = 'Finalizing pseudobulk outputs')
      })
      cond_values <- unique(pseudobulked$metadata[[input[["condition"]]]])
      shiny::updateSelectInput(session, 'variable1', choices = cond_values, selected = cond_values[1])
      if (length(cond_values) >= 2){
        shiny::updateSelectInput(session, 'variable2', choices = cond_values, selected = cond_values[2])
      } else {
        shiny::updateSelectInput(session, 'variable2', choices = cond_values)
      }
      
      return(pseudobulked)
    }) |> shiny::bindEvent(input[["pseudoBulk"]])
    
    output[['print_pseudobulk']]<-shiny::renderText({
      shiny::req(pseudobulk_samples())
      paste0('Conditions extracted: ', paste(unique(pseudobulk_samples()$metadata[[input[["condition"]]]]), collapse = ', '))
      }) |> shiny::bindEvent(input[["pseudoBulk"]])
    
    shiny::observe({
      condition.indices <- shared_data$updated.metadata[[input[["condition"]]]] %in% c(input[['variable1']], input[['variable2']])
      choices <- c("t-test", "Wilcox rank sum")
      shiny::updateSelectInput(session, 'pipeline', choices = choices)
    })
    
    
    DEresults <- shiny::reactive({
      shinyjs::disable("goDE")
      shiny::req(pseudobulk_samples())
      pseudobulk.metadata = pseudobulk_samples()$metadata
      pseudobulk.intensity.matrix = pseudobulk_samples()$intensity.matrix
      pseudobulk.intensity.matrix = pseudobulk.intensity.matrix[,pseudobulk.metadata[,1]]
      shiny::req(input[['variable1']], input[['variable2']], input[["condition"]])
      condition.indices <- pseudobulk.metadata[[input[["condition"]]]] %in% c(input[['variable1']], input[['variable2']])
      grp <- pseudobulk.metadata[[input[["condition"]]]][condition.indices]
      shiny::validate(shiny::need(sum(grp == input[['variable1']]) >= 2, 'Condition 1 must have at least 2 samples'))
      shiny::validate(shiny::need(sum(grp == input[['variable2']]) >= 2, 'Condition 2 must have at least 2 samples'))
      DEtable <- NULL
      shiny::withProgress(message = 'Running differential analysis...', value = 0, {
        shiny::incProgress(0.2, detail = 'Preparing inputs')
        DEtable <- bulk_utils_DA(
          intensity_matrix = pseudobulk.intensity.matrix[, condition.indices],
          condition = grp,
          var1 = input[['variable1']],
          var2 = input[['variable2']],
          test = input[["pipeline"]],
          anno = anno
        )
        shiny::incProgress(0.6, detail = 'Computing statistics')
        shiny::incProgress(0.2, detail = 'Filtering and ordering results')
      })
      DEtableSubset <- dplyr::arrange(
        dplyr::filter(DEtable, .data$pvalAdj < input[["pvalThreshold"]] & abs(.data$lfc) > input[['lfcThreshold']]),
        .data$pvalAdj
      )
      
      #the thresholds are returned here so that MA/volcano and table display
      #don't use new thresholds without the button being used
      shinyjs::enable("goDE")
      return(list('DEtable' = DEtable,
          "DEtableSubset" = DEtableSubset,
          'pvalThreshold' = input[["pvalThreshold"]],
          'lfcThreshold' = input[['lfcThreshold']]))
    }) |>
      shiny::bindCache(shared_data, input[["condition"]],
                input[['variable1']], input[['variable2']], input[["pipeline"]],
                input[["pvalThreshold"]],input[['lfcThreshold']]) |>
      shiny::bindEvent(input[["goDE"]])
    
    #Define output table (only DE peaks)
    dataTable <- shiny::reactive({
      DT::formatSignif(DT::datatable(DEresults()$DEtableSubset, selection = 'single'), columns = c('pval', 'pvalAdj','lfc','log2_intensity'), digits = 3)
    })
    
    output[['data']] <- DT::renderDataTable(dataTable())
    
    #DE data download
    output[['download']] <- shiny::downloadHandler(
      filename = function() {
        paste(input[['fileName']])
      },
      content = function(file) {
        utils::write.csv(x = DEresults()$DEtableSubset, file = file, row.names = FALSE)
      }
    )
    
    # Single peak selection: show boxplot for clicked peak
    peakBox <- shiny::reactive({
      # Get selected row index from DT
      selected_row <- input$data_rows_selected
      if (is.null(selected_row) || length(selected_row) != 1) return(NULL)
      # Get peak name
      peak_name <- DEresults()$DEtableSubset$m_z[selected_row]
      # Get pseudobulk data
      pseudobulk <- pseudobulk_samples()
      pseudobulk_int <- pseudobulk$intensity.matrix
      pseudobulk_meta <- pseudobulk$metadata

      # Defensive: check if peak exists in matrix
      if (!(peak_name %in% rownames(pseudobulk_int))) return(NULL)
      intensity_vec <- pseudobulk_int[peak_name, ]

      # If group_vec contains combined group_cluster (e.g. group1_cluster1), split it
      split_vals <- do.call(rbind, strsplit(as.character(pseudobulk_meta[[ input$condition ]]), '_', fixed = TRUE))
      group_vec <- stringr::word(as.character(pseudobulk_meta[[ input$condition ]]), sep = '_', 1, -2)
      cluster_vec <- stringr::word(as.character(pseudobulk_meta[[ input$condition ]]), sep = '_', -1)
      df <- data.frame(
        intensity = intensity_vec,
        group = group_vec,
        cluster = cluster_vec
      )
      if (nrow(df) == 0 || all(is.na(df$intensity))) {
        graphics::plot.new(); graphics::text(0.5, 0.5, "No data available for this peak.", cex = 1.2)
        return()
      }
      # Always color by cluster, and if only one cluster, use a default color
      # Always color by cluster, even for AllSamples
      n_clusters <- length(unique(df$cluster[!is.na(df$cluster)]))
      palette <- if (n_clusters == 1) {
        c("#3182bd")
      } else if (n_clusters <= 8) {
        RColorBrewer::brewer.pal(n_clusters, "Set2")
      } else {
        grDevices::rainbow(n_clusters)
      }
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$intensity, fill = as.factor(.data$cluster))) +
        ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8)) +
        ggplot2::scale_fill_manual(values = palette) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste('Peak', peak_name, 'across', ifelse(input$condition == 'AllSamples', 'all samples', 'groups and clusters')),
          x = ifelse(input$condition == 'AllSamples', 'Sample', input$condition),
          y = 'Intensity',
          fill = input$regionToGroupOn
        )
      p
    })
    
    output[['peakBoxplot']] <- shiny::renderPlot({
      peakBox()
    })

    output[['downloadBoxplot']] <- utils_create_download_plot_handler(
      plot_func = peakBox,
      filename_func = function() input[['boxplotFileName']],
      width_func = function() input[['boxplotWidth']],
      height_func = function() input[['boxplotHeight']]
    )
    
  })
}
