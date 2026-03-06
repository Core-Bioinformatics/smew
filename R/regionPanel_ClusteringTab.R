
#' Performs clustering and batch correction
#'
#' @description UI and server logic for spatial clustering analysis in the SMEW app. Supports multiple clustering algorithms (k-means, Louvain), preprocessing (PCA, optional Harmony batch correction), and interactive visualization.
#'
#' @details
#' \itemize{
#'  \item{Sample selection and dimensionality reduction}
#'  \item{Multiple clustering algorithms (k-means, Louvain)}
#'  \item{Batch correction with Harmony}
#'  \item{Interactive cluster visualization and persistent cluster assignment}
#'  \item{Download handlers for cluster plots}
#' }
#' @param id Shiny module id
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param full.metadata Data frame of full sample metadata
#' @param full.intensity.matrix Numeric matrix of expression data
#' @param show Logical; whether to render the panel (default: TRUE, UI only)
#' @param anno Data frame of peak annotations (server only)
#' @param shared_data Reactive or shared data object (server only)
#'
#' @return UI: A shiny::tabPanel object for the clustering panel. Server: None; called for side effects in Shiny module.
#'
#' @keywords internal
#' @name RegionPanel_ClusterTab
#' @rdname RegionPanel_ClusterTab
#' @export
RegionPanel_ClusterTabUI <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, show = TRUE) {
  ns <- shiny::NS(id)
  # add option to name cluster and retain it
  if(show){
    shiny::tabPanel(
      'Clustering',
      shiny::br(),
      bslib::accordion(
        bslib::accordion_panel(
          title = "Information",
          icon = bsicons::bs_icon("info-circle"),
          shiny::tags$ul(
            shiny::tags$li("Select samples and perform dimensionality reduction (PCA with or without Harmony batch correction) then the selected clustering appraoch to identify metabolic regions."),
            shiny::tags$li("Visualize cluster assignments interactively on spatial plots.")
          )
        ),
        bslib::accordion_panel(
          title = "Sample selection",
          icon = bsicons::bs_icon("gear"),
          shiny::selectInput(
            inputId = ns("samplesToFactor"),
            label = "Select samples for infer dimensionality reduction",
            choices = base::unique(bulk.metadata[,1]),
            selected = base::unique(bulk.metadata[,1]),
            multiple = TRUE, width = '100%'
          ),
          shiny::radioButtons(ns("logTransform"), "Log-transform data?", choices = c('Log intensity' = TRUE, 'Raw intensity' = FALSE), selected = TRUE),
          shiny::numericInput(inputId = ns("seed"),
            label = 'Set seed',
            value = 23,
            min = 1,
            max = 1000,
            step = 1)
        ),
        bslib::accordion_panel(
          title = "Downloads",
          icon = bsicons::bs_icon("download"),
          shiny::tags$div(
            shiny::tags$h4("Download output plots"),
            shiny::fluidRow(
              shiny::column(6,
                shiny::tags$strong("Cluster plot"),
                shiny::textInput(ns('clusterPlotFileName'), 'File name', value ='Clusters.png'),
                shiny::numericInput(ns('clusterPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('clusterPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadClusterPlot'), 'Download cluster plot')
              )
            )
          )
        ),
        id = ns("acc"),
        open = "Sample selection"
      ),
      shiny::fluidRow(
        shiny::column(12,
          bslib::card(
            class = "mb-3",
            style = "width:100%;",
            shiny::tags$h4("Clustering Settings"),
            shiny::fluidRow(
              shiny::column(6,
                shiny::sliderInput(ns('num_components'),label = 'Number of PCs',min = 5,max=50,value = 30),
                shiny::checkboxInput(ns('run_harmony'),label = 'Run harmony correction',value=F),
                shiny::conditionalPanel(condition = "input.run_harmony",
                  shiny::selectInput(ns('harmony_meta'),label = 'Correcting metadata',choices = colnames(bulk.metadata)),
                  shiny::numericInput(ns('harmony_theta'),label = 'Harmony theta value',value = 0.5,min = 0,max = 10,step = 0.1),
                  ns = ns),
                shiny::actionButton(
                  inputId = ns("run_preprocessing"),
                  label = "Run pre-processing",
                  icon = shiny::icon("play"),
                  style = "margin-top:10px;width:100%;"
                )
              ),
              shiny::column(6,
                shiny::selectInput(ns("clusteringMethod"),label = 'Clustering method', choices = c('k-means','Louvain')),
                shiny::conditionalPanel(condition = "input.clusteringMethod == 'k-means'",
                  shiny::numericInput(ns('numClusters'),label = 'Numbers of clusters',min = 1,max=30,step = 1,value=5),
                  ns=ns),
                shiny::conditionalPanel(condition = "input.clusteringMethod == 'Louvain'",
                  shiny::numericInput(inputId = ns("resolution"),'Resolution for clustering',value = 0.1,min = 0.01,max = 3,step = 0.01),
                  ns=ns),
                shiny::actionButton(
                  inputId = ns("run_clustering"),
                  label = "Run clustering",
                  icon = shiny::icon("play"),
                  style = "margin-top:10px;width:100%;"
                )
              )
            )
          )
        )
      ),
      shiny::plotOutput(ns('clusterPlot'),height = '600px'),
      shiny::textInput(ns('clusterName'),label = 'Name for clusters'),
      shiny::actionButton(ns('lock_cluster'),'Add clusters to object')
    )
  }
}

#' Server logic for Region-level Clustering Panel
#'
#' @rdname RegionPanel_ClusterTab
#' @param id Shiny module id
#' @param full.intensity.matrix Matrix of intensities (features x samples)
#' @param full.metadata Data frame of full sample metadata
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param anno Data frame of peak annotations
#' @param shared_data Reactive or shared data object
#' @return None; called for side effects in Shiny module
#' @export
RegionPanel_ClusterTabServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, anno, shared_data) {
  
  shiny::moduleServer(id, function(input, output, session) {
    
    get_subset_exp <- shiny::reactive({
      idx <- full.metadata$Sample %in% input[["samplesToFactor"]]
      current.intensity.matrix <- full.intensity.matrix[idx, ]
      if (input[['logTransform']]) {
        current.intensity.matrix <- log2(current.intensity.matrix + 1)
      }
      list(
        exp = current.intensity.matrix,
        meta = full.metadata[idx, ]
      )
    })
    
    # --- Preprocessing step (only runs on button press) ---
      preproc_obj <- shiny::reactive({
        shiny::withProgress(message = 'Running pre-processing...', value = 0, {
          dat <- get_subset_exp()
          current.intensity.matrix <- dat$exp
          current.metadata <- dat$meta
          shiny::incProgress(0.2)
           # Remove features (columns) with zero variance
           var_vec <- matrixStats::colVars(current.intensity.matrix)
           keep <- var_vec > 0
           filtered_matrix <- current.intensity.matrix[, keep, drop = FALSE]
           # Create SingleCellExperiment object
           spe <- SingleCellExperiment::SingleCellExperiment(
             assays = list(counts = t(filtered_matrix)),
             colData = current.metadata
           )
           shiny::incProgress(0.5)
           # PCA
           spe <- scater::runPCA(spe, ncomponents = input$num_components,
                               exprs_values = "counts", name = "PCA",
                               scale = TRUE, ntop = ncol(filtered_matrix))
          # Harmony batch correction (optional)
          if (input$run_harmony) {
            batch <- current.metadata[, input$harmony_meta]
            pca_mat <- SingleCellExperiment::reducedDims(spe)[['PCA']]
            harmony_mat <- harmony::RunHarmony(
              pca_mat,
              meta_data = data.frame(Batch = batch),
              vars_use = 'Batch',
              theta = 0.5
            )
            SingleCellExperiment::reducedDims(spe)[['Harmony']] <- harmony_mat
          }
          shiny::incProgress(0.9)
          list(metadata = current.metadata, spe_obj = spe)
        })
      }) |>
        shiny::bindCache(
          utils::head(full.intensity.matrix),
          input$samplesToFactor,
          input$num_components,
          input$run_harmony,
          input$harmony_meta
        ) |> shiny::bindEvent(input$run_preprocessing)
    
    # cache to safely store latest clustering result (nonreactive)
    cluster_cache <- shiny::reactiveVal(NULL)
    
    # --- Clustering step (same as before) ---
      clustering_obj <- shiny::reactive({
        shiny::withProgress(message = 'Running clustering...', value = 0, {
          obj <- preproc_obj()
          shiny::req(obj)
          current.metadata <- obj$metadata
          spe <- obj$spe_obj
          shiny::incProgress(0.3)
          # Choose reduction
          reduction <- if (input$run_harmony && !is.null(SingleCellExperiment::reducedDims(spe)[['Harmony']])) {
            SingleCellExperiment::reducedDims(spe)[['Harmony']]
          } else {
            SingleCellExperiment::reducedDims(spe)[['PCA']]
          }
          if (input$clusteringMethod == "Louvain") {
            g <- scran::buildSNNGraph(t(reduction), k = 20)
            cl <- igraph::cluster_louvain(g)
            current.metadata$cluster <- factor(cl$membership)
          } else {
            base::set.seed(input$seed)
            clusters <- stats::kmeans(reduction, centers = input$numClusters)
            current.metadata$cluster <- base::factor(clusters$cluster)
          }
          shiny::incProgress(0.9)
          cluster_cache(current.metadata)
          current.metadata
        })
      }) |> shiny::bindEvent(input$run_preprocessing, input$run_clustering)
    
    # --- Plot output ---
    output$clusterPlot <- shiny::renderPlot({
      df <- clustering_obj()
      ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, fill = .data$cluster, color = .data$cluster)) +
        ggplot2::geom_tile() +
        ggplot2::theme_classic() +
        ggplot2::facet_wrap(~.data$Sample, nrow = max(1, floor(sqrt(length(unique(df$Sample)) / 2)))) +
        ggplot2::theme_void() +
        ggplot2::theme(aspect.ratio = 1)
    })
    output[['downloadClusterPlot']] <- utils_create_download_plot_handler(
      plot_func = function(){
        df <- clustering_obj()
        ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, fill = .data$cluster, color = .data$cluster)) +
          ggplot2::geom_tile() +
          ggplot2::theme_classic() +
          ggplot2::facet_wrap(~.data$Sample, nrow = max(1, floor(sqrt(length(unique(df$Sample)) / 2)))) +
          ggplot2::theme_void() +
          ggplot2::theme(aspect.ratio = 1)
      },
      filename_func = function() input[['clusterPlotFileName']],
      width_func = function() input[['clusterPlotWidth']],
      height_func = function() input[['clusterPlotHeight']]
    )
    
    # --- Lock cluster (read from cache, NOT reactive) ---
    shiny::observeEvent(input$lock_cluster, {
      meta <- shared_data$updated.metadata
      clusterdf <- shiny::isolate(cluster_cache())  # << isolate to break dependency
      shiny::req(clusterdf)
      
      if (!input$clusterName %in% base::colnames(meta)) {
        meta[[input$clusterName]] <- NA
      }
      
      idx <- base::match(meta$pixel_id, clusterdf$pixel_id)
      to_update <- !base::is.na(idx)
      meta[[input$clusterName]][to_update] <- clusterdf$cluster[idx[to_update]]
      
      shared_data$updated.metadata <- meta
    })
    
  })
}

