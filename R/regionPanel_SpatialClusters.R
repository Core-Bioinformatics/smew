get_combined_cluster_colors <- function(comb) {
  combined_cluster_levels <- sort(unique(as.integer(comb)))
  stats::setNames(
    RColorBrewer::brewer.pal(max(3, length(combined_cluster_levels)), "Set2")[combined_cluster_levels],
    as.character(combined_cluster_levels)
  )
}

##' RegionSpatialClusterPanelUI
##'
##' UI for the region-level spatial clustering tab in the SMEW app.
##'
##' @param id Shiny module id
##' @param bulk.metadata Data frame of bulk sample metadata
##' @param full.metadata Data frame of full sample metadata
##' @param show Logical; whether to show the panel (default TRUE)
##' @return A shiny tabPanel object for the spatial clustering tab
##' @export
RegionSpatialClusterPanelUI <- function(id, bulk.metadata, full.metadata, show = TRUE) {
  ns <- shiny::NS(id)
  if (show) {
    shiny::tabPanel(
      'Spatial Clustering',
      shiny::br(),
      bslib::accordion(
        bslib::accordion_panel(
          title = "Information",
          icon = bsicons::bs_icon("info-circle"),
          shiny::tags$ul(
            shiny::tags$li("Run BayesSpace clustering on one or more selected samples. Choose the number of clusters and MCMC iterations."),
            shiny::tags$li("Visualize spatial cluster assignments for each sample, and add cluster labels to shared data for downstream analysis."),
            shiny::tags$li("Pseudobulk clusters across samples and perform PCA to explore similarities between clusters."),
            shiny::tags$li("View a similarity heatmap of pseudobulked clusters, with annotation bars for combined cluster assignments."),
            shiny::tags$li("Perform combined clustering across all samples to identify shared cluster patterns."),
            shiny::tags$li("Visualize combined clusters spatially across all samples in a unified plot."),
            shiny::tags$li("Download all output plots (spatial clusters, PCA, heatmap, combined spatial) in publication-ready format."),
            shiny::tags$li("Note: Clustering and pseudobulking may take several minutes depending on data size and number of samples.")
          )
        ),
        bslib::accordion_panel(
          title = "Sample selection",
          icon = bsicons::bs_icon("gear"),
          shiny::selectInput(
            inputId = ns("sampleToCluster"),
            label = "Select sample(s) for clustering",
            choices = unique(bulk.metadata[, 1]),
            selected = unique(bulk.metadata[1, 1]),
            multiple = TRUE,
            width = '100%'
          ),
          shiny::helpText("You can select one or more samples. Clustering will be run separately for each sample.")
        ),
        bslib::accordion_panel(
          title = "Downloads",
          icon = bsicons::bs_icon("download"),
          shiny::tags$div(
            shiny::tags$h4("Download output plots"),
            shiny::fluidRow(
              shiny::column(6,
                shiny::tags$strong("Spatial clusters"),
                shiny::textInput(ns('spatialPlotFileName'), 'File name', value = 'spatialClusters.png'),
                shiny::numericInput(ns('spatialPlotWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                shiny::numericInput(ns('spatialPlotHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                shiny::downloadButton(ns('downloadSpatial'), 'Download spatial clusters')
              ),
              shiny::column(6,
                shiny::tags$strong("PCA plot"),
                shiny::textInput(ns('pcaPlotFileName'), 'File name', value = 'pseudobulkPCA.png'),
                shiny::numericInput(ns('pcaPlotWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                shiny::numericInput(ns('pcaPlotHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                shiny::downloadButton(ns('downloadPCA'), 'Download PCA plot')
              )
            ),
            shiny::fluidRow(
              shiny::column(6,
                shiny::tags$strong("Heatmap"),
                shiny::textInput(ns('heatmapPlotFileName'), 'File name', value = 'pseudobulkHeatmap.png'),
                shiny::numericInput(ns('heatmapPlotWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                shiny::numericInput(ns('heatmapPlotHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                shiny::downloadButton(ns('downloadHeatmap'), 'Download heatmap')
              ),
              shiny::column(6,
                shiny::tags$strong("Combined spatial plot"),
                shiny::textInput(ns('combinedPlotFileName'), 'File name', value = 'combinedSpatial.png'),
                shiny::numericInput(ns('combinedPlotWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                shiny::numericInput(ns('combinedPlotHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                shiny::downloadButton(ns('downloadCombinedSpatial'), 'Download combined spatial plot')
              )
            )
          )
        ),
        id = ns("acc"),
        open = "Sample selection"
      ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput(ns('numClusters'), 'Number of spatial clusters to infer', multiple = FALSE, choices = 2:10),
          shiny::numericInput(ns('nRep'), 'Number of MCMC iterations for BayesSpace', min = 10, max = 10000, value = 100, step = 10),
          shiny::numericInput(ns('burnIn'), 'Burn-in period (iterations to exclude)', min = 1, max = 1000, value = 10, step = 10),
          shiny::actionButton(
            inputId = ns("run_clustering"),
            label = "Run BayesSpace clustering",
            icon = shiny::icon("play")
          ),
          shiny::hr(),
          shiny::textInput(ns('clusterName'), 'Name for sample cluster column', value = ''),
          shiny::actionButton(
            inputId = ns("add_to_shared"),
            label = "Add sample clusters to object",
            icon = shiny::icon("plus-circle")
          ),
          shiny::br(),
        ),
        shiny::mainPanel(
          shiny::plotOutput(ns('plotClusters'), height = '600px'),
        )),
        shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::h4('Combined clustering and pseudobulking'),
          shiny::conditionalPanel(
            condition = "input.sampleToCluster.length > 1",
            ns = ns,
            shiny::numericInput(ns('combined_n_clusters'), 'Number of combined clusters (across samples)', min = 2, max = 20, value = 4, step = 1),
            shiny::hr(),
            shiny::textInput(ns('combinedClusterName'), 'Name for combined cluster column', value = ''),
            shiny::actionButton(
              inputId = ns("add_combined_to_shared"),
              label = "Add combined clusters to object",
              icon = shiny::icon("plus-square")
            )
          )
        ),
          shiny::mainPanel(
            shiny::fluidRow(
              shiny::column(6,
                shiny::uiOutput(ns('pca_or_disabled'))
              ),
              shiny::column(6,
                shiny::uiOutput(ns('heatmap_or_disabled'))
              )
            ),
            shiny::uiOutput(ns('combined_spatial_or_disabled'))
          )
        )
    )
  } else {
    NULL
  }
}

##' RegionSpatialClusterPanelServer
##'
##' Server logic for the region-level spatial clustering tab in the SMEW app.
##'
##' @param id Shiny module id
##' @param bulk.metadata Data frame of bulk sample metadata
##' @param full.metadata Data frame of full sample metadata
##' @param full.intensity.matrix Matrix of intensities (features x samples)
##' @param anno Data frame of peak annotations
##' @param shared_data Reactive or shared data object
##' @return None; called for side effects in Shiny module
##' @export
RegionSpatialClusterPanelServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, anno = NULL, shared_data = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    # Store clustering result in a reactiveVal
    cluster_result <- shiny::reactiveVal(NULL)
    pseudobulk_result <- shiny::reactiveVal(NULL)

    shiny::observeEvent(input$run_clustering, {
      sample_col <- colnames(bulk.metadata)[1]
      sample_ids <- input$sampleToCluster
      all_results <- list()
      all_pseudobulk <- list()
      n_samples <- length(sample_ids)
      shiny::withProgress(message = 'Running BayesSpace clustering...', value = 0, {
        for (i in seq_along(sample_ids)) {
          sample_id <- sample_ids[i]
          intensity.sub <- full.intensity.matrix[full.metadata[, sample_col] == sample_id, ]
          metadata.sub <- full.metadata[full.metadata[, sample_col] == sample_id, ]
          metadata.sub$Sample <- metadata.sub[, sample_col]
          pos <- metadata.sub[, c('x', 'y')]
          rownames(pos) <- metadata.sub$pixel_id
          rownames(intensity.sub) <- metadata.sub$pixel_id
          metadata.sub$pixel_id <- paste0('pixel_', metadata.sub$pixel_id)
          rownames(intensity.sub) <- paste0('pixel_', rownames(intensity.sub))
          shiny::incProgress(0.1 / n_samples)

          # BayesSpace clustering only
          colData <- metadata.sub[, c('x_tf','y_tf')]
          colnames(colData)[1:2] <- c('col', 'row')
          shiny::incProgress(0.2 / n_samples)

          anno_df <- data.frame(anno)
          rownames(anno_df) <- paste0('mz_', anno_df$m_z)
          colnames(intensity.sub) <- paste0('mz_', colnames(intensity.sub))
          rowData <- anno_df[colnames(intensity.sub), ]
          shiny::incProgress(0.2 / n_samples)

          sce <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = methods::as(t(intensity.sub), "dgCMatrix")),
            rowData = rowData,
            colData = colData
          )
          S4Vectors::metadata(sce)$BayesSpace.data <- list(platform = 'ST', is.enhanced = FALSE)
          SingleCellExperiment::logcounts(sce) <- log2(SingleCellExperiment::counts(sce) + 1)
          shiny::incProgress(0.2 / n_samples)
          sce <- scater::runPCA(sce, subset_row = rownames(SingleCellExperiment::counts(sce)), ncomponents = 30,
            exprs_values = 'logcounts', BSPARAM = BiocSingular::ExactParam())
          set.seed(23)
          shiny::incProgress(0.1 / n_samples)
          sce <- BayesSpace::spatialCluster(
            sce,
            q = as.numeric(input[['numClusters']]),
            platform = "ST",
            d = 30,
            nrep = input[['nRep']],
            burn.in = input[['burnIn']],
            init.method = "mclust",
            model = "t",
            gamma = 2
          )
          shiny::incProgress(0.2 / n_samples)
          metadata.sub$cluster <- factor(sce$spatial.cluster, levels = 1:as.numeric(input[['numClusters']]))
          all_results[[sample_id]] <- metadata.sub

          # Pseudobulk: mean intensity for each cluster
          cluster_ids <- unique(metadata.sub$cluster)

          pb <- sapply(cluster_ids, function(cl) {
            mask <- metadata.sub$cluster == cl
            if (is.null(dim(intensity.sub))) {
              # Only one spot in cluster
              as.numeric(intensity.sub[mask])
            } else {
              colMeans(as.matrix(intensity.sub[mask, , drop = FALSE]), na.rm = TRUE)
            }
          })
          if (!is.matrix(pb)) pb <- matrix(pb, ncol = length(cluster_ids))
          colnames(pb) <- paste0(sample_id, ":C", cluster_ids)
          rownames(pb) <- colnames(intensity.sub)
          all_pseudobulk[[sample_id]] <- pb
        }
        cluster_result(all_results)
        pseudobulk_result(all_pseudobulk)
      })
    })

    cluster_plot <- shiny::reactive({
      meta_list <- cluster_result()
      shiny::req(meta_list)
      # If only one sample, plot as before
      if (is.data.frame(meta_list)) {
        meta <- meta_list
        sid <- unique(meta$Sample)
        meta$unique_cluster <- paste0(sid, '_', as.character(meta$cluster))
        ggplot2::ggplot(meta, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data$unique_cluster, fill = .data$unique_cluster)) +
          ggplot2::geom_tile() +
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
      } else if (is.list(meta_list)) {
        # Plot all samples in a grid
        plots <- lapply(names(meta_list), function(sid) {
          meta <- meta_list[[sid]]
          meta$unique_cluster <- paste0(sid, '_', as.character(meta$cluster))
          ggplot2::ggplot(meta, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data$unique_cluster, fill = .data$unique_cluster)) +
            ggplot2::geom_tile() +
            ggplot2::theme_classic() +
            ggplot2::ggtitle(sid) +
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
        })
        patchwork::wrap_plots(plots, nrow = max(1,floor(sqrt(length(input[['sampleToCluster']]))))) +
          patchwork::plot_annotation(title = "Spatial Clustering Results")
      }
    })
    output[['plotClusters']] <- shiny::renderPlot({
      cluster_plot()
    })

    # Disable combined clustering UI if only one sample selected
    output$pca_or_disabled <- shiny::renderUI({
      if (length(input$sampleToCluster) == 1) {
        shiny::tags$div(style = 'color: #888; margin-top: 40px;', 'Combined PCA is only available when multiple samples are selected.')
      } else {
        shiny::plotOutput(session$ns('pseudobulkPCA'))
      }
    }) |> shiny::bindEvent(input$run_clustering)

    output$heatmap_or_disabled <- shiny::renderUI({
      if (length(input$sampleToCluster) == 1) {
        shiny::tags$div(style = 'color: #888; margin-top: 40px;', 'Combined heatmap is only available when multiple samples are selected.')
      } else {
        plotly::plotlyOutput(session$ns('pseudobulkHeatmap'))
      }
    }) |> shiny::bindEvent(input$run_clustering)
    output$combined_spatial_or_disabled <- shiny::renderUI({
      if (length(input$sampleToCluster) == 1) {
        shiny::tags$div(style = 'color: #888; margin-top: 40px;', 'Combined spatial plot is only available when multiple samples are selected.')
      } else {
        shiny::plotOutput(session$ns('combinedSpatial'), height = '600px')
      }
    }) |> shiny::bindEvent(input$run_clustering)

    # Add clusters to shared_data with user-provided name
    # Add sample-specific clusters to shared_data
    shiny::observeEvent(input$add_to_shared, {
      meta <- cluster_result()
      cluster_name <- input$clusterName
      shiny::req(meta, cluster_name, nzchar(cluster_name))
      if (!is.null(shared_data)) {
        shared_data$updated.metadata <- shiny::isolate({
          m <- shared_data$updated.metadata
          sample_col <- colnames(bulk.metadata)[1]
          # Handle both single and multi-sample cases
          if (is.list(meta)) {
            for (sid in names(meta)) {
              mask <- m[, sample_col] == sid
              if (!cluster_name %in% colnames(m)) m[[cluster_name]] <- NA
              m[[cluster_name]][mask] <- as.character(meta[[sid]]$cluster)
            }
          } else {
            mask <- m[, sample_col] == meta$Sample
            if (!cluster_name %in% colnames(m)) m[[cluster_name]] <- NA
            m[[cluster_name]][mask] <- as.character(meta$cluster)
          }
          m
        })
      }
    })

    # Add combined clusters to shared_data
    shiny::observeEvent(input$add_combined_to_shared, {
      meta <- cluster_result()
      pb_list <- pseudobulk_result()
      comb <- combined_cluster_assignments()
      cluster_name <- input$combinedClusterName
      shiny::req(meta, pb_list, comb, cluster_name, nzchar(cluster_name))
      if (!is.null(shared_data)) {
        shared_data$updated.metadata <- shiny::isolate({
          m <- shared_data$updated.metadata
          sample_col <- colnames(bulk.metadata)[1]
          pb_mat <- do.call(cbind, pb_list)
          cluster_names <- colnames(pb_mat)
          comb_map <- stats::setNames(comb, cluster_names)
          comb_colname <- cluster_name
          if (!comb_colname %in% colnames(m)) m[[comb_colname]] <- NA
          # Map each spot to combined cluster
          if (is.list(meta)) {
            for (sid in names(meta)) {
              spot_cluster_names <- paste0(sid, ":C", as.character(meta[[sid]]$cluster))
              mask <- m[, sample_col] == sid
              m[[comb_colname]][mask] <- as.character(comb_map[spot_cluster_names])
            }
          } else {
            spot_cluster_names <- paste0(meta$Sample, ":C", as.character(meta$cluster))
            mask <- m[, sample_col] == meta$Sample
            m[[comb_colname]][mask] <- as.character(comb_map[spot_cluster_names])
          }
          m
        })
      }
    })

            # Helper: get unique cluster names per sample
        get_unique_cluster_names <- function(meta_list) {
          # Returns a list: sample -> vector of unique cluster names (e.g. sample1_1, sample1_2, ...)
          lapply(names(meta_list), function(sid) {
            cl <- as.character(meta_list[[sid]]$cluster)
            paste0(sid, '_', cl)
          })
        }

        pseudobulk_preproc_pca <- shiny::reactive({
          pb_list <- pseudobulk_result()
          comb <- combined_cluster_assignments()
          shiny::req(pb_list, comb)
          pb_mat <- do.call(cbind, pb_list)
          pb_mat <- pb_mat[matrixStats::rowVars(pb_mat) > 0, ]
          pca <- stats::prcomp(t(pb_mat), scale. = TRUE)
          pc_df <- as.data.frame(pca$x)
          pc_df$cluster <- colnames(pb_mat)
          pc_df$combined_cluster <- as.factor(comb[pc_df$cluster])
          ggplot2::ggplot(pc_df, ggplot2::aes(x = .data$PC1, y = .data$PC2, label = .data$cluster, color = .data$combined_cluster)) +
            ggplot2::geom_point(size = 3) +
            ggplot2::geom_text(size = 3, vjust = 1.5) +
            ggplot2::theme_classic() +
            ggplot2::labs(title = 'PCA of pseudobulked clusters', color = 'Combined cluster') +
            ggplot2::scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(pc_df$combined_cluster)), "Dark2"))

        })
        # PCA of pseudobulked clusters, coloured by combined cluster
        output$pseudobulkPCA <- shiny::renderPlot({
          pseudobulk_preproc_pca()
        })

        # Helper: compute similarity matrix and hc for pseudobulked clusters
        get_pseudobulk_hc <- shiny::reactive({
          pb_list <- pseudobulk_result()
          shiny::req(pb_list)
          pb_mat <- do.call(cbind, pb_list)
          pb_mat <- pb_mat[matrixStats::rowVars(pb_mat) > 0, ]
          sim_mat <- stats::cor(pb_mat, use = "pairwise.complete.obs", method = "pearson")
          d <- stats::as.dist(1 - sim_mat)
          hc <- stats::hclust(d)
          list(sim_mat = sim_mat, hc = hc)
        })

        # Similarity heatmap of pseudobulked clusters (correlation matrix), coloured by combined cluster
        output$pseudobulkHeatmap <- plotly::renderPlotly({
          pb_hc <- get_pseudobulk_hc()
          comb <- combined_cluster_assignments()
          shiny::req(pb_hc, comb)
          sim_mat <- pb_hc$sim_mat
          hc <- pb_hc$hc
          heatmaply::heatmaply(
            sim_mat,
            Colv = stats::as.dendrogram(hc),
            Rowv = stats::as.dendrogram(hc),
            scale = "none",
            labCol = colnames(sim_mat),
            labRow = rownames(sim_mat),
            main = "Similarity (correlation) between clusters",
            RowSideColors = data.frame('Cluster'=RColorBrewer::brewer.pal(length(unique(comb)), "Dark2")[as.numeric(comb)]),
            ColSideColors = data.frame('Cluster'=RColorBrewer::brewer.pal(length(unique(comb)), "Dark2")[as.numeric(comb)])
          )
        })

        # Combined clustering on pseudobulked clusters (use same hc as heatmap)
        combined_cluster_assignments <- shiny::reactive({
          pb_hc <- get_pseudobulk_hc()
          shiny::req(pb_hc)
          hc <- pb_hc$hc
          n_combined <- input$combined_n_clusters
          stats::cutree(hc, k = n_combined)
        })

        combined_spatial_preproc <- shiny::reactive({
                    meta_list <- cluster_result()
          pb_list <- pseudobulk_result()
          comb <- combined_cluster_assignments()
          shiny::req(meta_list, pb_list, comb)
          pb_mat <- do.call(cbind, pb_list)
          pb_mat <- pb_mat[matrixStats::rowVars(pb_mat) > 0, ]
          cluster_names <- colnames(pb_mat)
          comb_map <- stats::setNames(comb, cluster_names)
          meta_all <- do.call(rbind, lapply(names(meta_list), function(sid) {
            meta <- meta_list[[sid]]
            spot_cluster_names <- paste0(sid, ":C", as.character(meta$cluster))
            meta$combined_cluster <- factor(comb_map[spot_cluster_names], levels = 1:input$combined_n_clusters)
            meta$Sample <- sid
            meta
          }))
          ggplot2::ggplot(meta_all, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, fill = .data$combined_cluster, color = .data$combined_cluster)) +
            ggplot2::geom_tile() +
            ggplot2::facet_wrap(~.data$Sample, scales = 'free', nrow = max(1,floor(sqrt(length(input[['sampleToCluster']]))))) +
            ggplot2::theme_classic() +
            ggplot2::labs(title = 'Combined clusters (all samples)', fill = 'Combined cluster') +
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
            ) +
            ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(meta_all$combined_cluster)), "Dark2")) +
            ggplot2::scale_color_manual(values = RColorBrewer::brewer.pal(length(unique(meta_all$combined_cluster)), "Dark2"))
        })
        # Show combined clusters spatially (reactive)
        output$combinedSpatial <- shiny::renderPlot({
          combined_spatial_preproc()
        })

    # Heatmap download
    output[['downloadHeatmap']] <- shiny::downloadHandler(
      filename = function() { input[['heatmapPlotFileName']] },
      content = function(file) {
        pb_hc <- get_pseudobulk_hc()
        shiny::req(pb_hc)
        sim_mat <- pb_hc$sim_mat
        hc <- pb_hc$hc
        comb <- combined_cluster_assignments()
        shiny::req(comb)

        heatmaply::heatmaply(
          sim_mat,
          Colv = stats::as.dendrogram(hc),
          Rowv = stats::as.dendrogram(hc),
          scale = "none",
          labCol = colnames(sim_mat),
          labRow = rownames(sim_mat),
          main = "Similarity (correlation) between clusters",
          showticklabels = c(TRUE, TRUE),
          margins = c(60, 60),
          file = file,
          width = input[['heatmapPlotWidth']]*300,
          height = input[['heatmapPlotHeight']]*300,
          RowSideColors = data.frame('Cluster'=RColorBrewer::brewer.pal(length(unique(comb)), "Dark2")[as.numeric(comb)]),
          ColSideColors = data.frame('Cluster'=RColorBrewer::brewer.pal(length(unique(comb)), "Dark2")[as.numeric(comb)])

        )
      }
    )


    output[['downloadSpatial']] <- create_download_plot_handler(
      plot_func = cluster_plot,
      filename_func = function() { input[['spatialPlotFileName']] },
      width_func = function() { input[['spatialPlotWidth']] },
      height_func = function() { input[['spatialPlotHeight']] }
    )
        
    # PCA download
    output[['downloadPCA']] <- create_download_plot_handler(
      plot_func = pseudobulk_preproc_pca,
      filename_func = function() { input[['pcaPlotFileName']] },
      width_func = function() { input[['pcaPlotWidth']] },
      height_func = function() { input[['pcaPlotHeight']] }
    )

    # Combined spatial download
    output[['downloadCombinedSpatial']] <- create_download_plot_handler(
      plot_func = combined_spatial_preproc,
      filename_func = function() { input[['combinedPlotFileName']] },
      width_func = function() { input[['combinedPlotWidth']] },
      height_func = function() { input[['combinedPlotHeight']] }
    )

})
}

