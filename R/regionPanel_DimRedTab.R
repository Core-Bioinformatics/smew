#' Applies and visualises dimensionality reduction
#'
#' @description UI and server logic for the Dimensionality Reduction panel, enabling users to perform and visualise dimensionality reduction (PCA, NMF, UMAP) on spatial omics data at the region/pixel level. Supports sample selection, component/factor extraction, spatial and UMAP visualisation, and download of results.
#'
#' @details
#' \itemize{
#'   \item{Performs dimensionality reduction (PCA, NMF) on all pixels across selected samples, followed by UMAP projection into 2 dimensions.}
#'   \item{User selects the number of components/factors to compute.}
#'   \item{Dimensionality reduction is only performed when the 'Run dimensionality reduction' button is pressed.}
#'   \item{Distribution of factor/component values by metadata column can be visualised, and top contributing peaks extracted.}
#'   \item{Resulting components/factors can be visualised spatially, and a UMAP can be calculated to visualise the reduction in 2D.}
#'   \item{Pixels in the UMAP can be coloured by metadata or peak intensity.}
#'   \item{Regions can be created based on thresholding dimensionality reductions.}
#'   \item{All results and plots are available for download.}
#' }
#'
#' @rdname RegionPanel_DimRedTab
#' @name RegionPanel_DimRedTab
#' @param id Shiny module id
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param full.metadata Data frame of full sample metadata
#' @param full.intensity.matrix Matrix of intensities (features x samples)
#' @param show Logical; whether to show the panel (default TRUE)
#' @return A shiny tabPanel object for the dimensionality reduction tab
#' @export
RegionPanel_DimRedTabUI <- function(id, bulk.metadata, full.metadata, full.intensity.matrix,show = TRUE){
  ns <- shiny::NS(id)
  # add option to name cluster and retain it
  if(show){
      shiny::tabPanel(
      'Dimensionality Reduction',
      shiny::br(),
      bslib::accordion(
        bslib::accordion_panel(
          title = "Information",
          icon = bsicons::bs_icon("info-circle"),
          open = FALSE,
          shiny::tags$ul(
            shiny::tags$li("This tab will perform dimensionality reduction on all pixels across whichever samples are selected."),
            shiny::tags$li("NMF and PCA are currently available followed by UMAP projection into 2 dimensions."),
            shiny::tags$li("The number of dimensions (i.e. components/factors) computed is user-selected."),
            shiny::tags$li("The dimensionality reduction is only performed once 'Run dimensionality reduction' button is pressed as is a longer-running process."),
            shiny::tags$li("The distribution of factor/component values by metadata column can be visualised and the top contributing peaks extracted."),
            shiny::tags$li("The resulting components/factors can be visualised spatially and a UMAP can further be calculated to visualise the dimensionality reduction in 2 dimensions."),
            shiny::tags$li("Pixels in the UMAP can be coloured according to metadata information or peak intensity."),
            shiny::tags$li("Regions can be created based on thresholding these dimensionality reductions or individual peak intensities.")
          )
        ),
        bslib::accordion_panel(
          title = "Sample selection",
          icon = bsicons::bs_icon("gear"),
          open = TRUE,
          shiny::selectInput(
            inputId = ns("samplesToFactor"),
            label = "Select samples for infer dimensionality reduction",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[,1]),
            multiple = TRUE, width = '100%'
          ),
            shiny::radioButtons(ns("logTransform"), "Log-transform data?", choices = c('Log intensity' = TRUE, 'Raw intensity' = FALSE))
        ),
        bslib::accordion_panel(
          title = "Downloads",
          icon = bsicons::bs_icon("download"),
          open = FALSE,
          shiny::tags$div(
            shiny::tags$h4("Download output plots"),
            shiny::fluidRow(
              shiny::column(6,
                shiny::tags$strong("Dimensionality reduction metadata boxplot"),
                shiny::textInput(ns('perSampleFileName'), 'File name', value ='PerSampleDimRed.png'),
                shiny::numericInput(ns('perSampleWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('perSampleHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadPerSampleDimRed'), 'Download boxplot')
              ),
              shiny::column(6,
                shiny::tags$strong("Spatial dimensionality reduction"),
                shiny::textInput(ns('dimRedSpatialFileName'), 'File name', value ='DimRedSpatial.png'),
                shiny::numericInput(ns('dimRedSpatialWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('dimRedSpatialHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadDimRedSpatial'), 'Download spatial plot')
              )
            ),
            shiny::fluidRow(
              shiny::column(6,
                shiny::tags$strong("Top dimensionality reduction weights"),
                shiny::textInput(ns('dimRedWeightsFileName'), 'File name', value ='DimRedFeatureWeights.png'),
                shiny::numericInput(ns('dimRedWeightsWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('dimRedWeightsHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadDimRedWeights'), 'Download feature weights')
              ),
              shiny::column(6,
                shiny::tags$strong("UMAP scatter plot"),
                shiny::textInput(ns('umapFileName'), 'File name', value ='UMAP.png'),
                shiny::numericInput(ns('umapWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('umapHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadUMAP'), 'Download UMAP')
              )
            ),
            shiny::fluidRow(
              shiny::column(6,
                shiny::tags$strong("Spatial UMAP plot"),
                shiny::textInput(ns('spatialUMAPFileName'), 'File name', value ='SpatialUMAP.png'),
                shiny::numericInput(ns('spatialUMAPWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('spatialUMAPHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadSpatialUMAP'), 'Download spatial UMAP')
              ),
              shiny::column(6,
                shiny::tags$strong("Spatial clusters"),
                shiny::textInput(ns('clusterSpatialFileName'), 'File name', value ='ClusterSpatial.png'),
                shiny::numericInput(ns('clusterSpatialWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('clusterSpatialHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadClusterSpatial'), 'Download cluster spatial')
              )
            ),
            shiny::fluidRow(
              shiny::column(6,
                shiny::tags$strong("Cluster metadataboxplot"),
                shiny::textInput(ns('clusterBoxFileName'), 'File name', value ='ClusterBoxplot.png'),
                shiny::numericInput(ns('clusterBoxWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('clusterBoxHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadClusterBox'), 'Download cluster boxplot')
              )
            )
          )
        ),
        id = ns("acc"),
        open = "Sample selection"
      ),
      shiny::tabsetPanel(
        shiny::tabPanel(
          "Dimensionality Reduction",
          shiny::fluidRow(
            shiny::column(12,
              bslib::card(
                class = "mb-3",
                style = "width:100%;",
                shiny::h4("Dimensionality Reduction Settings"),
                shiny::fluidRow(
                  shiny::column(3,
                    shiny::radioButtons(ns("dimReduction"), "Method:", c('PCA','NMF'))
                  ),
                  shiny::column(3,
                    shiny::numericInput(inputId = ns("numDimensions"),
                      label = 'Number of dimensions',
                      value = 10,
                      min = 2,
                      max = 50,
                      step = 1),
                  shiny::numericInput(ns('seed'),label = 'Seed for reproducibility',value = 23)
                  ),
                  shiny::column(3,
                    shiny::actionButton(
                      inputId = ns("run_dimred"),
                      label = "Run",
                      icon = shiny::icon("play"),
                      style = "margin-top:25px;width:100%;"
                    )
                  ),
                  shiny::column(3,
                    shiny::selectInput(inputId=ns("focus_dimension"),
                      label = 'Dimension of interest',
                      choices = 1:10,selected = 1)
                  )
                )
              )
            )
          ),
          shiny::tags$h4("Spatial Plot"),
          shiny::plotOutput(ns('plotDimRed'),click=ns('dimRed_click'),height = '600px'),
          shiny::tags$h4("Top Contributing Peaks"),
          plotly::plotlyOutput(ns('plotDimRedFeatureWeights')),
          shiny::tags$h4("Dimension Distributions by Metadata"),
          shiny::sidebarLayout(
            shiny::sidebarPanel(
              shiny::selectInput(inputId = ns("groupingMetadataBox"),
                label = "Metadata to group boxplot on",
                choices = colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','x_tf','y_tf'))],
                selected = colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','x_tf','y_tf'))][1],
                multiple = FALSE)
            ),
            shiny::mainPanel(
              shiny::plotOutput(ns('plotPerSampleDimRed'))
            )
          )
        ),
        shiny::tabPanel(
          "UMAP Projection",
          shiny::fluidRow(
            shiny::column(12,
              bslib::card(
                class = "mb-3",
                style = "width:100%;",
                shiny::tags$h4("UMAP Settings"),
                shiny::fluidRow(
                  shiny::column(6,
                    shiny::actionButton(
                      inputId = ns("runUMAP"),
                      label = "Calculate UMAP",
                      icon = shiny::icon("play"),
                      style = "margin-top:25px;width:100%;"
                    )
                  ),
                  shiny::column(6,
                    shiny::selectizeInput(inputId = ns("colourUMAP"),
                      label = "Metadata to colour UMAP",
                      choices = c(),
                      multiple = FALSE)
                  ),
                )
              )
            )
          ),
          shiny::tags$h4("UMAP Projection from Dimensionality Reduction"),
          shiny::plotOutput(ns('DimRedUMAP'),height = '600px'),
          shiny::tags$h4("Spatial View of UMAP"),
          shiny::plotOutput(ns('DimRedUMAPSpatial'),click=ns('umap_click'),height = '600px'),
          shiny::div(style = "margin-top:10px")
        ),
        shiny::tabPanel(
          "Cluster Assignment",
          shiny::sidebarLayout(
            shiny::sidebarPanel(
              shiny::tags$h4("Clustering Settings"),
              shiny::radioButtons(ns('dimToCluster'),label = 'Dimensionality reduction to use for clustering',choices='UMAP'),
              shiny::selectInput(ns('clusterComponent'),label = 'Cluster component',choices = c()),
              shiny::sliderInput(ns('thresholdHigh'),label = 'Percentage of data to categorise as high intensity',value = 25,min = 1,max = 49,step = 1),
              shiny::sliderInput(ns('thresholdLow'),label = 'Percentage of data to categorise as low intensity',value = 25,min = 1,max = 49,step = 1),
              shiny::selectInput(ns('clusterBoxColour'),label = 'Metadata to split clusters on',choices = colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','x_tf','y_tf'))],selected = colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','x_tf','y_tf'))][1]),
              shiny::textInput(ns('clusterName'),label = 'Name for clusters'),
              shiny::actionButton(ns('lock_cluster'),'Add clusters to object')
            ),
            shiny::mainPanel(
              shiny::tags$h4("Cluster Plot"),
              shiny::plotOutput(ns('clusterPlot')),
              shiny::tags$h4("Cluster Boxplot"),
              shiny::plotOutput(ns('clusterBoxPlot'))
            )
          )
        )
      ),
  )
  }
  }

#' @rdname RegionPanel_DimRedTab
#' @param id Shiny module id
#' @param full.intensity.matrix Matrix of intensities (features x samples)
#' @param full.metadata Data frame of full sample metadata
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param anno Data frame of peak annotations
#' @param shared_data Reactive or shared data object
#' @return None; called for side effects in Shiny module
#' @export
RegionPanel_DimRedTabServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, anno, shared_data){

  shiny::moduleServer(id, function(input, output, session){

    shiny::updateSelectizeInput(session, "colourUMAP",
      choices = c(colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','x_tf','y_tf'))],colnames(full.intensity.matrix)),
      selected = colnames(full.metadata)[!(colnames(full.metadata)%in%c('pixel_id','x','y','x_tf','y_tf'))][1],
      server = TRUE)
    shiny::observeEvent(input$numDimensions, {
      shiny::updateSelectInput(session, "focus_dimension", choices=1:input$numDimensions)
    })

    shiny::observeEvent(list(input$runUMAP, input$run_dimred), {

      # React only when either button is clicked
      btn <- c(input$run_dimred,input$runUMAP)

      # Use another input to decide behavior
      mode <- input$dimReduction
      choices = c()
      if (input$run_dimred > 0){
        choices = c(choices,input$dimReduction)
      }
      if (input$runUMAP > 0){
        choices = c(choices,'UMAP')
      }
      shiny::updateRadioButtons(
        session,
        "dimToCluster",
        choices = choices,
        selected = choices[1]
      )
    })

    shiny::observeEvent(input$dimToCluster, {
      if (input$dimToCluster == 'UMAP'){
        choices = paste0('UMAP',1:2)
      } else {
        choices = paste0(input$dimToCluster,'_',seq(1,input$numDimensions))
      }
      shiny::updateSelectInput(session, 'clusterComponent',choices = choices)

    })


    get_subset_exp <- shiny::reactive({
      idx <- full.metadata$Sample %in% input[['samplesToFactor']]
      current.intensity.matrix <- full.intensity.matrix[idx, ]
      if (input[['logTransform']]) {
        current.intensity.matrix <- log2(current.intensity.matrix + 1)
      }
      current.metadata <- full.metadata[idx, ]
      return(list('exp' = current.intensity.matrix, 'meta' = current.metadata))
    })


    shiny::observeEvent(input$dimRed_click, {
      ns <- session$ns
      shiny::showModal(
        shiny::modalDialog(
          shiny::plotOutput(ns("plotDimRedZoom")),
          easyClose = TRUE,
          footer = NULL
        )
      )
    })

    get_nmf <- shiny::reactive({
      shiny::withProgress(message = 'Computing NMF...', value = 0, {
        dat <- get_subset_exp()
        current.intensity.matrix <- dat$exp
        current.metadata <- dat$meta
        shiny::incProgress(0.1)

        # Remove features with zero variance
        var_vec <- matrixStats::colVars(current.intensity.matrix)
        keep <- var_vec > 0
        filtered_matrix <- current.intensity.matrix[, keep, drop = FALSE]

        shiny::incProgress(0.1)

        # Create SCE object
        spe <- SingleCellExperiment::SingleCellExperiment(
          assays = list(counts = t(filtered_matrix)),
          colData = current.metadata
        )

        shiny::incProgress(0.1)

        # Run NMF externally
        nmf.factors <- RcppML::nmf(as.matrix(filtered_matrix), seed = input[['seed']], k = input[['numDimensions']])
        shiny::incProgress(0.4)

        # Add NMF to reducedDims
        nmf_embeddings <- nmf.factors$w
        colnames(nmf_embeddings) <- paste0('NMF_', seq_len(ncol(nmf_embeddings)))
        SingleCellExperiment::reducedDims(spe)[['NMF']] <- nmf_embeddings

        # Store feature weights (loadings)
        nmf.factor.features <- data.frame(nmf.factors$h)
        colnames(nmf.factor.features) <- colnames(filtered_matrix)
        rownames(nmf.factor.features) <- paste0('NMF_', seq_len(nrow(nmf.factor.features)))

        shiny::incProgress(0.2)
        message('Computed NMF factors')

        return.list <- list(
          'spe_obj' = spe,
          'feature_weights' = nmf.factor.features,
          'dimRed' = input[['dimReduction']],
          'numDimensions' = input[['numDimensions']]
        )
        return(return.list)
      })
    })  |> shiny::bindEvent(input[["run_dimred"]])

    get_pca <- shiny::reactive({
      shiny::withProgress(message = 'Running PCA...', value = 0, {
        message('Running PCA')
        dat <- get_subset_exp()
        current.intensity.matrix <- dat$exp
        current.metadata <- dat$meta
        shiny::incProgress(0.1)

        # Remove features with zero variance
        var_vec <- matrixStats::colVars(current.intensity.matrix)
        keep <- var_vec > 0
        filtered_matrix <- current.intensity.matrix[, keep, drop = FALSE]

        shiny::incProgress(0.1)

        # Create SCE object
        spe <- SingleCellExperiment::SingleCellExperiment(
          assays = list(counts = t(filtered_matrix)),
          colData = current.metadata
        )

        shiny::incProgress(0.1)

        # Run PCA using scater
        spe <- scater::runPCA(
          spe,
          ncomponents = input[['numDimensions']],
          exprs_values = "counts",
          name = "PCA",
          scale = TRUE,
          ntop = ncol(filtered_matrix)
        )

        shiny::incProgress(0.5)

        # Extract feature weights (loadings)
        pca_attr <- attributes(SingleCellExperiment::reducedDims(spe)[['PCA']])
        if (!is.null(pca_attr$rotation)) {
          pc.component.features <- data.frame(t(pca_attr$rotation))
          colnames(pc.component.features) <- colnames(filtered_matrix)
          rownames(pc.component.features) <- paste0('PCA_', seq_len(nrow(pc.component.features)))
        } else {
          # Fallback if rotation not stored
          pc.component.features <- matrix(NA, nrow = input[['numDimensions']], ncol = ncol(filtered_matrix))
          pc.component.features <- data.frame(pc.component.features)
          colnames(pc.component.features) <- colnames(filtered_matrix)
          rownames(pc.component.features) <- paste0('PCA_', seq_len(input[['numDimensions']]))
        }

        message('PCA complete')
        shiny::incProgress(0.2)

        return.list <- list(
          'spe_obj' = spe,
          'feature_weights' = pc.component.features,
          'dimRed' = input[['dimReduction']],
          'numDimensions' = input[['numDimensions']]
        )
        return(return.list)
      })
    })  |> shiny::bindEvent(input[["run_dimred"]])

    get_dimred <- shiny::reactive({
      if (input[["dimReduction"]]=="NMF"){
        return(get_nmf())
      } else {
        return(get_pca())
      }
    })  |> shiny::bindEvent(input[["run_dimred"]])

    nmf_plot <- shiny::reactive({
      obj <- get_dimred()
      spe <- obj$spe_obj

      # Extract embeddings from SCE
      embeddings <- SingleCellExperiment::reducedDims(spe)[[obj$dimRed]]
      current.metadata <- as.data.frame(SingleCellExperiment::colData(spe))
      current.metadata[[paste0(obj$dimRed, '_', input[['focus_dimension']])]] <- embeddings[, as.numeric(input[['focus_dimension']])]

      message('Computed dimensional reduction metadata')
      if (get_dimred()$dimRed=='PCA'){
      my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data[[paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])]], fill = .data[[paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])]])) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low='darkblue',high='darkred',mid='white',name=paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))+
        ggplot2::scale_color_gradient2(low='darkblue',high='darkred',mid='white',name=paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))+
        ggplot2::facet_wrap(~current.metadata$Sample,
                            nrow = max(1,floor(sqrt(length(unique(current.metadata$Sample))/1.5))), scales = 'free') +
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
      } else {
        my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data[[paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])]], fill = .data[[paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])]])) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient2(low='darkblue',high='brown',mid='lightgrey',name=paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))+
          ggplot2::scale_color_gradient2(low='darkblue',high='brown',mid='lightgrey',name=paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))+
          ggplot2::facet_wrap(~current.metadata$Sample,
                              nrow = max(1,floor(sqrt(length(unique(current.metadata$Sample))/1.5))), scales = 'free') +
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
      }
      return(list('plot'=my_plot,'upper_lim'=max(current.metadata[,paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])]),'lower_lim'=min(current.metadata[,paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])])))
    })

    nmf_plot_zoom <- shiny::reactive({
      obj <- get_dimred()
      spe <- obj$spe_obj

      # Extract embeddings from SCE
      embeddings <- SingleCellExperiment::reducedDims(spe)[[obj$dimRed]]
      current.metadata <- as.data.frame(SingleCellExperiment::colData(spe))
      current.metadata[[paste0(obj$dimRed, '_', input[['focus_dimension']])]] <- embeddings[, as.numeric(input[['focus_dimension']])]

      full_plot <- nmf_plot()
      current.metadata <- current.metadata[current.metadata$Sample == input$dimRed_click$panelvar1, ]

      message('Generating zoomed dimensional reduction plot')

      dim_col <- paste0(obj$dimRed, '_', input[['focus_dimension']])
      if (get_dimred()$dimRed=='PCA'){
      my_plot <- ggplot2::ggplot(current.metadata, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data[[dim_col]], fill = .data[[dim_col]])) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = 'darkblue', high = 'darkred', mid = 'white', name = dim_col, limits = c(full_plot$lower_lim, full_plot$upper_lim)) +
        ggplot2::scale_color_gradient2(low = 'darkblue', high = 'darkred', mid = 'white', name = dim_col, limits = c(full_plot$lower_lim, full_plot$upper_lim)) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank()) +
        ggplot2::theme(aspect.ratio = 1)
      } else {
        my_plot <- ggplot2::ggplot(current.metadata, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data[[dim_col]], fill = .data[[dim_col]])) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_gradient2(low = 'darkblue', high = 'brown', mid = 'lightgrey', name = dim_col, limits = c(full_plot$lower_lim, full_plot$upper_lim)) +
          ggplot2::scale_color_gradient2(low = 'darkblue', high = 'brown', mid = 'lightgrey', name = dim_col, limits = c(full_plot$lower_lim, full_plot$upper_lim)) +
          ggplot2::theme_classic() +
          ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                         axis.text.x = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(),
                         axis.line.x = ggplot2::element_blank(),
                         axis.title.y = ggplot2::element_blank(),
                         axis.text.y = ggplot2::element_blank(),
                         axis.ticks.y = ggplot2::element_blank(),
                         axis.line.y = ggplot2::element_blank()) +
          ggplot2::theme(aspect.ratio = 1)
      }
      return(my_plot)
    })

    nmf_persample <- shiny::reactive({
      obj <- get_dimred()
      spe <- obj$spe_obj

      # Extract embeddings from SCE
      embeddings <- SingleCellExperiment::reducedDims(spe)[[obj$dimRed]]
      current.metadata <- as.data.frame(SingleCellExperiment::colData(spe))
      current.metadata[[paste0(obj$dimRed, '_', input[['focus_dimension']])]] <- embeddings[, as.numeric(input[['focus_dimension']])]

      message('Computed per-sample dimensional reduction metadata')

      dim_col <- paste0(obj$dimRed, '_', input[['focus_dimension']])
      my_plot <- ggplot2::ggplot(current.metadata, ggplot2::aes(y = .data[[dim_col]], x = .data[[input[['groupingMetadataBox']]]], fill = .data[[input[['groupingMetadataBox']]]])) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggplot2::ylab(dim_col) +
        ggplot2::scale_fill_discrete(name = input[['groupingMetadataBox']]) +
        ggplot2::xlab(input[['groupingMetadataBox']])
      return(my_plot)
    })

    nmf_featureweights <- shiny::reactive({
      pca = get_dimred()
      nmf.weights = data.frame(t(pca$feature_weights))
      nmf.weights = nmf.weights[order(-abs(nmf.weights[,paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])])),]
      nmf.weights = utils::head(nmf.weights,20)
      nmf.weights$m_z = rownames(nmf.weights)
      nmf.weights$m_z = factor(nmf.weights$m_z,levels=rev(nmf.weights$m_z))
      nmf.weights = merge(nmf.weights,anno)
      nmf.weights$name = stringr::str_wrap(nmf.weights$name,width=30)
      ggplot2::ggplot(nmf.weights,ggplot2::aes(y=.data$m_z,x=.data[[paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])]],label=.data$name))+
        ggplot2::geom_bar(stat='identity')+
        ggplot2::theme_classic()+
        ggplot2::xlab(paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))+
        ggplot2::theme(aspect.ratio = 1)
    })

    nmf_umap_prep <- shiny::reactive({
      shiny::withProgress(message = 'Calculating UMAP...', value = 0, {
        obj <- get_dimred()
        spe <- obj$spe_obj

        shiny::incProgress(0.2)

        # Run UMAP using scater on the dimensionality reduction
        spe <- scater::runUMAP(
          spe,
          dimred = obj$dimRed,
          name = "UMAP",
          n_neighbors = 30,
          n_dimred = obj$numDimensions
        )

        shiny::incProgress(0.5)

        # Extract UMAP and metadata
        umap_embeddings <- SingleCellExperiment::reducedDims(spe)[['UMAP']]
        current.metadata <- as.data.frame(SingleCellExperiment::colData(spe))
        current.metadata$UMAP_1 <- umap_embeddings[, 1]
        current.metadata$UMAP_2 <- umap_embeddings[, 2]

        # Rescale for color blending
        current.metadata$UMAP_1_scaled <- scales::rescale(current.metadata$UMAP_1)
        current.metadata$UMAP_2_scaled <- scales::rescale(current.metadata$UMAP_2)
        colours <- region_utils_colour_blender(current.metadata[, c('UMAP_1_scaled', 'UMAP_2_scaled')], channels.use = c("red", "blue"))
        current.metadata$my.colour <- colours

        shiny::incProgress(0.3)
        return(list('metadata' = current.metadata, 'spe_obj' = spe))
      })
    }) |> shiny::bindEvent(input[['runUMAP']])

    nmf_umap <- shiny::reactive({
      umap_obj <- nmf_umap_prep()
      current.metadata <- umap_obj$metadata

      if (input[['colourUMAP']] %in% colnames(get_subset_exp()$exp)){
        current.metadata$selectedPeak = get_subset_exp()$exp[,input[['colourUMAP']]]
        current.metadata$selectedPeak = pmin(stats::quantile(current.metadata$selectedPeak,0.95),current.metadata$selectedPeak)
        current.metadata$selectedPeak = pmax(stats::quantile(current.metadata$selectedPeak,0.05),current.metadata$selectedPeak)
        current.metadata = current.metadata[order(current.metadata$selectedPeak),]
        colour.function = ggplot2::scale_color_gradient(name=input[['colourUMAP']],low = "lightgrey", high = "brown")

      } else {
        current.metadata$selectedPeak = current.metadata[,input[['colourUMAP']]]
        colour.function = ggplot2::scale_color_discrete(name=input[['colourUMAP']])
      }
      return(list('SpatialView'=ggplot2::ggplot(current.metadata,ggplot2::aes(x=.data$x_tf,y=.data$y_tf, color=.data$my.color,fill=.data$my.color))+
                    ggplot2::geom_tile()+
                    ggplot2::theme_classic() +
                    ggplot2::scale_color_identity() +
                    ggplot2::scale_fill_identity() +
                    ggplot2::facet_wrap(~current.metadata$Sample,
                               nrow = max(1,floor(sqrt(length(unique(current.metadata$Sample))/1.5))),
                               scales='free')+
                    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                                   axis.text.x=ggplot2::element_blank(),
                                   axis.ticks.x=ggplot2::element_blank(),
                                   axis.line.x = ggplot2::element_blank(),
                                   axis.title.y=ggplot2::element_blank(),
                                   axis.text.y=ggplot2::element_blank(),
                                   axis.ticks.y=ggplot2::element_blank(),
                                   axis.line.y = ggplot2::element_blank())+
                    ggplot2::theme(aspect.ratio = 1),
                  'UMAP'=ggplot2::ggplot(current.metadata,ggplot2::aes(x=.data$UMAP_1,y=.data$UMAP_2,color=.data$selectedPeak))+
                    ggplot2::geom_point()+
                    ggplot2::theme_classic()+
                    colour.function+
                    ggplot2::theme(aspect.ratio = 1)))

    })

    nmf_umap_zoom <- shiny::reactive({
      # Guard: require a valid panel click context
      if (is.null(input$umap_click) || is.null(input$umap_click$panelvar1)) return(NULL)
      # Get full UMAP metadata
      umap_obj <- nmf_umap_prep()
      current.metadata <- umap_obj$metadata
      # Subset to the clicked sample
      current.metadata <- current.metadata[current.metadata$Sample == input$umap_click$panelvar1, , drop = FALSE]
      # Recalculate UMAP colors for the subset
      if (nrow(current.metadata) > 0) {
        current.metadata$UMAP_1_scaled <- scales::rescale(current.metadata$UMAP_1)
        current.metadata$UMAP_2_scaled <- scales::rescale(current.metadata$UMAP_2)
        current.metadata$my.colour <- region_utils_colour_blender(current.metadata[,c('UMAP_1_scaled','UMAP_2_scaled')], channels.use = c("red","blue"))
      }
      # Recalculate selectedPeak and colour function for the subset using safe alignment by pixel_id
      if (input[['colourUMAP']] %in% colnames(get_subset_exp()$exp)) {
        exp_mat <- get_subset_exp()$exp
        # Prefer matching by pixel_id to avoid rowname-based out-of-bounds
        if ('pixel_id' %in% colnames(current.metadata) && !is.null(rownames(exp_mat))) {
          idx <- match(current.metadata$pixel_id, rownames(exp_mat))
        } else {
          # Fallback to rownames if available
          idx <- match(rownames(current.metadata), rownames(exp_mat))
        }
        sel_vec <- rep(NA_real_, nrow(current.metadata))
        valid <- !is.na(idx) & idx > 0 & idx <= nrow(exp_mat)
        if (any(valid)) {
          sel_vec[valid] <- exp_mat[idx[valid], input[['colourUMAP']], drop = TRUE]
        }
        current.metadata$selectedPeak <- sel_vec
        if (sum(!is.na(current.metadata$selectedPeak)) > 0) {
          q95 <- stats::quantile(current.metadata$selectedPeak, 0.95, na.rm = TRUE)
          q05 <- stats::quantile(current.metadata$selectedPeak, 0.05, na.rm = TRUE)
          current.metadata$selectedPeak <- pmin(q95, current.metadata$selectedPeak)
          current.metadata$selectedPeak <- pmax(q05, current.metadata$selectedPeak)
          current.metadata <- current.metadata[order(current.metadata$selectedPeak), , drop = FALSE]
        }
        colour.function <- ggplot2::scale_color_gradient(name=input[['colourUMAP']], low = "lightgrey", high = "brown")
      } else {
        current.metadata$selectedPeak <- current.metadata[,input[['colourUMAP']]]
        colour.function <- ggplot2::scale_color_discrete(name=input[['colourUMAP']])
      }
      ggplot2::ggplot(current.metadata, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data$my.colour, fill = .data$my.colour)) +
        ggplot2::geom_tile() +
        ggplot2::theme_classic() +
        ggplot2::scale_color_identity() +
        ggplot2::scale_fill_identity() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank()) +
        ggplot2::theme(aspect.ratio = 1)
    })

    get_clusters <- shiny::reactive({
      if (input$dimToCluster == 'UMAP') {
        umap_obj <- nmf_umap_prep()
        clusterdf <- umap_obj$metadata
        my_component <- input$clusterComponent
      } else if (input$dimToCluster %in% c('PCA', 'NMF')) {
        obj <- get_dimred()
        spe <- obj$spe_obj
        embeddings <- SingleCellExperiment::reducedDims(spe)[[obj$dimRed]]
        clusterdf <- as.data.frame(SingleCellExperiment::colData(spe))
        # Add all dimension columns
        for (i in seq_len(ncol(embeddings))) {
          clusterdf[[paste0(obj$dimRed, '_', i)]] <- embeddings[, i]
        }
        my_component <- input$clusterComponent
      } else {
        message(paste0('Clustering by selected peak: ', input$selectedPeak))
        clusterdf = cbind(get_subset_exp()$meta,
                          get_subset_exp()$exp[,input$selectedPeak,drop=F])
        my_component = input$selectedPeak
      }
      my_component_level <- clusterdf[, my_component]
      quantiles <- stats::quantile(my_component_level, prob = c(input[['thresholdLow']] / 100, (100 - input[['thresholdHigh']]) / 100), type = 1)
      clusterdf$cluster <- ifelse(my_component_level <= quantiles[1], 'Low',
                     ifelse(my_component_level >= quantiles[2], 'High', 'Medium'))
        return(clusterdf)
      })

    cluster_plot <- shiny::reactive({
      clusterdf = get_clusters()
      ggplot2::ggplot(clusterdf,ggplot2::aes(x=.data$x_tf,y=.data$y_tf,fill=.data$cluster,color=.data$cluster)) +
        ggplot2::geom_tile() +
        ggplot2::facet_wrap(ggplot2::vars(clusterdf$Sample),
                            nrow = max(1,floor(sqrt(length(unique(clusterdf$Sample))/1.5))), scales = 'free') +
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

    cluster_boxplot <- shiny::reactive({
      clusterdf = get_clusters()
      ggplot2::ggplot(clusterdf,ggplot2::aes(x=.data[[input$clusterBoxColour]],fill=.data$cluster))+
        ggplot2::geom_bar(position='fill') + ggplot2::theme_classic()

    })
    output[['plotDimRed']] <- shiny::renderPlot({
      nmf_plot()$plot
    })

    output[['clusterPlot']] <- shiny::renderPlot({
      cluster_plot()
    })

    output[['clusterBoxPlot']] <- shiny::renderPlot({
      cluster_boxplot()
    })

    output[['plotDimRedZoom']] <- shiny::renderPlot({
      nmf_plot_zoom()
    })

    output[['plotDimRedFeatureWeights']] <- plotly::renderPlotly(
      nmf_featureweights()
    )
    output[['plotPerSampleDimRed']] <- shiny::renderPlot(
      nmf_persample()
    )
    # output[['plotPerSampleDimRed']] <- plotly::renderPlotly(
    #   plotly::ggplotly(nmf_persample()) |>
    #     layout(boxmode = "group")
#    )

    output[['DimRedUMAPSpatial']] <- shiny::renderPlot({
      nmf_umap()$SpatialView
    })

    output[['DimRedUMAPSpatialZoom']] <- shiny::renderPlot({
      nmf_umap_zoom()
    })

    output[['downloadSpatialUMAP']] <- utils_create_download_plot_handler(
      plot_func = function() nmf_umap()$SpatialView,
      filename_func = function() input[['spatialUMAPFileName']],
      width_func = function() input[['spatialUMAPWidth']],
      height_func = function() input[['spatialUMAPHeight']]
    )

    output[['DimRedUMAP']] <- shiny::renderPlot(
      nmf_umap()$UMAP
    )

    output[['downloadUMAP']] <- utils_create_download_plot_handler(
      plot_func = function() nmf_umap()$SpatialView,
      filename_func = function() input[['umapFileName']],
      width_func = function() input[['umapWidth']],
      height_func = function() input[['umapHeight']]
    )

    shiny::observeEvent(input$umap_click, {
      ns <- session$ns
      shiny::showModal(
        shiny::modalDialog(
          shiny::plotOutput(ns('DimRedUMAPSpatialZoom')),
          easyClose = TRUE,
          footer = NULL
        )
      )
    })

    shiny::observeEvent(input$lock_cluster, {
      shiny::withProgress(message = 'Applying cluster assignment...', value = 0, {
        message('Locking cluster assignment')
        meta <- shared_data$updated.metadata
        clusterdf = get_clusters()
        shiny::incProgress(0.3)
        message(paste0('meta rows: ', nrow(meta), '; clusterdf rows: ', nrow(clusterdf)))
        message(paste0('cluster name: ', input$clusterName))
        if (!input$clusterName %in% colnames(meta)) {
          meta[[input$clusterName]] <- NA
        }
        idx <- match(meta$pixel_id, clusterdf$pixel_id)
        to_update <- !is.na(idx)
        meta[[input$clusterName]][to_update] <- clusterdf$cluster[idx[to_update]]
        shiny::incProgress(0.6)
        shared_data$updated.metadata <- shiny::isolate(meta)   # update shared dataset
      })
    })

    output[['downloadPerSampleDimRed']] <- utils_create_download_plot_handler(
      plot_func = nmf_persample,
      filename_func = function() input[['perSampleFileName']],
      width_func = function() input[['perSampleWidth']],
      height_func = function() input[['perSampleHeight']]
    )
    output[['downloadDimRedSpatial']] <- utils_create_download_plot_handler(
      plot_func = function() nmf_plot()$plot,
      filename_func = function() input[['dimRedSpatialFileName']],
      width_func = function() input[['dimRedSpatialWidth']],
      height_func = function() input[['dimRedSpatialHeight']]
    )
    output[['downloadDimRedWeights']] <- utils_create_download_plot_handler(
      plot_func = nmf_featureweights,
      filename_func = function() input[['dimRedWeightsFileName']],
      width_func = function() input[['dimRedWeightsWidth']],
      height_func = function() input[['dimRedWeightsHeight']]
    )
    output[['downloadClusterSpatial']] <- utils_create_download_plot_handler(
      plot_func = cluster_plot,
      filename_func = function() input[['clusterSpatialFileName']],
      width_func = function() input[['clusterSpatialWidth']],
      height_func = function() input[['clusterSpatialHeight']]
    )
    output[['downloadClusterBox']] <- utils_create_download_plot_handler(
      plot_func = cluster_boxplot,
      filename_func = function() input[['clusterBoxFileName']],
      width_func = function() input[['clusterBoxWidth']],
      height_func = function() input[['clusterBoxHeight']]
    )

  })
}


