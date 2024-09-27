#' @rdname RegionSpatialClusterPanel
#' @export
RegionSpatialClusterPanelUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- NS(id)
  # add option to name cluster and retain it
  if(show){
    tabPanel(
      'Spatial Clustering',
      sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(
          # samples to include
          # number of clusters
          # any other parameters
          dropMenu(
            circleButton(ns("info_clustering"), icon = icon("info"),status = "success"),
            tags$div(
              tags$h3("Clustering"),
              tags$ul(
                tags$li("Here you can apply spatially-resolved clustering methods to individual samples. The methodd currently implemented are BayesSpace and a manual lasso user-selection of points but we hope to add more approaches."),
                tags$li("For BayesSpace You need to select the number of clusters you want to select as well as the number of iterations to run the model for. The authors suggest using at least 10,000 iterations but we have set the default to a low number here for the sake of runtime. You can then set the burn-in period, which is a number of samples to exclude at the beginning before the Markov Chain may have reached its equilibrium distribution."),
                tags$li("You can then further smooth these clusters then find markers for the different clusters."),
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
          selectInput(
            inputId = ns("sampleToCluster"),
            label = "Select samples for clustering",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[1,1]),
            multiple = FALSE
          ),
          selectInput(inputId = ns("clusteringApproach"),
                      "Clustering method",
                      choices = c("BayesSpace","Lasso selection")),
          conditionalPanel(
            id = ns('bayesSpaceOptions'),
            ns=ns,
            condition = "input.clusteringApproach == 'BayesSpace'",
            selectInput(ns('numClusters'), 'Number of spatial clusters to infer', multiple = FALSE, choices = 2:10),
            numericInput(ns('nRep'),'The number of MCMC iterations for BayesSpace',min = 10,max = 10000,value = 100,step = 10),
            numericInput(ns('burnIn'),'The number of MCMC iterations to exclude as burn-in period for BayesSpace (must be less than number of iterations above)',min = 1,max = 1000,value = 10,step = 10),
            # button to start clustering
            tags$p("Warning: this analysis can take a few minutes to run."),
            actionButton(
              inputId = ns("run_clustering"),
              label = "Find clusters",
              icon = icon("play")
            ),
          ),
          # conditionalPanel(
          #   id = ns('lassoOptions'),
          #   ns=ns,
          #   condition = "input.clusteringApproach == 'Lasso selection'"
          # ),

          div(style = "margin-top:10px"),
          dropMenu(
            circleButton(ns("downloads"), icon = icon("download"),status = "success"),
            tags$div(
              tags$h3("Downloads"),
              fluidRow(
                column(5,offset=0,
                       tags$h4("Spatial visualisation"),
                       textInput(ns('spatialPlotFileName'),'File name for download', value ='spatialClusters.png', placeholder = 'spatialClusters.png'),
                       numericInput(ns('spatialPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('spatialPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadSpatial'), 'Download clusters')),
                column(5,offset=0,
                       tags$h4("Smoothed spatial visualisation"),
                       textInput(ns('spatialSmoothPlotFileName'),'File name for download', value ='spatialSmoothClusters.png', placeholder = 'spatialSmoothClusters.png'),
                       numericInput(ns('spatialSmoothPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('spatialSmoothPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadSmoothSpatial'), 'Download clusters'))),
              # fluidRow(
              #   column(5,offset=0,
              #          tags$h4("Bar plot"),
              #          textInput(ns('barPlotFileName'),'File name for download', value ='bar.png', placeholder = 'bar.png'),
              #          numericInput(ns('barPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              #          numericInput(ns('barPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
              #          downloadButton(ns('downloadBarPlot'), 'Download barplot')),
              #
              #   column(5,offset=0,
              #          tags$h4("Box plot"),
              #          textInput(ns('boxPlotFileName'),'File name for download', value ='box.png', placeholder = 'box.png'),
              #          numericInput(ns('boxPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              #          numericInput(ns('boxPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
              #          downloadButton(ns('downloadBoxPlot'), 'Download boxplot')),
              #   ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )),



        ),

        #Main panel for displaying table of enriched pathways
        mainPanel(
          conditionalPanel(
            id = ns('bayesSpaceShowClusters'),
            ns=ns,
            condition = "input.clusteringApproach == 'BayesSpace'",

            plotOutput(ns('plotClusters'))),
          conditionalPanel(
            id = ns('lassoShow'),
            ns=ns,
            condition = "input.clusteringApproach == 'Lasso selection'",
            plotly::plotlyOutput(ns('lassoSelection'),height=800))
      )
    ),
    conditionalPanel(
      id = ns('bayesSpaceSmooth'),
      ns=ns,
      condition = "input.clusteringApproach == 'BayesSpace'",
      sidebarLayout(
        # Sidebar panel for inputs ----
        sidebarPanel(
          actionButton(
            inputId = ns("run_smoothing"),
            label = "Spatially smooth clusters",
            icon = icon("play")
          ),

        ),
        mainPanel(
          plotOutput(ns('plotSmoothedClusters')),
        ))
    ),
    # conditionalPanel(
    #   id = ns('bayesSpaceDE'),
    #   ns=ns,
    #   condition = "input.clusteringApproach == 'BayesSpace'",
    sidebarLayout(
      # Sidebar panel for inputs ----
      sidebarPanel(
        selectInput(ns('clusters'), 'Clusters to use:', c('cluster','smoothed_cluster'),
                    selected = 'smoothed_cluster'),

        # Input: Selector clusters to compare
        selectInput(ns('variable1'), 'Cluster set 1:', choices = 1:2,selected = 1,multiple=T),
        selectInput(ns('variable2'), 'Cluster set 2:', choices = 1:2,selected = 2,multiple=T),
        selectInput(ns('pipeline'), 'DE pipeline:', c("t-test", "Wilcox rank sum"),selected = 'Wilcox rank sum'),

        #DE thresholds

        sliderInput(ns('pvalThreshold'), label = 'Adjusted p-value threshold',
                    min = 0, value = 0.05, max = 1, step = 0.005),
        sliderInput(ns('lfcThreshold'), label = 'log2 fold change threshold',
                    min = 0, value = 1, max = 5, step = 0.1),

        #Only start DE when button is pressed
        actionButton(ns('goDE'), label = 'Start differential intensity analysis'),

        #download file name and button
        textInput(ns('fileNameTable'),'File name for download', value ='DIAset.csv', placeholder = 'DIAset.csv'),
        downloadButton(ns('downloadTable'), 'Download Table'),
        hr(),
        tags$b("Peak selection"),
        div("\nSelect peaks of interest by clicking on the corresponds rows in the table\n"),
        div(style="margin-bottom:10px"),
        actionButton(ns('resetSelection'), label = "Reset row selection"),
        div(style="margin-bottom:10px"),
        actionButton(ns('selectTop50'), label = "Select top 50 peaks")
      ),

      #Main panel for displaying table of DE peaks
      mainPanel(
        DT::DTOutput(ns('data'))
      )
    ))
    #)

  }else{
    NULL
  }
}

#' @rdname RegionSpatialClusterPanel
#' @export
RegionSpatialClusterPanelServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, anno, gcd.values){

  moduleServer(id, function(input, output, session){

    cluster.names <- reactive({
      if (input[['clusteringApproach']]=='BayesSpace'){
        return(c('smoothed_cluster','cluster'))
      } else if (input[['clusteringApproach']]=='Lasso selection'){
        return('selected')
      }
    })

    cluster.options <- reactive({
      if (input[['clusteringApproach']]=='BayesSpace'){
        print(as.character(1:as.numeric(input[['numClusters']])))
        return(as.character(1:as.numeric(input[['numClusters']])))
      } else if (input[['clusteringApproach']]=='Lasso selection'){
        return(c('Selected','NotSelected'))
      }    })

    observe({
    updateSelectInput(session,'variable1', 'Cluster set 1:', choices = cluster.options(),selected = cluster.options()[1])
    updateSelectInput(session,'variable2', 'Cluster set 2:', choices = cluster.options(),selected = cluster.options()[2])
    updateSelectInput(session,inputId = 'clusters', choices = cluster.names(), selected = cluster.names()[1])
    })

    get_clusters <- reactive({

      intensity.sub <- full.intensity.matrix[full.metadata[,colnames(bulk.metadata)[1]] == input[['sampleToCluster']],]
      metadata.sub <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] == input[['sampleToCluster']],]
      metadata.sub$Sample = metadata.sub[,colnames(bulk.metadata)[1]]
      pos <- metadata.sub[,c('x','y')]
      rownames(pos)=metadata.sub$spot_id
      rownames(intensity.sub)=metadata.sub$spot_id
      metadata.sub$spot_id = paste0('spot_',metadata.sub$spot_id)
      rownames(intensity.sub)=paste0('spot_',rownames(intensity.sub))

      if (input[['clusteringApproach']]=='BayesSpace'){
        colData = metadata.sub[,2:3]
        colnames(colData)[1:2]=c('col','row')
        colData$row = colData$row - min(colData$row)
        gcd = gcd.values[[input[['sampleToCluster']]]]
        colData$row = round((colData$row/gcd)+1)
        colData$col = colData$col - min(colData$col)
        colData$col = round((colData$col/gcd)+1)
        colData$imagerow = colData$row
        colData$imagecol = colData$col
        anno = data.frame(anno)
        rownames(anno)=paste0('mz_',anno$m_z)
        colnames(intensity.sub)=paste0('mz_',colnames(intensity.sub))
        rowData = anno[colnames(intensity.sub),]

        sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=as(t(intensity.sub), "dgCMatrix")),
                                    rowData=rowData,
                                    colData=colData)
        print('Created singlecellexperiment object')
        S4Vectors::metadata(sce)$BayesSpace.data <- list()
        S4Vectors::metadata(sce)$BayesSpace.data$platform <- 'ST'
        S4Vectors::metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
        SingleCellExperiment::logcounts(sce)=log2(SingleCellExperiment::counts(sce)+1)
        print('Performed log transform')

        sce <- scater::runPCA(sce, subset_row=rownames(SingleCellExperiment::counts(sce)), ncomponents=30,
                              exprs_values='logcounts', BSPARAM=BiocSingular::ExactParam())
        print('Ran PCA')
        # use a select input with 'Run optimisation' as an option or using a specific number too
#        sce <- qTune(sce, qs=seq(2, 10), platform="ST")
        set.seed(23)
        sce <- BayesSpace::spatialCluster(sce, q=as.numeric(input[['numClusters']]), platform="ST", d=30,
                              nrep=input[['nRep']], burn.in=input[['burnIn']],
                              init.method="mclust", model="t", gamma=2)
        metadata.sub$cluster = factor(sce$spatial.cluster,levels=1:as.numeric(input[['numClusters']]))
        print('Performed BayesSpace clustering')
        return(metadata.sub)
      }
    })  %>% bindEvent(input[["run_clustering"]])

    selection_image <- reactive({
      intensity.sub <- full.intensity.matrix[full.metadata[,colnames(bulk.metadata)[1]] == input[['sampleToCluster']],]
      metadata.sub <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] == input[['sampleToCluster']],]
      metadata.sub$Sample = metadata.sub[,colnames(bulk.metadata)[1]]
      print(head(metadata.sub))
      return(ggplot(metadata.sub,aes(x=x,y=y,color=annotations,fill=annotations,key=spot_id))+geom_point(size=0.8)+theme_classic()+
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       aspect.ratio = 1))
    })

    selected_image <- reactive({
      intensity.sub <- full.intensity.matrix[full.metadata[,colnames(bulk.metadata)[1]] == input[['sampleToCluster']],]
      metadata.sub <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] == input[['sampleToCluster']],]
      metadata.sub$Sample = metadata.sub[,colnames(bulk.metadata)[1]]
      metadata.sub$selected = 'NotSelected'
      select_data <- event_data("plotly_selected")
      if (!is.null(select_data)) {
        metadata.sub[metadata.sub$spot_id %in% select_data$key, "selected"] <- 'Selected'
      }
      return(metadata.sub)
    })

    smoothed_clusters <- reactive({
      current.metadata = get_clusters()
      current.metadata$smoothed_cluster = smoothed.cluster(current.metadata,gcd.values = gcd.values)
      return(current.metadata)
    }) %>% bindEvent(input[["run_smoothing"]])

    cluster_plot <- reactive({
      current.metadata = get_clusters()
      my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(x = x, y = y, color = cluster, fill = cluster)) +
        ggplot2::geom_tile() +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       aspect.ratio = 1)

      return(my_plot)
    }) %>% bindEvent(input[["run_clustering"]])

    smoothed_cluster_plot <- reactive({
      current.metadata = smoothed_clusters()
      my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(x = x, y = y, color = smoothed_cluster, fill = smoothed_cluster)) +
        ggplot2::geom_tile() +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       aspect.ratio = 1)

      return(my_plot)
    }) %>% bindEvent(input[["run_smoothing"]])

    DEresults <- reactive({
      shinyjs::disable("goDE")
      if (input[['clusteringApproach']]=='BayesSpace'){
      current.metadata <- smoothed_clusters()
      } else if (input[['clusteringApproach']]=='Lasso selection'){
        current.metadata <- selected_image()
      }
      condition.indices <- current.metadata[,input[["clusters"]]] %in% c(input[['variable1']], input[['variable2']])
      condition = current.metadata[condition.indices,input[["clusters"]]]
      condition = ifelse(condition %in% input[['variable1']],'clusterSet1','clusterSet2')
      intensity.sub = full.intensity.matrix[full.metadata[,colnames(bulk.metadata)[1]] == input[['sampleToCluster']],]
      intensity.sub <- intensity.sub[condition.indices,]
      rownames(intensity.sub)=current.metadata[condition.indices,'spot_id']
      # Need to add error if any group has <2 samples
      DEtable <- DEanalysis(
        intensity.matrix = t(intensity.sub),
        condition = condition,
        var1 = 'clusterSet1',
        var2 = 'clusterSet2',
        test = input[["pipeline"]],
        anno = anno
      )

      DEtableSubset <- DEtable %>%
        dplyr::filter(.data$pvalAdj < input[["pvalThreshold"]] & abs(.data$lfc) > input[['lfcThreshold']]) |>
        dplyr::arrange(.data$pvalAdj)

      #the thresholds are returned here so that MA/volcano and table display
      #don't use new thresholds without the button being used
      shinyjs::enable("goDE")
      return(list('DEtable' = DEtable,
                  "DEtableSubset" = DEtableSubset,
                  'pvalThreshold' = input[["pvalThreshold"]],
                  'lfcThreshold' = input[['lfcThreshold']]))
    }) %>%
      bindCache(input[["condition"]],
                input[['variable1']], input[['variable2']], input[["pipeline"]],
                input[["pvalThreshold"]],input[['lfcThreshold']]) |>
      bindEvent(input[["goDE"]])


    #Define output table (only DE peaks)
    dataTable <- reactive({
      DEresults()$DEtableSubset %>%
        DT::datatable() %>%
        DT::formatSignif(columns = c('pval', 'pvalAdj','lfc','log2_intensity'), digits = 3)
    })

    output[['data']] <- DT::renderDataTable(dataTable())

    output[['plotClusters']] <- renderPlot({
      cluster_plot()
    })

    output[['lassoSelection']] <- plotly::renderPlotly(
      plotly::ggplotly(selection_image(),width=800,height=700) %>% layout(dragmode = "lasso")
    )

    output[['lassoSelected']] <- renderPlot({
      selected_image()},width=800,height=700)

    output[['plotSmoothedClusters']] <- renderPlot({
      smoothed_cluster_plot()
    })

    output[['downloadSpatial']] <- downloadHandler(
      filename = function() { input[['spatialPlotFileName']] },
      content = function(file) {
        ggsave(file, plot = cluster_plot(), width=input[['spatialPlotWidth']],height=input[['spatialPlotHeight']],units = 'in')
      }
    )

    output[['downloadSmoothSpatial']] <- downloadHandler(
      filename = function() { input[['spatialSmoothPlotFileName']] },
      content = function(file) {
        ggsave(file, plot = smoothed_cluster_plot(), width=input[['spatialSmoothPlotWidth']],height=input[['spatialSmoothPlotHeight']],units = 'in')
      }
    )

    #DE data download
    output[['downloadTable']] <- downloadHandler(
      filename = function() {
        paste(input[['fileNameTable']])
      },
      content = function(file) {
        utils::write.csv(x = DEresults()$DEtableSubset, file = file, row.names = FALSE)
      }
    )

    return(reactive(return_object()))

  })
}


