#' @rdname RegionClusterPanel
#' @export
RegionClusterPanelUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- NS(id)
  # add option to name cluster and retain it
  if(show){
    tabPanel(
      'Clustering',
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
                tags$li("Regions/clusters can be identified using k-means clustering or intensity threshold-based segmentation (more methods to be added soon)."),
                tags$li("Select which samples to use for clusters - note that the intensity thresholding method will only calculate a threshold based on the samples you include."),
                tags$li("After clustering, the spatial distribution of the clusters is shown, followed by a barplot showing the proportion of pixels in your selected metadata column in each cluster and a boxplot showing the distribution of proportion of pixels per sample in each cluster grouped by your selected metadata column."),
                tags$li("The most recent clusters identified are passed onto the region-level differential intensity analysis tab"),
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
          selectInput(
            inputId = ns("samplesToCluster"),
            label = "Select samples for clustering",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[,1]),
            multiple = TRUE
          ),
          selectInput(inputId = ns("clusteringApproach"),
                      "Clustering method",
                      choices = c("k-means","Intensity thresholding")),
          conditionalPanel(
            id = ns('thresholdingOptions'),
            ns=ns,
            condition = "input.clusteringApproach == 'Intensity thresholding'",
            selectInput(ns('peakToThreshold'), 'Peak to use for intensity thresholding', , multiple = FALSE, choices = character(0)),
            sliderInput(ns('thresholdHigh'),label = 'Percentage of data to categorise as high intensity',value = 25,min = 1,max = 49,step = 1),
            sliderInput(ns('thresholdLow'),label = 'Percentage of data to categorise as low intensity',value = 25,min = 1,max = 49,step = 1)
          ),
          conditionalPanel(
            id = ns('kmeansOptions'),
            ns=ns,
            condition = "input.clusteringApproach == 'k-means'",
            numericInput(inputId = ns("numClusters"),
                         label = "Select number of clusters",
                         min = 2,
                         max = 10,
                         value = 3,
                         step = 1)),

          # button to start clustering
          actionButton(
            inputId = ns("run_clustering"),
            label = "Find clusters",
            icon = icon("play")
          ),

          selectInput(inputId = ns("groupingMetadataBarPlot"),
                      label = "Metadata to group barplot on",
                      choices = colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y'))],
                      selected = colnames(bulk.metadata)[1],
                      multiple = FALSE),
          selectInput(inputId = ns("groupingMetadataBox"),
                      label = "Metadata to per-sample cluster proportion boxplot on",
                      choices = colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata)[1]))],
                      selected = colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata)[1]))][1],
                      multiple = FALSE),
          div(style = "margin-top:10px"),
          dropMenu(
            circleButton(ns("downloads"), icon = icon("download"),status = "success"),
            tags$div(
              tags$h3("Downloads"),
              fluidRow(
                       column(5,offset=0,
                              tags$h4("Barplot"),
                              textInput(ns('spatialPlotFileName'),'File name for download', value ='spatialClusters.png', placeholder = 'spatialClusters.png'),
                              numericInput(ns('spatialPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                              numericInput(ns('spatialPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                              downloadButton(ns('downloadSpatial'), 'Download clusters')),
                       column(5,offset=0,
                              tags$h4("Upset"),
                              textInput(ns('barPlotFileName'),'File name for download', value ='bar.png', placeholder = 'bar.png'),
                              numericInput(ns('barPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                              numericInput(ns('barPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                              downloadButton(ns('downloadBarPlot'), 'Download barplot'))),
                fluidRow(
                        column(10,offset=0,
                               textInput(ns('boxPlotFileName'),'File name for download', value ='box.png', placeholder = 'box.png'),
                               numericInput(ns('boxPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                               numericInput(ns('boxPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                               downloadButton(ns('downloadBoxPlot'), 'Download boxplot'))),
              theme = "light-border",
            placement = "right",
            arrow = FALSE
          )),



        ),

        #Main panel for displaying table of enriched pathways
        mainPanel(
          fluidRow(column=10,plotOutput(ns('plotClusters'),height = 600)),
          plotOutput(ns('plotClusterProps')),
          plotOutput(ns('plotClusterPropsPerSample')))
        )
      )
  }else{
    NULL
  }
}

#' @rdname RegionClusterPanel
#' @export
RegionClusterPanelServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, anno){

  moduleServer(id, function(input, output, session){
    updateSelectizeInput(session, "peakToThreshold", choices = anno$display_name, server = TRUE, selected = anno$display_name[1])

    get_clusters <- reactive({
      current.intensity.matrix <- full.intensity.matrix[,full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']]]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']],]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      if (input[['clusteringApproach']]=='Intensity thresholding'){
      my_peak = anno[anno$display_name==input[['peakToThreshold']],]
      my_peak_expression <- t(current.intensity.matrix[my_peak$m_z,])
      quantiles = quantile(my_peak_expression, prob=c(input[['thresholdLow']]/100,(100-input[['thresholdHigh']])/100), type=1)
      current.metadata$cluster = factor(ifelse(my_peak_expression<=quantiles[1],'Low',
                                        ifelse(my_peak_expression>=quantiles[2],'High','Medium')),levels=c('Low','Medium','High'))
      }
      if (input[['clusteringApproach']]=='k-means'){

 #     set.seed(23)
      current.intensity.matrix <- scale(x = current.intensity.matrix,center = T,scale = T)
      current.intensity.matrix <- current.intensity.matrix[complete.cases(current.intensity.matrix),]
      clusters <- kmeans(t(current.intensity.matrix),
                         centers = input[['numClusters']],
                         nstart = 1)
      current.metadata$cluster = factor(clusters$cluster,levels=1:input[['numClusters']])
      }
      return(current.metadata)
    })  %>% bindEvent(input[["run_clustering"]])

    cluster_plot <- reactive({
      current.metadata = get_clusters()
      my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(x = x, y = y, color = cluster, fill = cluster)) +
        ggplot2::geom_tile() +
        ggplot2::facet_wrap(~current.metadata$Sample, nrow = floor(sqrt(length(input[['samplesToCluster']]))), scales = 'free') +
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

    cluster_props <- reactive({
      current.metadata = get_clusters()
      current.metadata$SelectedMetadata = current.metadata[,input[['groupingMetadataBarPlot']]]
      ggplot2::ggplot(current.metadata,ggplot2::aes(y=SelectedMetadata,fill=cluster)) +
        ggplot2::geom_bar(position = 'fill') +
        ggplot2::ylab(input[['groupingMetadataBarPlot']]) +
        ggplot2::xlab('Proportion of spots')+
        ggplot2::theme_classic()
    })

    cluster_props_persample <- reactive({
      current.metadata = get_clusters()
      current.metadata$SelectedMetadata = current.metadata[,input[['groupingMetadataBox']]]
      current.metadata.count = current.metadata |>
        dplyr::group_by(Sample, cluster) |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::mutate(freq = n / sum(n))
      current.metadata.count = merge(data.frame(current.metadata.count),unique(current.metadata[,c('Sample','SelectedMetadata')]))
      ggplot2::ggplot(current.metadata.count,ggplot2::aes(fill=SelectedMetadata,y=freq,x=cluster)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_fill_discrete(name=input[['groupingMetadataBox']]) +
        ggplot2::ylab('Proportion of spots per sample') +
        ggplot2::theme_classic()
    })

    return_object <- reactive({
      rownames(full.metadata)<-full.metadata$spot_id
      merged.metadata = merge(full.metadata,get_clusters(),all.x=T,sort=F)
      rownames(merged.metadata)=merged.metadata$spot_id
      merged.metadata = merged.metadata[rownames(full.metadata),]
      merged.metadata$cluster = as.character(merged.metadata$cluster)
      merged.metadata <- tidyr::replace_na(merged.metadata, list(cluster = 'None'))
      return(merged.metadata)
    })

    output[['plotClusters']] <- renderPlot({
      cluster_plot()
    },height=600)

    output[['plotClusterProps']] <- renderPlot({
      cluster_props()
    })

    output[['plotClusterPropsPerSample']] <- renderPlot({
      cluster_props_persample()
    })
    output[['downloadSpatial']] <- downloadHandler(
      filename = function() { input[['spatialPlotFileName']] },
      content = function(file) {
        ggsave(file, plot = cluster_plot(), width=input[['spatialPlotWidth']],height=input[['spatialPlotHeight']],units = 'in')
      }
    )

    output[['downloadBarPlot']] <- downloadHandler(
      filename = function() { input[['barPlotFileName']] },
      content = function(file) {
        ggsave(file, plot = cluster_props(), width=input[['barPlotWidth']],height=input[['barPlotHeight']],units = 'in')
      }
    )

    output[['downloadBoxPlot']] <- downloadHandler(
      filename = function() { input[['boxPlotFileName']] },
      content = function(file) {
        ggsave(file, plot = cluster_props_persample(), width=input[['boxPlotWidth']],height=input[['boxPlotHeight']],units = 'in')
      }
    )

    return(reactive(return_object()))

  })
}


