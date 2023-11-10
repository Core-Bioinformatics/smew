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
          
          
          
        ),
        
        #Main panel for displaying table of enriched pathways
        mainPanel(
          plotOutput(ns('plotClusters')),
          plotOutput(ns('plotClusterProps')),
          plotOutput(ns('plotClusterPropsPerSample')))
        )
      )
  }else{
    NULL
  }
}

#' @rdname DEsummaryPanel
#' @export
RegionClusterPanelServer <- function(id, full.expression.matrix, full.metadata, bulk.metadata, anno){
  
  moduleServer(id, function(input, output, session){
    updateSelectizeInput(session, "peakToThreshold", choices = anno$display_name, server = TRUE, selected = anno$display_name[1])
    
    get_clusters <- reactive({
      current.expression.matrix <- full.expression.matrix[,full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']]]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']],]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      if (input[['clusteringApproach']]=='Intensity thresholding'){
      my_peak = anno[anno$display_name==input[['peakToThreshold']],]
      my_peak_expression <- t(current.expression.matrix[my_peak$m_z,])
      quantiles = quantile(my_peak_expression, prob=c(input[['thresholdLow']]/100,(100-input[['thresholdHigh']])/100), type=1)
      current.metadata$cluster = factor(ifelse(my_peak_expression<=quantiles[1],'Low',
                                        ifelse(my_peak_expression>=quantiles[2],'High','Medium')),levels=c('Low','Medium','High'))
      }
      if (input[['clusteringApproach']]=='k-means'){
        
 #     set.seed(23)
      current.expression.matrix <- scale(x = current.expression.matrix,center = T,scale = T)
      current.expression.matrix <- current.expression.matrix[complete.cases(current.expression.matrix),]
      clusters <- kmeans(t(current.expression.matrix),
                         centers = input[['numClusters']],
                         nstart = 1)
      current.metadata$cluster = factor(clusters$cluster,levels=1:input[['numClusters']])
      }
      return(current.metadata)
    })  %>% bindEvent(input[["run_clustering"]])
    
    cluster_plot <- reactive({
      current.metadata = get_clusters()
      my_plot <- ggplot(current.metadata,aes(x = x, y = y, color = cluster, fill = cluster)) +
                    geom_tile() +
                    facet_wrap(~current.metadata$Sample, nrow = floor(sqrt(length(input[['samplesToCluster']]))), scales = 'free') +
                    theme_classic() +
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank(),
                          axis.line.x = element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          axis.line.y = element_blank())
      
      return(my_plot)
    }) %>% bindEvent(input[["run_clustering"]])
    
    cluster_props <- reactive({
      current.metadata = get_clusters()
      current.metadata$SelectedMetadata = current.metadata[,input[['groupingMetadataBarPlot']]]
      ggplot(current.metadata,aes(y=SelectedMetadata,fill=cluster)) +
        geom_bar(position = 'fill') +
        ylab(input[['groupingMetadataBarPlot']]) +
        xlab('Proportion of spots')+
        theme_classic()
    })
    
    cluster_props_persample <- reactive({
      current.metadata = get_clusters()
      current.metadata$SelectedMetadata = current.metadata[,input[['groupingMetadataBox']]]
      current.metadata.count = current.metadata %>%
        group_by(Sample, cluster) %>%
        summarise(n = n()) %>%
        mutate(freq = n / sum(n))
      current.metadata.count = merge(data.frame(current.metadata.count),unique(current.metadata[,c('Sample','SelectedMetadata')]))
      ggplot(current.metadata.count,aes(fill=SelectedMetadata,y=freq,x=cluster)) +
        geom_boxplot() +
        scale_fill_discrete(name=input[['groupingMetadataBox']]) +
        ylab('Proportion of spots per sample') +
        theme_classic()
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
    })
    
    output[['plotClusterProps']] <- renderPlot({
      cluster_props()
    })
    
    output[['plotClusterPropsPerSample']] <- renderPlot({
      cluster_props_persample()
    })
    
    return(reactive(return_object()))
           
  })
}


