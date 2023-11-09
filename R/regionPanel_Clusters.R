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
          
          numericInput(inputId = ns("numClusters"),
                       label = "Select number of clusters",
                       min = 2,
                       max = 10,
                       value = 3,
                       step = 1),
          
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
    
    get_clusters <- reactive({
      print(full.expression.matrix[1:5,1:5])
      print(summary(full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']]))
      print(length(full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']]))
      print(head(full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']]))
      print(ncol(full.expression.matrix))
      current.expression.matrix <- full.expression.matrix[,full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']]]
      print(dim(current.expression.matrix))
      current.expression.matrix <- scale(x = current.expression.matrix,center = T,scale = T)
      print(dim(current.expression.matrix))
      print(current.expression.matrix[1:5,1:5])
      print(current.expression.matrix[(nrow(current.expression.matrix)-5):(nrow(current.expression.matrix)),(ncol(current.expression.matrix)-5):(ncol(current.expression.matrix))])
      print(colnames(current.expression.matrix)[colSums(is.na(current.expression.matrix)) > 0])
      print(head(current.expression.matrix[224,]))
      current.expression.matrix <- current.expression.matrix[complete.cases(current.expression.matrix),]
      print(which(rowSums(is.na(current.expression.matrix)) > 0))
 #     set.seed(23)
      clusters <- kmeans(t(current.expression.matrix),
                         centers = input[['numClusters']],
                         nstart = 1)
      print(str(clusters))
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToCluster']],]
      print(head(current.metadata))
      print(input[['groupingMetadataBarPlot']])
      print(head(current.metadata))
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      current.metadata$cluster = factor(clusters$cluster,levels=1:input[['numClusters']])
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
      print(head(current.metadata.count))
      ggplot(current.metadata.count,aes(fill=SelectedMetadata,y=freq,x=cluster)) +
        geom_boxplot() +
        scale_fill_discrete(name=input[['groupingMetadataBox']]) +
        ylab('Proportion of spots per sample') +
        theme_classic()
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
    
    return(reactive(get_clusters()[,c('spot_id','cluster')]))
  })
}


