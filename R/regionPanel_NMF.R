#' @rdname RegionNMFPanel
#' @export
RegionNMFPanelUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- NS(id)
  # add option to name cluster and retain it
  if(show){
    tabPanel(
      'NMF',
      sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(
          # samples to include
          # number of clusters
          # any other parameters
          selectInput(
            inputId = ns("samplesToFactor"),
            label = "Select samples for factorising",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[,1]),
            multiple = TRUE
          ),
          numericInput(inputId = ns("numFactors"),
                       label = 'Number of factors',
                       value = 10,
                       min = 2,
                       max = 50,
                       step = 1),
          numericInput(inputId = ns("seed"),
                       label = 'Set seed',
                       value = 23,
                       min = 1,
                       max = 1000,
                       step = 1),
          # button to start clustering
          actionButton(
            inputId = ns("run_nmf"),
            label = "Find NMF factors",
            icon = icon("play")
          ),
          selectInput(inputId=ns("focus_NMF"),
                      label = 'Select NMF factor',
                      choices = 1:10,selected = 1),
          selectInput(inputId = ns("groupingMetadataBox"),
                      label = "Metadata to group boxplot on",
                      choices = colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y'))],
                      selected = colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y'))][1],
                      multiple = FALSE),
          selectInput(inputId = ns("colourMetadataBox"),
                      label = "Metadata to colour boxplots by",
                      choices = colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y'))],
                      selected = colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y'))][1],
                      multiple = FALSE)



        ),

        #Main panel for displaying table of enriched pathways
        mainPanel(
          plotOutput(ns('plotNMF')),
          plotOutput(ns('plotPerSampleNMF')),
          plotOutput(ns('plotNMFFeatureWeights')))
          # plotOutput(ns('plotClusterProps')),
          # plotOutput(ns('plotClusterPropsPerSample')))
      )
    )
  }else{
    NULL
  }
}

#' @rdname RegionNMFPanel
#' @export
RegionNMFPanelServer <- function(id, full.expression.matrix, full.metadata, bulk.metadata, anno){

  moduleServer(id, function(input, output, session){
    observe(
      updateSelectizeInput(session, "focus_NMF", choices = 1:input[['numFactors']], server = TRUE, selected = 1)
    )


    get_nmf <- reactive({
      current.expression.matrix <- full.expression.matrix[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToFactor']],]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToFactor']],]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      # current.expression.matrix <- scale(x = current.expression.matrix,center = T,scale = T)
      # current.expression.matrix <- current.expression.matrix[complete.cases(current.expression.matrix),]
      nmf.factors <- RcppML::nmf(as.matrix(current.expression.matrix),seed = input[['seed']],k = input[['numFactors']])
      nmf.factor.weights = data.frame(nmf.factors$w)
      nmf.factor.features = data.frame(nmf.factors$h)
      colnames(nmf.factor.weights)=paste0('NMF_',gsub('X','',colnames(nmf.factor.weights)))
      current.metadata = cbind(current.metadata,nmf.factor.weights)
      colnames(nmf.factor.features)=colnames(current.expression.matrix)
      rownames(nmf.factor.features)=paste0('NMF_',gsub('X','',rownames(nmf.factor.features)))
      return.list = list('metadata'=current.metadata,'feature_weights'=nmf.factor.features)
      return(return.list)
    })  %>% bindEvent(input[["run_nmf"]])

    nmf_plot <- reactive({
      current.metadata = get_nmf()$metadata
      my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(x = x, y = y, color = get(paste0('NMF_',input[['focus_NMF']])), fill = get(paste0('NMF_',input[['focus_NMF']])))) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient(low='white',high='darkred',name=paste0('NMF_',input[['focus_NMF']]))+
        ggplot2::scale_color_gradient(low='white',high='darkred',name=paste0('NMF_',input[['focus_NMF']]))+
        ggplot2::facet_wrap(~current.metadata$Sample, nrow = floor(sqrt(length(input[['samplesToFactor']]))), scales = 'free') +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank())

      return(my_plot)
    }) %>% bindEvent(input[["run_nmf"]])

    nmf_persample <- reactive({
      current.metadata = get_nmf()$metadata
      my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(y = get(paste0('NMF_',input[['focus_NMF']])), x = get(input[['groupingMetadataBox']]), fill = get(input[['colourMetadataBox']])))+geom_boxplot()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab(paste0('NMF_',input[['focus_NMF']]))+
        scale_fill_discrete(name=input[['colourMetadataBox']])+
        scale_fill_discrete(name=input[['groupingMetadataBox']])+
        xlab(input[['groupingMetadataBox']])
      return(my_plot)
    })

    nmf_featureweights <- reactive({
      nmf.weights = data.frame(t(get_nmf()$feature_weights))
      nmf.weights = nmf.weights[order(-nmf.weights[,paste0('NMF_',input[['focus_NMF']])]),]
      nmf.weights = head(nmf.weights,20)
      nmf.weights$peak = rownames(nmf.weights)
      nmf.weights$peak = factor(nmf.weights$peak,levels=rev(nmf.weights$peak))
      ggplot(nmf.weights,aes(y=peak,x=get(paste0('NMF_',input[['focus_NMF']]))))+geom_bar(stat='identity')+theme_classic()+xlab(paste0('NMF_',input[['focus_NMF']]))
    })
    # cluster_props <- reactive({
    #   current.metadata = get_clusters()
    #   current.metadata$SelectedMetadata = current.metadata[,input[['groupingMetadataBarPlot']]]
    #   ggplot2::ggplot(current.metadata,ggplot2::aes(y=SelectedMetadata,fill=cluster)) +
    #     ggplot2::geom_bar(position = 'fill') +
    #     ggplot2::ylab(input[['groupingMetadataBarPlot']]) +
    #     ggplot2::xlab('Proportion of spots')+
    #     ggplot2::theme_classic()
    # })
    #
    # cluster_props_persample <- reactive({
    #   current.metadata = get_clusters()
    #   current.metadata$SelectedMetadata = current.metadata[,input[['groupingMetadataBox']]]
    #   current.metadata.count = current.metadata |>
    #     dplyr::group_by(Sample, cluster) |>
    #     dplyr::summarise(n = dplyr::n()) |>
    #     dplyr::mutate(freq = n / sum(n))
    #   current.metadata.count = merge(data.frame(current.metadata.count),unique(current.metadata[,c('Sample','SelectedMetadata')]))
    #   ggplot2::ggplot(current.metadata.count,ggplot2::aes(fill=SelectedMetadata,y=freq,x=cluster)) +
    #     ggplot2::geom_boxplot() +
    #     ggplot2::scale_fill_discrete(name=input[['groupingMetadataBox']]) +
    #     ggplot2::ylab('Proportion of spots per sample') +
    #     ggplot2::theme_classic()
    # })
    #
    # return_object <- reactive({
    #   rownames(full.metadata)<-full.metadata$spot_id
    #   merged.metadata = merge(full.metadata,get_clusters(),all.x=T,sort=F)
    #   rownames(merged.metadata)=merged.metadata$spot_id
    #   merged.metadata = merged.metadata[rownames(full.metadata),]
    #   merged.metadata$cluster = as.character(merged.metadata$cluster)
    #   merged.metadata <- tidyr::replace_na(merged.metadata, list(cluster = 'None'))
    #   return(merged.metadata)
    # })
    #
    output[['plotNMF']] <- renderPlot({
      nmf_plot()
    })

    output[['plotPerSampleNMF']] <- renderPlot({
      nmf_persample()
    })

    output[['plotNMFFeatureWeights']] <- renderPlot({
      nmf_featureweights()
    })

    # output[['plotClusterProps']] <- renderPlot({
    #   cluster_props()
    # })
    #
    # output[['plotClusterPropsPerSample']] <- renderPlot({
    #   cluster_props_persample()
    # })
    #
    # return(reactive(return_object()))

  })
}


