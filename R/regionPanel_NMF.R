ColorBlender <- function (
    data,
    channels.use = NULL
) {
  rgb.order <- setNames(1:3, c("red", "green", "blue"))
  if (!length(channels.use) == ncol(data)) {
    stop(paste0("channels.use must be same length as number of features or dimensions"))
  } else if (!all(channels.use %in% names(rgb.order))) {
    stop("Invalid color names in channels.use. Valid options are: 'red', 'green' and 'blue'")
  } else if (sum(duplicated(channels.use))){
    stop("Duplicate color names are not allowed in channels.use")
  }
  col.order <- rgb.order[channels.use]

  if (ncol(data) == 2) {
    first_vec <- data[, 1]
    second_vec <- data[, 2]
    data <- matrix(data = 0, nrow = nrow(data), ncol = 3)
    data[, col.order[1]] <- first_vec; data[, col.order[2]] <- second_vec
  } else if (ncol(data) == 3) {
    data <- data[, col.order]
  }
  color.codes <- rgb(data)
}


#' @rdname RegionDimRedPanel
#' @export
RegionDimRedPanelUI <- function(id, bulk.metadata, full.metadata, full.intensity.matrix,show = TRUE){
  ns <- NS(id)
  # add option to name cluster and retain it
  if(show){
    tabPanel(
      'Dimensionality reduction',
      sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(
          # samples to include
          # number of clusters
          # any other parameters
          dropMenu(
            circleButton(ns("info_dim_reduction"), icon = icon("info"),status = "success"),
            tags$div(
              tags$h3("Dimensionality reduction"),
              tags$ul(
                tags$li("This tab will perform dimensionality reduction on all pixels across whichever samples are selected."),
                tags$li("NMF and PCA are currently available."),
                tags$li("The number of dimensions (i.e. components/factors) computed is user-selected and a random seed can also be selected."),
                tags$li("The dimensionality reduction is only performed once 'Run dimensionality reduction' button is pressed."),
                tags$li("The distribution of factor/component values by metadata column can be visualised and the top contributing peaks extracted."),
                tags$li("The resulting components/factors can be visualised spatially and a UMAP can further be calculated to visualise the dimensionality reduction in 2 dimensions."),
                tags$li("Pixels in the UMAP can be coloured according to metadata information or peak intensity."),
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),

          selectInput(inputId = ns("dimReduction"),label = "Dimensionality reduction",choices = c('NMF','PCA'),selected = 'NMF'),
          selectInput(
            inputId = ns("samplesToFactor"),
            label = "Select samples for infer dimensionality reduction",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[,1]),
            multiple = TRUE
          ),
          numericInput(inputId = ns("numDimensions"),
                       label = 'Number of dimensions to infer',
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
            inputId = ns("run_dimred"),
            label = "Run dimensionality reduction",
            icon = icon("play")
          ),
          selectInput(inputId=ns("focus_dimension"),
                      label = 'Select dimension of interest',
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
                      multiple = FALSE),
          selectInput(inputId = ns("colourUMAP"),
                      label = "Metadata to colour UMAP",
                      choices = c(colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y'))],colnames(full.intensity.matrix)),
                      selected = colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y'))][1],
                      multiple = FALSE)



        ),

        #Main panel for displaying table of enriched pathways
        mainPanel(
          plotOutput(ns('plotDimRed')),
          plotOutput(ns('plotPerSampleDimRed')),
          plotOutput(ns('plotDimRedFeatureWeights')),
          plotOutput(ns('DimRedUMAPSpatial')),
          plotOutput(ns('DimRedUMAP')))
        # plotOutput(ns('plotClusterProps')),
        # plotOutput(ns('plotClusterPropsPerSample')))
      )
    )
  }else{
    NULL
  }
}

#' @rdname RegionDimRedPanel
#' @export
RegionDimRedPanelServer <- function(id, full.intensity.matrix, full.metadata, bulk.metadata, anno){

  moduleServer(id, function(input, output, session){
    observe(
      updateSelectizeInput(session, "focus_DimRed", choices = 1:input[['numDimensions']], server = TRUE, selected = 1)
    )
    get_subset_exp <- reactive({
      current.intensity.matrix <- full.intensity.matrix[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToFactor']],]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToFactor']],]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      return(list('exp'=current.intensity.matrix,'meta'=current.metadata))
    }) %>% bindEvent(input[["run_dimred"]])

    get_nmf <- reactive({
      current.intensity.matrix <- get_subset_exp()$exp
      current.metadata <- get_subset_exp()$meta      # current.intensity.matrix <- scale(x = current.intensity.matrix,center = T,scale = T)
      # current.intensity.matrix <- current.intensity.matrix[complete.cases(current.intensity.matrix),]
      nmf.factors <- RcppML::nmf(as.matrix(current.intensity.matrix),seed = input[['seed']],k = input[['numDimensions']])
      nmf.factor.weights = data.frame(nmf.factors$w)
      nmf.factor.features = data.frame(nmf.factors$h)
      colnames(nmf.factor.weights)=paste0('NMF_',gsub('X','',colnames(nmf.factor.weights)))
      current.metadata = cbind(current.metadata,nmf.factor.weights)
      colnames(nmf.factor.features)=colnames(current.intensity.matrix)
      rownames(nmf.factor.features)=paste0('NMF_',gsub('X','',rownames(nmf.factor.features)))
      return.list = list('metadata'=current.metadata,'feature_weights'=nmf.factor.features,'dimRed'=input[['dimReduction']],'numDimensions'=input[['numDimensions']])
      return(return.list)
    })  %>% bindEvent(input[["run_dimred"]])

    get_pca <- reactive({
      current.intensity.matrix <- get_subset_exp()$exp
      current.metadata <- get_subset_exp()$meta
      current.intensity.matrix <- current.intensity.matrix[,apply(current.intensity.matrix, 2, var, na.rm=TRUE) != 0]
      # current.intensity.matrix <- scale(x = current.intensity.matrix,center = T,scale = T)
      # current.intensity.matrix <- current.intensity.matrix[complete.cases(current.intensity.matrix),]
      print('I am running')
      pc.components <- prcomp(as.matrix(current.intensity.matrix),center = TRUE,scale.=TRUE,rank. = input[["numDimensions"]])
      print('I have calculate PCs')
      pc.component.weights = data.frame(pc.components$x)
      pc.component.features = data.frame(t(pc.components$rotation))
      colnames(pc.component.weights)=paste0('PCA_',1:input[["numDimensions"]])
      current.metadata = cbind(current.metadata,pc.component.weights)
      colnames(pc.component.features)=colnames(current.intensity.matrix)
      rownames(pc.component.features)=paste0('PCA_',1:input[["numDimensions"]])
      print('I am here')
      return.list = list('metadata'=current.metadata,'feature_weights'=pc.component.features,'dimRed'=input[['dimReduction']],'numDimensions'=input[['numDimensions']])
      return(return.list)
    })  %>% bindEvent(input[["run_dimred"]])

    get_dimred <- reactive({
      if (input[["dimReduction"]]=="NMF"){
        return(get_nmf())
      } else {
        return(get_pca())
      }
    })
    nmf_plot <- reactive({
      current.metadata = get_dimred()$metadata
      my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(x = x, y = y, color = get(paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])), fill = get(paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])))) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low='darkblue',high='darkred',mid='white',name=paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))+
        ggplot2::scale_color_gradient2(low='darkblue',high='darkred',mid='white',name=paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))+
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
    }) %>% bindEvent(input[["run_dimred"]])

    nmf_persample <- reactive({
      current.metadata = get_dimred()$metadata
      my_plot <- ggplot2::ggplot(current.metadata,ggplot2::aes(y = get(paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])), x = get(input[['groupingMetadataBox']]), fill = get(input[['colourMetadataBox']])))+geom_boxplot()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab(paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))+
        scale_fill_discrete(name=input[['colourMetadataBox']])+
        scale_fill_discrete(name=input[['groupingMetadataBox']])+
        xlab(input[['groupingMetadataBox']])
      return(my_plot)
    })

    nmf_featureweights <- reactive({
      nmf.weights = data.frame(t(get_dimred()$feature_weights))
      nmf.weights = nmf.weights[order(-nmf.weights[,paste0(get_dimred()$dimRed,'_',input[['focus_dimension']])]),]
      nmf.weights = head(nmf.weights,20)
      nmf.weights$peak = rownames(nmf.weights)
      nmf.weights$peak = factor(nmf.weights$peak,levels=rev(nmf.weights$peak))
      ggplot(nmf.weights,aes(y=peak,x=get(paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))))+geom_bar(stat='identity')+theme_classic()+xlab(paste0(get_dimred()$dimRed,'_',input[['focus_dimension']]))
    })

    nmf_umap_prep <- reactive({
      current.metadata = get_dimred()$metadata
      dim.red.weights = current.metadata[,paste0(get_dimred()$dimRed,'_',1:get_dimred()$numDimensions)]
      my.umap = uwot::umap(dim.red.weights,n_neighbors = 30,n_components = 2)
      my.umap = as.data.frame(my.umap)
      colnames(my.umap)=c('UMAP_1','UMAP_2')
      current.metadata = cbind(current.metadata,my.umap)
      current.metadata$UMAP_1 = scales::rescale(current.metadata$UMAP_1)
      current.metadata$UMAP_2 = scales::rescale(current.metadata$UMAP_2)
      colors=ColorBlender(current.metadata[,c('UMAP_1','UMAP_2')],channels.use = c("red","blue"))
      current.metadata$my.color = colors
      return(current.metadata)})

    nmf_umap <- reactive({
      current.metadata = nmf_umap_prep()
      if (input[['colourUMAP']]%in%colnames(get_subset_exp()$exp)){
        current.metadata$selectedPeak = get_subset_exp()$exp[,input[['colourUMAP']]]
        current.metadata$selectedPeak = pmin(quantile(current.metadata$selectedPeak,0.95),current.metadata$selectedPeak)
        current.metadata$selectedPeak = pmax(quantile(current.metadata$selectedPeak,0.05),current.metadata$selectedPeak)
        current.metadata = current.metadata[order(current.metadata$selectedPeak),]
        colour.function = ggplot2::scale_fill_gradient(name=input[['colourUMAP']],low = "lightgrey", high = "brown")

      } else {
        current.metadata$selectedPeak = current.metadata[,input[['colourUMAP']]]
        colour.function = ggplot2::scale_color_discrete(name=input[['colourUMAP']])
      }
      return(list('SpatialView'=ggplot(current.metadata,aes(x=x,y=y))+
                    geom_tile(color=current.metadata$my.color,fill=current.metadata$my.color)+
                    theme_classic()+facet_wrap(~current.metadata$Group,scales='free')+
                    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                                   axis.text.x=ggplot2::element_blank(),
                                   axis.ticks.x=ggplot2::element_blank(),
                                   axis.line.x = ggplot2::element_blank(),
                                   axis.title.y=ggplot2::element_blank(),
                                   axis.text.y=ggplot2::element_blank(),
                                   axis.ticks.y=ggplot2::element_blank(),
                                   axis.line.y = ggplot2::element_blank()),
                  'UMAP'=ggplot(current.metadata,aes(x=UMAP_1,y=UMAP_2,color=selectedPeak))+
                    geom_point()+
                    theme_classic()+colour.function))

    })

    output[['plotDimRed']] <- renderPlot({
      nmf_plot()
    })

    output[['plotPerSampleDimRed']] <- renderPlot({
      nmf_persample()
    })

    output[['plotDimRedFeatureWeights']] <- renderPlot({
      nmf_featureweights()
    })

    output[['DimRedUMAPSpatial']] <- renderPlot({
      nmf_umap()$SpatialView
    })

    output[['DimRedUMAP']] <- renderPlot({
      nmf_umap()$UMAP
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


