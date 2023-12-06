#' @rdname IntroSpatialVisPanel
#' @export
PixelSVMPanelUI <- function(id, bulk.metadata, full.metadata, full.expression.matrix, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Spatial visualisation',
      selectInput(
        inputId = ns("samplesToInclude"),
        label = "Select samples to include",
        choices = unique(bulk.metadata[,1]),
        selected = unique(bulk.metadata[,1])[1],
        multiple = TRUE,width = '100%'
      ),
      actionButton(ns("run_SVM"),'Identify spatially variable metabolites'),
      sidebarLayout(
        sidebarPanel(
          numericInput(ns("topNumber"),label = 'Pick top n spatially variable metabolites',value = 20,min = 1,max = 100,step = 1),
          conditionalPanel(id=ns('upsetPlotPanel'),
                           ns = ns,
                           condition = "output.check > 1",
                           selectInput(ns("metadataBarplot"), "Metadata to split barplot:", multiple = FALSE, choices = colnames(bulk.metadata),selected=colnames(bulk.metadata)[length(colnames(bulk.metadata))])),
          selectInput(ns("peakName"), "Peaks to display:", multiple = FALSE, choices = c()),
        ),
        mainPanel(
          conditionalPanel(id=ns('upsetPlotPanel'),
            ns = ns,
            condition = "output.check > 1",
            plotOutput(ns('SVMUpset')),
            plotOutput(ns('SVMBarPlot'))),
          conditionalPanel(id=ns('tablePanel'),
            ns = ns,
            condition = "output.check == 1",
            DT::DTOutput(ns('SVMTable'))))),
          plotOutput(ns('plotPeak')),
    )
  }else{
    NULL
  }
}

#' @rdname IntroSpatialVisPanel
#' @export
PixelSVMPanelServer <- function(id, bulk.metadata, full.metadata, full.expression.matrix, anno){

  moduleServer(id, function(input, output, session){

    get_subset_exp <- reactive({
      current.expression.matrix <- full.expression.matrix[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToInclude']],]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToInclude']],]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      return(list('exp'=current.expression.matrix,'meta'=current.metadata))
    }) %>% bindEvent(input[["run_SVM"]])

    run_svm <- reactive({
      current.expression.matrix <- get_subset_exp()$exp
      current.metadata <- get_subset_exp()$meta
      coords = current.metadata[,c('spot_id','x','y',colnames(bulk.metadata)[1])]
      colnames(coords)=c('barcode','x','y','sampleID')
      coords$barcode = as.character(coords$barcode)
      coords = as.tibble(coords)
      coords$sampleID = as.numeric(as.factor(coords$sampleID))
      spatnet = semla::GetSpatialNetwork(coords)
      spatgenes = semla::CorSpatialFeatures(current.expression.matrix,spatnet)
      top.svm = lapply(spatgenes,FUN = function(x)head(x$gene,input[['topNumber']]))
      names(top.svm)=unique(get_subset_exp()$meta[,colnames(bulk.metadata)[1]])
      return(list('top.peaks'=top.svm,'spatialPeaks'=spatgenes))
    })

    svm_upset <- reactive({
      top.svm <- run_svm()$top.peaks
      upset.plot <- UpSetR::upset(UpSetR::fromList(top.svm),nsets = length(names(top.svm)))
      return(upset.plot)
    })
    svm_barplot <- reactive({
 #     updateSelectInput(session, 'peakName', choices = svm_barplot()$names)
      top.svm <- run_svm()$top.peaks
      top.svm = bind_rows(top.svm, .id = "column_label")
      top.svm = tidyr::pivot_longer(top.svm,cols=colnames(top.svm))
      colnames(top.svm)[1]=colnames(bulk.metadata)[1]
      top.svm = merge(top.svm,bulk.metadata)
      plot = ggplot(top.svm,aes(x=forcats::fct_infreq(value),fill=get(input[['metadataBarplot']])))+
        geom_bar()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab('Peak')+
        scale_fill_discrete(name=input[['metadataBarplot']])
      return(list('plot'=plot,'names'=unique(forcats::fct_infreq(top.svm$value))))

    })

    observe({
      updateSelectInput(session, 'peakName', choices = svm_barplot()$names)
    })
    show_peak <- reactive({
      print(input[['peakName']])
      svm_results <- run_svm()
      my_peak = anno[anno$m_z==input[['peakName']],]
      print(my_peak)
      current.expression.matrix <- get_subset_exp()$exp
      current.metadata <- get_subset_exp()$meta
      caps = quantile(current.expression.matrix[,my_peak$m_z],probs=c(0.05,0.95))
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      current.metadata$peak = current.expression.matrix[,my_peak$m_z]
      current.metadata$peak = pmin(caps[2],current.metadata$peak)
      current.metadata$peak = pmax(caps[1],current.metadata$peak)
      spatial.plot = ggplot2::ggplot(current.metadata,ggplot2::aes(x=x,y=y,color=peak,fill=peak))+geom_tile()+
        ggplot2::facet_wrap(~current.metadata$Sample, scales = 'free',ncol=floor(2*sqrt(length(input[['samplesToInclude']]))))  +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(input[['peakName']]))+
        ggplot2::scale_fill_gradient(name=input[['peakName']],low = "lightgrey", high = "brown")+
        ggplot2::scale_color_gradient(name=input[['peakName']],low = "lightgrey", high = "brown")+
        ggplot2::theme(aspect.ratio = 1)
      return(spatial.plot)
    })

    output[['plotPeak']] <- renderPlot({
      show_peak()},height=600)
    output[['SVMUpset']] <- renderPlot({
      svm_upset()
    })%>% bindEvent(input[["run_SVM"]])
    output[['SVMTable']] <- DT::renderDT({
      svm.table = run_svm()$spatialPeaks[[1]]
      colnames(svm.table)=c('m_z','SVM_corr')
      svm.table = merge(svm.table,anno,all.x=T)
      svm.table = svm.table[order(-svm.table$SVM_corr),]
      svm.table
    })%>% bindEvent(input[["run_SVM"]])
    output[['SVMBarPlot']] <- renderPlot({
      svm_barplot()$plot
    })
    output$check <- reactive({
      length(input$samplesToInclude)
    })%>% bindEvent(input[["run_SVM"]])
    outputOptions(output, 'check', suspendWhenHidden=FALSE)
  })
}
