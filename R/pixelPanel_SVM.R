#' @rdname IntroSpatialVisPanel
#' @export
PixelSVMPanelUI <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Spatial visualisation',
      dropMenu(
        circleButton(ns("info_svm"), icon = icon("info"),status = "success"),
        tags$div(
          tags$h3("Spatially variable metabolite identification"),
          tags$ul(
            tags$li("Spatially variable peaks/metabolites can be identified using the approach implemented in semla."),
            tags$li("The identification is run separately on each sample selected then the results are aggregated across samples to show the intersection of peaks which were among the top most spatially variable across samples in upset plots and barplots (where the bar can also be split by colour according to metadata information)."),
            tags$li("If differential intensity analysis has been run (in the Pseudobulk tab) then any significantly changing peaks are highlighted by stars above the barplot (if multiple samples are included) or in the table (when one sample only is included)."),
            tags$li("The top spatially variable peaks (the number of which is user-selected) can then be visualised spatially"),
            tags$li("Finally, the spatially variable peaks are clustered to identified modules of peaks which have similar intensity patterns"),
          )
        ),
        theme = "light-border",
        placement = "right",
        arrow = FALSE
      ),
      div(style = "margin-top:10px"),
      dropMenu(
        circleButton(ns("downloads"), icon = icon("download"),status = "success"),
        tags$div(
          tags$h3("Downloads"),
          fluidRow(
            column(5,offset=0,
                   tags$h4("Autocorrelation heatmap"),
                   textInput(ns('heatmapFileName'),'File name for download', value ='autocorPerSample.png', placeholder = 'autocorPerSample.png'),
                   numericInput(ns('heatmapWidth'),value = 500,label = 'Width (px)',min = 50,max = 5000,step = 10),
                   numericInput(ns('heatmapHeight'),value = 800,label = 'Height (px)',min = 50,max = 5000,step = 10),
                   downloadButton(ns('downloadHeatmap'), 'Download autocorrelation heatmap')
            ),
            column(5,offset=1,
                   tags$h4("Upset plot"),
                   textInput(ns('upsetFileName'),'File name for download', value ='upset.png', placeholder = 'autocorUpset.png'),
                   numericInput(ns('upsetWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                   numericInput(ns('upsetHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                   downloadButton(ns('downloadUpset'), 'Download upset plot')

            ),
          ),
          fluidRow(
            column(5,offset=0,
                   tags$h4("Barplot"),
                   textInput(ns('barFileName'),'File name for download', value ='autocorBar.png', placeholder = 'autocorBar.png'),
                   numericInput(ns('barWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                   numericInput(ns('barHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                   downloadButton(ns('downloadBar'), 'Download bar plot')
            ),
            column(5,offset=1,
                   tags$h4("Spatial cross-correlation heatmap"),
                   textInput(ns('crossCorFileName'),'File name for download', value ='crossCorHeatmap.png', placeholder = 'crossCorHeatmap.png'),
                   numericInput(ns('crossCorWidth'),value = 800,label = 'Width (px)',min = 50,max = 5000,step = 10),
                   numericInput(ns('crossCorHeight'),value = 800,label = 'Height (px)',min = 50,max = 5000,step = 10),
                   downloadButton(ns('downloadCrossCor'), 'Download cross-correlation heatmap')

            ),
          ),
          fluidRow(
            column(10,offset=0,
                   tags$h4("Spatial peak modules"),
                   textInput(ns('moduleFileName'),'File name for download', value ='spatialModules.png', placeholder = 'spatialModules.png'),
                   numericInput(ns('moduleWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                   numericInput(ns('moduleHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                   downloadButton(ns('downloadModule'), 'Download spatial modules')
            ),
          )
          ),
        theme = "light-border",
        placement = "right",
        arrow = FALSE
      ),

      sidebarLayout(
        sidebarPanel(
      numericInput(ns("topNumPeaksHeatmap"),label = 'Pick top n spatially variable metabolites',value = 20,min = 1,max = 100,step = 1),),
      mainPanel(
        tags$h3('Auto-correlation levels for top peaks'),
        plotly::plotlyOutput(ns('AutoCorHeatmap'),height='1000px'))),

      sidebarLayout(
        sidebarPanel(

          selectInput(
            inputId = ns("samplesToIncludeTopN"),
            label = "Select samples to show upset plot",
            choices = unique(bulk.metadata[,1]),
            selected = unique(bulk.metadata[,1]),
            multiple = TRUE,width = '100%'
          ),

          numericInput(ns("topNumber"),label = 'Pick top n spatially variable metabolites',value = 20,min = 1,max = 100,step = 1),
          # conditionalPanel(id=ns('upsetPlotPanel'),
          #                  ns = ns,
          #                  condition = "output.check > 1",
          selectInput(ns("metadataBarplot"), "Metadata to split barplot:", multiple = FALSE, choices = colnames(bulk.metadata),selected=colnames(bulk.metadata)[length(colnames(bulk.metadata))]),
          # ),
#          selectInput(ns("peakName"), "Peaks to display:", multiple = FALSE, choices = c()),
        ),
        mainPanel(
          tags$h3('Similarity in top peaks between selected samples'),
          # conditionalPanel(id=ns('upsetPlotPanel'),
          #   ns = ns,
          #   condition = "output.check > 1",
            plotOutput(ns('SVMUpset'),height=600),
            plotOutput(ns('SVMBarPlot'),height=600))),
          # conditionalPanel(id=ns('tablePanel'),
          #   ns = ns,
          #   condition = "output.check == 1",

    sidebarLayout(
      sidebarPanel(
      selectInput(
        inputId = ns("tableSample"),
        label = "Select samples to show top peaks",
        choices = unique(bulk.metadata[,1]),
        selected = unique(bulk.metadata[,1])[1],
        multiple = FALSE,width = '100%'
      )),
    mainPanel(
      tags$h3('Individual sample top peaks'),
      DT::DTOutput(ns('SVMTable')))),

    sidebarLayout(
      sidebarPanel(
        numericInput(ns('numClusters'),'Number of peak modules',min=2,max=20,step=1,value=10),
        numericInput(ns('selectedCluster'),'Cluster to show',min=2,max=20,step=1,value=10)),
      mainPanel(
        tags$h3('Spatial cross-correlation to define peak modules'),
        plotly::plotlyOutput(ns('CoexpHeatmap'),height='1000px'),
        textOutput(ns('selectedPeaks')),
        plotOutput(ns('CoexpSpatial'),height=600)))
#          fluidRow(column=10,plotOutput(ns('plotPeak'),height = 600)),
#          plotOutput(ns('plotHClust')),
    )
  }else{
    NULL
  }
}

#' @rdname IntroSpatialVisPanel
#' @export
PixelSVMPanelServer <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno, DEresults, svm_identification, spatial.cross.cor){

  moduleServer(id, function(input, output, session){

    autocor_heatmap <- reactive({
      svm.peaks = Reduce(union,lapply(FUN = function(x)head(x$gene,input[['topNumPeaksHeatmap']]),X=svm_identification))

      svm_identification_sub = lapply(svm_identification,function(x)x[x$gene%in%svm.peaks,])
      svm_identification_sub = dplyr::bind_rows(svm_identification_sub,.id = 'sample')
      svm_identification_sub_square = as.data.frame(tidyr::pivot_wider(svm_identification_sub,id_cols = 'sample',names_from = 'gene',values_from = 'cor'))
      rownames(svm_identification_sub_square)=svm_identification_sub_square$sample
      svm_identification_sub_square = as.data.frame(t(svm_identification_sub_square[,2:ncol(svm_identification_sub_square)]))
      peaks = rownames(svm_identification_sub_square)
      peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
      mat <- svm_identification_sub_square
      mat[] <- peaks
      return(heatmaply::heatmaply_cor(svm_identification_sub_square,limits = c(min(svm_identification_sub_square),1),custom_hovertext=mat))
    })
    # get_subset_exp <- reactive({
    #   current.intensity.matrix <- full.intensity.matrix[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToInclude']],]
    #   current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToInclude']],]
    #   current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
    #   rownames(current.intensity.matrix)=current.metadata[,1]
    #   return(list('exp'=current.intensity.matrix,'meta'=current.metadata))
    # }) #%>% bindEvent(input[["run_SVM"]])

    svm_results <- reactive({
      load('svm_identification.rda')
      spatgenes <- svm_identification[input[['samplesToInclude']]]
      return(spatgenes)
    }) #%>% bindEvent(input[["run_SVM"]])

    run_svm <- reactive({
      spatgenes <- svm_results()
      top.svm = lapply(spatgenes,FUN = function(x)head(x$gene,input[['topNumber']]))
      return(top.svm)
    })

    svm_upset <- reactive({
      svm.peaks = Reduce(union,lapply(FUN = function(x)head(x$gene,input[['topNumber']]),X=svm_identification[input[['samplesToIncludeTopN']]]))

      top.svm = lapply(FUN = function(x)x[x$gene%in%head(x$gene,input[['topNumber']]),],X=svm_identification[input[['samplesToIncludeTopN']]])
      top.svm = lapply(top.svm,function(x)x$gene[x$gene%in%svm.peaks])
      upset.plot <- UpSetR::upset(UpSetR::fromList(top.svm),nsets = length(names(top.svm)))
      return(upset.plot)
    })

    svm_barplot <- reactive({
 #     updateSelectInput(session, 'peakName', choices = svm_barplot()$names)
#      top.svm <- run_svm()
      svm.peaks = Reduce(union,lapply(FUN = function(x)head(x$gene,input[['topNumPeaksHeatmap']]),X=svm_identification[input[['samplesToIncludeTopN']]]))
      top.svm = lapply(FUN = function(x)x[x$gene%in%head(x$gene,input[['topNumPeaksHeatmap']]),],X=svm_identification[input[['samplesToIncludeTopN']]])
      top.svm = lapply(top.svm,function(x)x$gene[x$gene%in%svm.peaks])
      top.svm = bind_rows(top.svm, .id = "column_label")
      top.svm = tidyr::pivot_longer(top.svm,cols=colnames(top.svm))
      colnames(top.svm)[1]=colnames(bulk.metadata)[1]
      top.svm = merge(top.svm,bulk.metadata)
      if (DEresults()$runDE==1){
        deres = DEresults()$DE()$DEtable[,c('m_z','pvalAdj','lfc')]
        deres$DE = deres$m_z %in% DEresults()$DE()$DEtableSubset$m_z
        colnames(deres)[1]='value'
        top.svm = merge(top.svm,deres,all.x=T)
        top.svm$label = ifelse(top.svm$DE==TRUE,'*','')
      } else {
        top.svm$label=''
      }
      top.svm$metadataBarplot = top.svm[,input[['metadataBarplot']]]
      top.svm = top.svm %>% group_by(value,metadataBarplot,label) %>% dplyr::summarise(n=n()) %>% arrange(desc(n))
      top.svm.grouped = top.svm %>% group_by(value,label) %>% dplyr::summarise(n=sum(n)) %>% arrange(desc(n))
      top.svm$value = factor(top.svm$value,levels=top.svm.grouped$value)
      plot = ggplot(data=top.svm,aes(x=value, fill=metadataBarplot,y=n)) +
        geom_bar(stat='identity')+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab('Peak')+
        scale_fill_discrete(name=input[['metadataBarplot']])+
        ylab('Number of samples')+
        geom_text(data = top.svm.grouped,aes(x=value,y=n,label = label,fill=NULL),size=10)

      return(list('plot'=plot,'names'=top.svm.grouped$value))
    })

    coexp_heatmap <- reactive({
      peaks = rownames(spatial.cross.cor)
      peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
      peaks.combo = expand.grid(peaks,peaks,stringsAsFactors = F)
      peaks.combo = paste0('row:',peaks.combo$Var1,',\ncol:',peaks.combo$Var2)
      mat <- spatial.cross.cor
      mat[] <- peaks.combo
      return(heatmaply::heatmaply_cor(spatial.cross.cor,custom_hovertext = mat,k_col=input[['numClusters']]))
    })

    cexp_clusters <- reactive({
      d <- dist(spatial.cross.cor, method = "euclidean")
      clusters = cutree(hclust(d), k = input[['numClusters']])
      dend <- as.dendrogram(hclust(d, method = "complete"))
      dend <- dendextend::seriate_dendrogram(dend, d)
      clusters = clusters[rev(rownames(spatial.cross.cor[order.dendrogram(dend),]))]
      cluster.names = names(clusters)
      clusters = as.numeric(factor(clusters,levels=unique(clusters)))
      names(clusters)=cluster.names
      return(clusters)
    })
    selected_peaks <- reactive({
      clusters = cexp_clusters()
      return(names(clusters[clusters==input[['selectedCluster']]]))
    })

    output[['selectedPeaks']] <- renderText({paste('Selected peaks:',paste(selected_peaks(),collapse = ', '))})

    cexp_cluster_spatial <- reactive({
      clusters = cexp_clusters()
      cluster.peaks = names(clusters[clusters==input[['selectedCluster']]])
      current.metadata = full.metadata
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      current.metadata$combined_cluster = rowMeans(scale(full.intensity.matrix[,cluster.peaks]))
      spatial.plot = ggplot2::ggplot(current.metadata,ggplot2::aes(x=x,y=y,color=combined_cluster,fill=combined_cluster))+geom_tile()+
        ggplot2::facet_wrap(~current.metadata$Sample, scales = 'free',ncol=floor(2*sqrt(length(unique(current.metadata$Sample)))))  +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank())+
        ggplot2::scale_fill_gradient2(name=paste0('Cluster ',input[['selectedCluster']]))+
        ggplot2::scale_color_gradient2(name=paste0('Cluster ',input[['selectedCluster']]))+
        ggplot2::theme(aspect.ratio = 1)
      return(spatial.plot)
    })

    observe({
      updateSelectInput(session, 'peakName', choices = svm_barplot()$names)
      updateNumericInput(session, 'selectedCluster', min = 1,max = input[['numClusters']],value = 1,step=1)
    })  #%>% bindEvent(input[["run_SVM"]])

    svmTable <- reactive({
      svm.table = svm_identification[[input[['tableSample']]]]
      colnames(svm.table)=c('m_z','SVM_corr')
      svm.table = merge(svm.table,anno[,c('m_z','name')],all.x=T)
      if (DEresults()$runDE==1){
        deres = DEresults()$DE()$DEtable[,c('m_z','pvalAdj','lfc')]
        deres$DE = deres$m_z %in% DEresults()$DE()$DEtableSubset$m_z
        svm.table = merge(svm.table,deres,all.x=T)
      }
      svm.table = svm.table[order(-svm.table$SVM_corr),]
      svm.table
    }) #%>% bindEvent(input[["run_SVM"]])

    output[['AutoCorHeatmap']] <- plotly::renderPlotly({
      autocor_heatmap()
    })

    output[['CoexpHeatmap']] <- plotly::renderPlotly({
      coexp_heatmap()
    })
    output[['SVMUpset']] <- renderPlot({
      svm_upset()
    },height=600)

    output[['SVMTable']] <- DT::renderDT({
      svmTable()
    })

    output[['SVMBarPlot']] <- renderPlot({
      svm_barplot()$plot
    },height=600)


    output[['CoexpSpatial']] <- renderPlot({
      cexp_cluster_spatial()},height=600)

    output[['downloadBar']] <- downloadHandler(
      filename = function() { input[['barFileName']] },
      content = function(file) {
        ggsave(file, plot = svm_barplot()$plot, dpi = 300,
               width=input[['barWidth']],height=input[['barHeight']])
      }
    )

    output[['downloadModule']] <- downloadHandler(
      filename = function() { input[['moduleFileName']] },
      content = function(file) {
        ggsave(file, plot = cexp_cluster_spatial(), dpi = 300,
               width=input[['moduleWidth']],height=input[['moduleHeight']])
      }
    )


    output[['downloadUpset']] <- downloadHandler(
      filename = function() { input[['upsetFileName']] },
      content = function(file) {
        if (base::strsplit(input[['upsetFileName']], split="\\.")[[1]][-1] == 'pdf'){
          grDevices::pdf(file, width = input[['upsetWidth']], height = input[['upsetHeight']])
          print(svm_upset())
          grDevices::dev.off()
        } else if (base::strsplit(input[['upsetFileName']], split="\\.")[[1]][-1] == 'svg'){
          grDevices::svg(file, width = input[['upsetWidth']], height = input[['upsetHeight']])
          print(svm_upset())
          grDevices::dev.off()
        } else {
          grDevices::png(file, width = input[['upsetWidth']], height = input[['upsetHeight']], units = "in",
                         res = 300, bg = "white")
          print(svm_upset())
          grDevices::dev.off()
        }
      }
    )

    output[['downloadHeatmap']] <- downloadHandler(
      # for this phantomjs has to be available
      filename = function() { input[['heatmapFileName']] },
      content = function(file) {
        svm.peaks = Reduce(union,lapply(FUN = function(x)head(x$gene,input[['topNumPeaksHeatmap']]),X=svm_identification))

        svm_identification_sub = lapply(svm_identification,function(x)x[x$gene%in%svm.peaks,])
        svm_identification_sub = dplyr::bind_rows(svm_identification_sub,.id = 'sample')
        svm_identification_sub_square = as.data.frame(tidyr::pivot_wider(svm_identification_sub,id_cols = 'sample',names_from = 'gene',values_from = 'cor'))
        rownames(svm_identification_sub_square)=svm_identification_sub_square$sample
        svm_identification_sub_square = as.data.frame(t(svm_identification_sub_square[,2:ncol(svm_identification_sub_square)]))
        peaks = rownames(svm_identification_sub_square)
        peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
        mat <- svm_identification_sub_square
        mat[] <- peaks
        print(head(mat))
        my.plot = heatmaply::heatmaply_cor(
          svm_identification_sub_square,
          limits = c(min(svm_identification_sub_square),1),
          custom_hovertext=mat,
          height = input[['heatmapHeight']],
          width = input[['heatmapWidth']],
          file = file
          )
        rm(my.plot)
        }
    )

    output[['downloadCrossCor']] <- downloadHandler(
      # for this phantomjs has to be available
      filename = function() { input[['crossCorFileName']] },
      content = function(file) {
        peaks = rownames(spatial.cross.cor)
        peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
        peaks.combo = expand.grid(peaks,peaks,stringsAsFactors = F)
        peaks.combo = paste0('row:',peaks.combo$Var1,',\ncol:',peaks.combo$Var2)
        mat <- spatial.cross.cor
        mat[] <- peaks.combo
        my.plot = heatmaply::heatmaply_cor(spatial.cross.cor,
                                           custom_hovertext = mat,
                                           k_col=input[['numClusters']],
                                           height = input[['crossCorHeight']],
                                           width = input[['crossCorWidth']],
                                           file = file)

        rm(my.plot)
      }
    )



  })

}
