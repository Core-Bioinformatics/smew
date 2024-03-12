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

      selectInput(
        inputId = ns("samplesToInclude"),
        label = "Select samples to include",
        choices = unique(bulk.metadata[,1]),
        selected = unique(bulk.metadata[,1])[1],
        multiple = TRUE,width = '100%'
      ),
      tags$p("Warning: this analysis can take a few minutes to run and this will increase with the number of samples included."),
      actionButton(ns("run_SVM"),'Identify spatially variable metabolites'),
      sidebarLayout(
        sidebarPanel(
          numericInput(ns("topNumber"),label = 'Pick top n spatially variable metabolites',value = 20,min = 1,max = 100,step = 1),
          conditionalPanel(id=ns('upsetPlotPanel'),
                           ns = ns,
                           condition = "output.check > 1",
                           selectInput(ns("metadataBarplot"), "Metadata to split barplot:", multiple = FALSE, choices = colnames(bulk.metadata),selected=colnames(bulk.metadata)[length(colnames(bulk.metadata))]),
          ),
          selectInput(ns("peakName"), "Peaks to display:", multiple = FALSE, choices = c()),
          div(style = "margin-top:10px"),
          dropMenu(
            circleButton(ns("downloads"), icon = icon("download"),status = "success"),
            tags$div(
              tags$h3("Downloads"),
              fluidRow(
                conditionalPanel(id=ns('multipleSamplePlots'),
                                 ns = ns,
                                 condition = "output.check > 1",
                                 column(5,offset=0,
                                        tags$h4("Barplot"),
                                        textInput(ns('barPlotFileName'),'File name for download', value ='SVMbar.png', placeholder = 'SVMbar.png'),
                                        numericInput(ns('barPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                                        numericInput(ns('barPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                                        downloadButton(ns('downloadBar'), 'Download ORA table')),
                                 column(5,offset=0,
                                        tags$h4("Upset"),
                                        textInput(ns('upsetFileName'),'File name for download', value ='SVMupset.png', placeholder = 'SVMupset.png'),
                                        numericInput(ns('upsetPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                                        numericInput(ns('upsetPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                                        downloadButton(ns('downloadUpsetPlot'), 'Download ORA table'))),

                 conditionalPanel(id=ns('oneSamplePlots'),
                                  ns = ns,
                                  condition = "output.check == 1",
                                  column(10,offset=0,
                                         textInput(ns('tableFileName'),'File name for download', value ='ORATable.csv', placeholder = 'ORATable.csv'),
                                         downloadButton(ns('downloadTable'), 'Download ORA table')))),
              fluidRow(
                column(5,offset=0,
                       tags$h4("Spatial distribution"),
                       textInput(ns('spatialFileName'),'File name for download', value ='spatial.png', placeholder = 'spatial.png'),
                       numericInput(ns('spatialWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('spatialHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadSpatial'), 'Download spatial figure')
                ),
                column(5,offset=1,
                       tags$h4("Hierarchical Gene Clustering"),
                       textInput(ns('HClustFileName'),'File name for download', value ='hclust.png', placeholder = 'hclust.png'),
                       numericInput(ns('HClustWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('HClustHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadHClust'), 'Download clustering')

                ),

              )),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
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
            DT::DTOutput(ns('SVMTable')),
          ))),
          fluidRow(column=10,plotOutput(ns('plotPeak'),height = 600)),
          plotOutput(ns('plotHClust')),
    )
  }else{
    NULL
  }
}

#' @rdname IntroSpatialVisPanel
#' @export
PixelSVMPanelServer <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno, DEresults){

  moduleServer(id, function(input, output, session){

    get_subset_exp <- reactive({
      current.intensity.matrix <- full.intensity.matrix[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToInclude']],]
      current.metadata <- full.metadata[full.metadata[,colnames(bulk.metadata)[1]] %in% input[['samplesToInclude']],]
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      rownames(current.intensity.matrix)=current.metadata[,1]
      return(list('exp'=current.intensity.matrix,'meta'=current.metadata))
    }) %>% bindEvent(input[["run_SVM"]])

    svm_results <- reactive({
      load('svm_identification.rda')
      spatgenes <- svm_identification[input[['samplesToInclude']]]
      return(spatgenes)
    }) %>% bindEvent(input[["run_SVM"]])

    run_svm <- reactive({
      spatgenes <- svm_results()
      top.svm = lapply(spatgenes,FUN = function(x)head(x$gene,input[['topNumber']]))
      return(top.svm)
    })

    svm_upset <- reactive({
      top.svm <- run_svm()
      upset.plot <- UpSetR::upset(UpSetR::fromList(top.svm),nsets = length(names(top.svm)))
      return(upset.plot)
    })

    svm_barplot <- reactive({
 #     updateSelectInput(session, 'peakName', choices = svm_barplot()$names)
      top.svm <- run_svm()
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

    observe({
      updateSelectInput(session, 'peakName', choices = svm_barplot()$names)
    })  %>% bindEvent(input[["run_SVM"]])

    show_peak <- reactive({
      svm_results <- run_svm()
      my_peak = anno[anno$m_z==input[['peakName']],]
      current.intensity.matrix <- get_subset_exp()$exp
      current.metadata <- get_subset_exp()$meta
      caps = quantile(current.intensity.matrix[,my_peak$m_z],probs=c(0.05,0.95))
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      current.metadata$peak = current.intensity.matrix[,my_peak$m_z]
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
                       axis.line.y = ggplot2::element_blank())+
        ggplot2::scale_fill_gradient(name=input[['peakName']],low = "lightgrey", high = "brown")+
        ggplot2::scale_color_gradient(name=input[['peakName']],low = "lightgrey", high = "brown")+
        ggplot2::theme(aspect.ratio = 1)
      return(spatial.plot)
    })

    hclust.spatgenes <- reactive({
      top.svm <- run_svm()
      top.svm = bind_rows(top.svm, .id = "column_label")
      top.svm = tidyr::pivot_longer(top.svm,cols=colnames(top.svm))
      top.svm = unique(top.svm$value)
      full.matrix = get_subset_exp()$exp[,top.svm]
      full.matrix = scale(full.matrix,center = T,scale=T)
      d = dist(t(full.matrix))
      hclust(d)
    })

    output[['SVMUpset']] <- renderPlot({
      svm_upset()
    })

    output[['downloadUpsetPlot']] <- downloadHandler(
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

    svmTable <- reactive({
      svm.table = svm_results()[[input[['samplesToInclude']]]]
      colnames(svm.table)=c('m_z','SVM_corr')
      svm.table = merge(svm.table,anno[,c('m_z','name')],all.x=T)
      if (DEresults()$runDE==1){
        deres = DEresults()$DE()$DEtable[,c('m_z','pvalAdj','lfc')]
        deres$DE = deres$m_z %in% DEresults()$DE()$DEtableSubset$m_z
        svm.table = merge(svm.table,deres,all.x=T)
      }
      svm.table = svm.table[order(-svm.table$SVM_corr),]
      svm.table
    }) %>% bindEvent(input[["run_SVM"]])

    output[['SVMTable']] <- DT::renderDT({
      svmTable()
    })

    output[['downloadTable']] <- downloadHandler(
      filename = function() {
        paste(input[['tableFileName']])
      },
      content = function(file) {
        utils::write.csv(x = svmTable(), file = file, row.names = FALSE)
      }
    )

    output[['SVMBarPlot']] <- renderPlot({
      svm_barplot()$plot
    })

    output[['downloadBar']] <- downloadHandler(
      filename = function() { input[['barPlotFileName']] },
      content = function(file) {
        ggsave(file, plot = svm_barplot()$plot, dpi = 300,
               width=input[['barPlotWidth']],height=input[['barPlotHeight']])
      }
    )

    output[['plotPeak']] <- renderPlot({
      show_peak()},height=600)

    output[['downloadSpatial']] <- downloadHandler(
      filename = function() { input[['spatialFileName']] },
      content = function(file) {
          ggsave(file, plot = show_peak(), dpi = 300,
                 width=input[['spatialWidth']],height=input[['spatialHeight']],
                 device=base::strsplit(input[['spatialFileName']], split="\\.")[[1]][-1])
      }
    )
    output$check <- reactive({
      length(input$samplesToInclude)
    })%>% bindEvent(input[["run_SVM"]])

    outputOptions(output, 'check', suspendWhenHidden=FALSE)

    output[['plotHClust']] <- renderPlot({
      plot(hclust.spatgenes())})

    output[['downloadHClust']] <- downloadHandler(
      filename = function() { input[['HClustFileName']] },
      content = function(file) {
        if (base::strsplit(input[['HClustFileName']], split="\\.")[[1]][-1] == 'pdf'){
          grDevices::pdf(file, width = input[['HClustWidth']], height = input[['HClustHeight']])
          plot(hclust.spatgenes())
          grDevices::dev.off()
        } else if (base::strsplit(input[['HClustFileName']], split="\\.")[[1]][-1] == 'svg'){
          grDevices::svg(file, width = input[['HClustWidth']], height = input[['HClustHeight']])
          plot(hclust.spatgenes())
          grDevices::dev.off()
        } else {
          grDevices::png(file, width = input[['HClustWidth']], height = input[['HClustHeight']], units = "in",
                         res = 300, bg = "white")
          plot(hclust.spatgenes())
          grDevices::dev.off()
        }
      }
    )


  })

}
