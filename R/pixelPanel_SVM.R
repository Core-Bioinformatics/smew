#' Identifies modules of spatially variable and colocalised metabolite peaks
#'
#' @description UI and server logic for spatially variable peak analysis using SVM and spatial correlation in the SMEW app. Visualizes spatial autocorrelation, cross-correlation, and top spatially variable peaks across samples.
#'
#' @details
#' \itemize{
#'  \item{Visualize spatial autocorrelation heatmaps for top spatially variable metabolites}
#'  \item{Upset plots and barplots for summarizing spatially variable peaks across samples}
#'  \item{Table of top spatially variable peaks per sample}
#'  \item{Cross-correlation modules for peaks with similar spatial patterns}
#'  \item{Download handlers for all plots}
#'  \item{Interactive controls for sample and peak selection}
#' }
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param full.metadata Data frame of full sample metadata
#' @param full.intensity.matrix Matrix of intensities (features x samples)
#' @param show Logical; whether to show the panel (default TRUE)
#' @param spatial.cross.cor List or data structure with spatial cross-correlation results
#'
#' @return UI: A shiny::tabPanel object for the SVM panel. Server: None; called for side effects in Shiny module.
#' @export
#' @name PixelPanel_SVMTab
#' @rdname PixelPanel_SVMTab
PixelPanel_SVMTabUI <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, show = TRUE){
  ns <- shiny::NS(id)
  
  if(show){
    shiny::tabPanel(
      'Spatial Auto/Cross-Correlation',
      shiny::tags$h2('Spatial auto-correlation heatmap'),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_autocor"), icon = shiny::icon("info"), status = "info"),
              shiny::tags$div(
                shiny::tags$h3("Auto-correlation heatmap info"),
                shiny::tags$ul(
                  shiny::tags$li("Shows the spatial autocorrelation for the top n spatially variable metabolites across all samples."),
                  shiny::tags$li("Use the controls to select the number of peaks.")
                )
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloadsHeatmap"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download autocorrelation heatmap"),
                shiny::textInput(ns('heatmapFileName'), 'File name for download', value = 'autocorPerSample.png'),
                shiny::numericInput(ns('heatmapWidth'), value = 500, label = 'Width (px)', min = 50, max = 5000, step = 10),
                shiny::numericInput(ns('heatmapHeight'), value = 800, label = 'Height (px)', min = 50, max = 5000, step = 10),
                shiny::downloadButton(ns('downloadHeatmap'), 'Download autocorrelation heatmap')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::numericInput(ns("topNumPeaksHeatmap"), label = 'Pick top n spatially variable metabolites', value = 20, min = 1, max = 100, step = 1)
        ),
        shiny::mainPanel(
          plotly::plotlyOutput(ns('AutoCorHeatmap'), height = '1000px')
        )
      ),
      
      shiny::tags$h2('Summarising spatially auto-correlated peaks across samples'),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_upset"), icon = shiny::icon("info"), status = "info"),
              shiny::tags$div(
                shiny::tags$h3("Upset plot & barplot info"),
                shiny::tags$ul(
                  shiny::tags$li("The upset plot shows the intersection of top spatially variable peaks across selected samples."),
                  shiny::tags$li("The intersection points are ordered by number of overlapping samples."),
                  shiny::tags$li("The barplot shows the number of samples where each peak is significant. The Barplot can be split by metadata and differential analysis information from bulk tabs is overlaid.")
                )
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloadsUpsetBar"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download upset plot & barplot"),
                shiny::textInput(ns('upsetFileName'), 'File name for upset plot', value = 'upset.png'),
                shiny::numericInput(ns('upsetWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                shiny::numericInput(ns('upsetHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                shiny::downloadButton(ns('downloadUpset'), 'Download upset plot'),
                shiny::textInput(ns('barFileName'), 'File name for barplot', value = 'autocorBar.png'),
                shiny::numericInput(ns('barWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                shiny::numericInput(ns('barHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                shiny::downloadButton(ns('downloadBar'), 'Download bar plot')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::selectInput(
            inputId = ns("samplesToIncludeTopN"),
            label = "Select samples to show upset plot",
            choices = unique(bulk.metadata[, 1]),
            selected = unique(bulk.metadata[, 1]),
            multiple = TRUE, width = '100%'
          ),
          shiny::numericInput(ns("topNumber"), label = 'Pick top n spatially variable metabolites', value = 20, min = 1, max = 100, step = 1),
          shiny::selectInput(ns("metadataBarplot"), "Metadata to split barplot:", multiple = FALSE, choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[length(colnames(bulk.metadata))])
        ),
        shiny::mainPanel(
          shiny::plotOutput(ns('SVMUpset'), height = 600),
          shiny::plotOutput(ns('SVMBarPlot'), height = 600)
        )
      ),
      
      shiny::tags$h2('Top auto-correlated peaks per sample'),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_table"), icon = shiny::icon("info"), status = "info"),
              shiny::tags$div(
                shiny::tags$h3("Individual sample top peaks info"),
                shiny::tags$ul(
                  shiny::tags$li("Shows the top spatially variable peaks for the selected sample.")
                )
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::selectInput(
            inputId = ns("tableSample"),
            label = "Select samples to show top peaks",
            choices = unique(bulk.metadata[, 1]),
            selected = unique(bulk.metadata[, 1])[1],
            multiple = FALSE, width = '100%'
          )
        ),
        shiny::mainPanel(
          DT::DTOutput(ns('SVMTable'))
        )
      ),
      shiny::tags$h2('Spatial cross-correlation modules'),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_crosscor"), icon = shiny::icon("info"), status = "info"),
              shiny::tags$div(
                shiny::tags$h3("Spatial cross-correlation info"),
                shiny::tags$ul(
                  shiny::tags$li("Shows modules of peaks with similar spatial patterns across samples.")
                )
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloadsCrossCor"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download cross-correlation heatmap & modules"),
                shiny::textInput(ns('crossCorFileName'), 'File name for cross-correlation heatmap', value = 'crossCorHeatmap.png'),
                shiny::numericInput(ns('crossCorWidth'), value = 800, label = 'Width (px)', min = 50, max = 5000, step = 10),
                shiny::numericInput(ns('crossCorHeight'), value = 800, label = 'Height (px)', min = 50, max = 5000, step = 10),
                shiny::downloadButton(ns('downloadCrossCor'), 'Download cross-correlation heatmap'),
                shiny::textInput(ns('moduleFileName'), 'File name for spatial modules', value = 'spatialModules.png'),
                shiny::numericInput(ns('moduleWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                shiny::numericInput(ns('moduleHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                shiny::downloadButton(ns('downloadModule'), 'Download spatial modules')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::numericInput(ns('numClusters'), 'Number of peak modules', min = 2, max = 20, step = 1, value = 10),
          shiny::numericInput(ns('selectedCluster'), 'Cluster to show', min = 1, max = 10, step = 1, value = 1)
        ),
        shiny::mainPanel(
          plotly::plotlyOutput(ns('CoexpHeatmap'), height = '1000px'),
          shiny::plotOutput(ns('CoexpSpatial'), height = 600),
          shiny::tags$h4('Peaks in selected module'),
          shiny::uiOutput(ns('selectedModulePeaks'))
        )
      ),
      shiny::tags$h2('Spatial cross-correlation network'),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_network"), icon = shiny::icon("info"), status = "info"),
              shiny::tags$div(
                shiny::tags$h3("Spatial cross-correlation network info"),
                shiny::tags$ul(
                  shiny::tags$li("Shows a network of peaks with strong spatial cross-correlation. Node color highlights the selected module if enabled.")
                )
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloadsNetwork"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download network plot"),
                shiny::textInput(ns('networkFileName'), 'File name for network plot', value = 'crosscorNetwork.html'),
                shiny::downloadButton(ns('downloadNetwork'), 'Download network (HTML)')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::radioButtons(ns('edgeMode'), 'Edge selection mode',
            choices = c('Top edges overall' = 'overall', 'Top edges per node' = 'pernode'),
            selected = 'overall', inline = TRUE),
          shiny::conditionalPanel(
            sprintf("input['%s'] == 'overall'", ns('edgeMode')),
            shiny::numericInput(ns('numNetworkEdgesOverall'), 'Number of top edges to show', min = 10, max = 200, value = 40, step = 1)
          ),
          shiny::conditionalPanel(
            sprintf("input['%s'] == 'pernode'", ns('edgeMode')),
            shiny::numericInput(ns('numNetworkEdgesPerNode'), 'Number of top edges per node', min = 1, max = 20, value = 3, step = 1)
          ),
          shiny::checkboxInput(ns('highlightModule'), 'Highlight selected module', value = TRUE)
        ),
        shiny::mainPanel(
          visNetwork::visNetworkOutput(ns('ccNet'), height = '900px')
        )
      )
    )
  }else{
    NULL
  }
}

#' @param id Shiny module id
#' @param anno Data frame of peak annotations
#' @param de.results Reactive expression or function returning DE results
#' @param svm.identification List of SVM identification results per sample
#' @param spatial.cross.cor Matrix of spatial cross-correlation values
#' @return None; called for side effects in Shiny module
#' @export
#' @rdname PixelPanel_SVMTab
#' @export
PixelPanel_SVMTabServer <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno, de.results, svm.identification, spatial.cross.cor){
  
  shiny::moduleServer(id, function(input, output, session){
    
    autocor_heatmap <- shiny::reactive({
      svm.peaks = Reduce(union,lapply(FUN = function(x)utils::head(x$peak,input[['topNumPeaksHeatmap']]),X=svm.identification))
      
      svm.identification_sub = lapply(svm.identification,function(x)x[x$peak%in%svm.peaks,])
      svm.identification_sub = dplyr::bind_rows(svm.identification_sub,.id = 'sample')
      svm.identification_sub_square = as.data.frame(tidyr::pivot_wider(svm.identification_sub,id_cols = 'sample',names_from = 'peak',values_from = 'cor'))
      rownames(svm.identification_sub_square)=svm.identification_sub_square$sample
      svm.identification_sub_square = as.data.frame(t(svm.identification_sub_square[,2:ncol(svm.identification_sub_square)]))
      peaks = rownames(svm.identification_sub_square)
      peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
      mat <- svm.identification_sub_square
      mat[] <- peaks
      return(heatmaply::heatmaply_cor(svm.identification_sub_square,limits = c(min(svm.identification_sub_square),1),custom_hovertext=mat))
    })
    
    svm_results <- shiny::reactive({
      load('svm.identification.rda')
      spatgenes <- svm.identification[input[['samplesToInclude']]]
      return(spatgenes)
    }) #|> bindEvent(input[["run_SVM"]])
    
    run_svm <- shiny::reactive({
      spatgenes <- svm_results()
      top.svm = lapply(spatgenes,FUN = function(x)utils::head(x$peak,input[['topNumber']]))
      return(top.svm)
    })
    
    svm_upset <- shiny::reactive({
      svm.peaks = Reduce(union,lapply(FUN = function(x)utils::head(x$peak,input[['topNumber']]),X=svm.identification[input[['samplesToIncludeTopN']]]))
      
      top.svm = lapply(FUN = function(x)x[x$peak%in%utils::head(x$peak,input[['topNumber']]),],X=svm.identification[input[['samplesToIncludeTopN']]])
      top.svm = lapply(top.svm,function(x)x$peak[x$peak%in%svm.peaks])
      upset.plot <- UpSetR::upset(UpSetR::fromList(top.svm),nsets = length(names(top.svm)))
      return(upset.plot)
    })
    
    svm_barplot <- shiny::reactive({
      svm.peaks = Reduce(union,lapply(FUN = function(x)utils::head(x$peak,input[['topNumPeaksHeatmap']]),X=svm.identification[input[['samplesToIncludeTopN']]]))
      top.svm = lapply(FUN = function(x)x[x$peak%in%utils::head(x$peak,input[['topNumPeaksHeatmap']]),],X=svm.identification[input[['samplesToIncludeTopN']]])
      top.svm = lapply(top.svm,function(x)x$peak[x$peak%in%svm.peaks])
      top.svm = dplyr::bind_rows(top.svm, .id = "column_label")
      top.svm = tidyr::pivot_longer(top.svm,cols=colnames(top.svm))
      colnames(top.svm)[1]=colnames(bulk.metadata)[1]
      top.svm = merge(top.svm,bulk.metadata)
      if (de.results()$runDE==1){
        deres = de.results()$DE()$DEtable[,c('m_z','pvalAdj','lfc')]
        deres$DE = deres$m_z %in% de.results()$DE()$DEtableSubset$m_z
        colnames(deres)[1]='value'
        top.svm = merge(top.svm,deres,all.x=T)
        top.svm$label = ifelse(top.svm$DE==TRUE,'*','')
      } else {
        top.svm$label=''
      }
      top.svm$metadataBarplot = top.svm[,input[['metadataBarplot']]]
      top.svm = top.svm |> dplyr::group_by(.data$value,.data$metadataBarplot,.data$label) |> dplyr::summarise(n=dplyr::n(), .groups = 'drop') |> dplyr::arrange(dplyr::desc(.data$n))
      top.svm.grouped = top.svm |> dplyr::group_by(.data$value,.data$label) |> dplyr::summarise(n=sum(.data$n), .groups = 'drop') |> dplyr::arrange(dplyr::desc(.data$n))
      top.svm$value = factor(top.svm$value,levels=top.svm.grouped$value)
      plot = ggplot2::ggplot(data=top.svm,ggplot2::aes(x=.data$value, fill=.data$metadataBarplot,y=.data$n)) +
        ggplot2::geom_bar(stat='identity')+
        ggplot2::theme_classic()+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))+
        ggplot2::xlab('Peak')+
        ggplot2::scale_fill_discrete(name=input[['metadataBarplot']])+
        ggplot2::ylab('Number of samples')+
        ggplot2::geom_text(data = top.svm.grouped,ggplot2::aes(x=.data$value,y=.data$n,label = .data$label,fill=NULL),size=10)
      
      return(list('plot'=plot,'names'=top.svm.grouped$value))
    })
    
    coexp_heatmap <- shiny::reactive({
      peaks = rownames(spatial.cross.cor)
      peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
      peaks.combo = expand.grid(peaks,peaks,stringsAsFactors = F)
      peaks.combo = paste0('row:',peaks.combo$Var1,',\ncol:',peaks.combo$Var2)
      mat <- spatial.cross.cor
      mat[] <- peaks.combo
      return(heatmaply::heatmaply_cor(spatial.cross.cor,custom_hovertext = mat,k_col=input[['numClusters']]))
    })
    
    cexp_clusters <- shiny::reactive({
      d <- stats::dist(spatial.cross.cor, method = "euclidean")
      clusters = stats::cutree(stats::hclust(d), k = input[['numClusters']])
      dend <- stats::as.dendrogram(stats::hclust(d, method = "complete"))
      dend <- dendextend::seriate_dendrogram(dend, d)
      clusters = clusters[rev(rownames(spatial.cross.cor[stats::order.dendrogram(dend),]))]
      cluster.names = names(clusters)
      clusters = as.numeric(factor(clusters,levels=unique(clusters)))
      names(clusters)=cluster.names
      return(clusters)
    })
    
    crosscor_network <- shiny::reactive({
      edges = spatial.cross.cor
      edges$from = rownames(spatial.cross.cor)
      edges = tidyr::pivot_longer(edges,cols = 1:(ncol(edges)-1))
      colnames(edges)=c('from','to','value')
      edges = edges |> dplyr::filter(.data$from<.data$to)
      if (input$edgeMode == 'overall') {
        n_edges <- input$numNetworkEdgesOverall
        edges = edges |> dplyr::arrange(dplyr::desc(.data$value)) |> utils::head(n_edges)
      } else {
        n_edges <- input$numNetworkEdgesPerNode
        edges = edges |> dplyr::group_by(.data$from) |> dplyr::arrange(dplyr::desc(.data$value)) |> dplyr::slice_head(n = n_edges) |> dplyr::ungroup()
      }
      all_peaks <- rownames(spatial.cross.cor)
      node_df = data.frame(id=all_peaks, label=all_peaks, shape='box')
      clusters = cexp_clusters()
      sel_mod <- input$selectedCluster
      if (isTRUE(input$highlightModule)) {
        in_selected <- clusters[as.character(node_df$id)] == sel_mod
        node_df$color.background = ifelse(in_selected, '#4A90E2', '#e0e0e0')
        node_df$color.border = ifelse(in_selected, 'black', '#b0b0b0')
        node_df$font.color = ifelse(in_selected, 'black', '#b0b0b0')
      } else {
        node_df$color.background = '#e0e0e0'
        node_df$color.border = '#b0b0b0'
        node_df$font.color = 'black'
      }
      if (nrow(edges) > 0) edges$width <- 1 else edges$width <- numeric(0)
      set.seed(1234) # Ensures deterministic layout
      visNetwork::visNetwork(node_df,edges) |>
        visNetwork::visIgraphLayout(layout = 'layout_nicely') |>
        visNetwork::visPhysics(enabled = FALSE) |>
        visNetwork::visNodes(font = list(size = 28))
    })
    
    output[['ccNet']] = visNetwork::renderVisNetwork(crosscor_network())
    selected_peaks <- shiny::reactive({
      clusters = cexp_clusters()
      return(names(clusters[clusters==input[['selectedCluster']]]))
    })
    
    output$selectedModulePeaks <- shiny::renderUI({
      clusters <- cexp_clusters()
      sel_mod <- input$selectedCluster
      if (is.null(clusters) || is.null(sel_mod)) return(NULL)
      peaks <- names(clusters)[clusters == sel_mod]
      if (length(peaks) == 0) return(shiny::tags$em('No peaks in this module.'))
      shiny::tags$div(
        style = 'padding: 10px; background: #f8f8f8; border-radius: 8px; border: 1px solid #ddd; margin-bottom: 20px; max-height: 100px; overflow-y: auto; font-size: 1.05em; color: #333;',
        paste(peaks, collapse = ', ')
      )
    })
    
    cexp_cluster_spatial <- shiny::reactive({
      clusters = cexp_clusters()
      cluster.peaks = names(clusters[clusters==input[['selectedCluster']]])
      current.metadata = full.metadata
      current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
      current.metadata$combined_cluster = rowMeans(scale(full.intensity.matrix[,cluster.peaks]))
      spatial.plot = ggplot2::ggplot(current.metadata,ggplot2::aes(x=.data$x_tf,y=.data$y_tf,color=.data$combined_cluster,fill=.data$combined_cluster))+ggplot2::geom_tile()+
        ggplot2::facet_wrap(~current.metadata$Sample,ncol=floor(1.5*sqrt(length(unique(current.metadata$Sample)))))  +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank())+
        ggplot2::scale_fill_viridis_c(option = 'A', name=paste0('Cluster ',input[['selectedCluster']]))+
        ggplot2::scale_color_viridis_c(option = 'A', name=paste0('Cluster ',input[['selectedCluster']]))+
        ggplot2::theme(aspect.ratio = 1)
      return(spatial.plot)
    })
    
    shiny::observe({
      shiny::updateSelectInput(session, 'peakName', choices = svm_barplot()$names)
      shiny::updateNumericInput(session, 'selectedCluster', min = 1, max = input[['numClusters']], value = min(input[['selectedCluster']], input[['numClusters']]))
    })  #|> bindEvent(input[["run_SVM"]])
    
    svmTable <- shiny::reactive({
      svm.table = svm.identification[[input[['tableSample']]]]
      colnames(svm.table)=c('m_z','SVM_corr')
      svm.table = merge(svm.table,anno[,c('m_z','name')],all.x=T)
      if (de.results()$runDE==1){
        deres = de.results()$DE()$DEtable[,c('m_z','pvalAdj','lfc')]
        deres$DE = deres$m_z %in% de.results()$DE()$DEtableSubset$m_z
        svm.table = merge(svm.table,deres,all.x=T)
      }
      svm.table = svm.table[order(-svm.table$SVM_corr),]
      svm.table
    }) #|> bindEvent(input[["run_SVM"]])
    
    output[['AutoCorHeatmap']] <- plotly::renderPlotly({
      autocor_heatmap()
    })
    
    output[['CoexpHeatmap']] <- plotly::renderPlotly({
      coexp_heatmap()
    })
    output[['SVMUpset']] <- shiny::renderPlot({
      svm_upset()
    },height=600)
    
    output[['SVMTable']] <- DT::renderDT({
      svmTable()
    })
    
    output[['SVMBarPlot']] <- shiny::renderPlot({
      svm_barplot()$plot
    },height=600)
    
    
    output[['CoexpSpatial']] <- shiny::renderPlot({
      cexp_cluster_spatial()},height=600)
    
    output[['downloadBar']] <- shiny::downloadHandler(
      filename = function() { input[['barFileName']] },
      content = function(file) {
        ggplot2::ggsave(file, plot = svm_barplot()$plot, dpi = 300,
               width=input[['barWidth']],height=input[['barHeight']])
      }
    )
    
    output[['downloadModule']] <- shiny::downloadHandler(
      filename = function() { input[['moduleFileName']] },
      content = function(file) {
        ggplot2::ggsave(file, plot = cexp_cluster_spatial(), dpi = 300,
               width=input[['moduleWidth']],height=input[['moduleHeight']])
      }
    )
    
    
    output[['downloadUpset']] <- shiny::downloadHandler(
      filename = function() { input[['upsetFileName']] },
      content = function(file) {
        ext <- tolower(tools::file_ext(input[['upsetFileName']]))
        ggplot2::ggsave(filename = file, plot = ggplotify::as.ggplot(svm_upset()), device = ext,
               width = input[['upsetWidth']], height = input[['upsetHeight']],
               units = "in", dpi = 300)
      }
    )
    
    output[['downloadHeatmap']] <- shiny::downloadHandler(
      # for this phantomjs has to be available
      filename = function() { input[['heatmapFileName']] },
      content = function(file) {
        svm.peaks = Reduce(union,lapply(FUN = function(x)utils::head(x$peak,input[['topNumPeaksHeatmap']]),X=svm.identification))
        
        svm.identification_sub = lapply(svm.identification,function(x)x[x$peak%in%svm.peaks,])
        svm.identification_sub = dplyr::bind_rows(svm.identification_sub,.id = 'sample')
        svm.identification_sub_square = as.data.frame(tidyr::pivot_wider(svm.identification_sub,id_cols = 'sample',names_from = 'peak',values_from = 'cor'))
        rownames(svm.identification_sub_square)=svm.identification_sub_square$sample
        svm.identification_sub_square = as.data.frame(t(svm.identification_sub_square[,2:ncol(svm.identification_sub_square)]))
        peaks = rownames(svm.identification_sub_square)
        peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
        mat <- svm.identification_sub_square
        mat[] <- peaks
        my.plot = heatmaply::heatmaply_cor(
          svm.identification_sub_square,
          limits = c(min(svm.identification_sub_square),1),
          custom_hovertext=mat,
          height = input[['heatmapHeight']],
          width = input[['heatmapWidth']],
          file = file
        )
        rm(my.plot)
      }
    )
    
    output[['downloadCrossCor']] <- shiny::downloadHandler(
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
    
    
    output[['downloadNetwork']] <- shiny::downloadHandler(
      filename = function() {
        fn <- input[['networkFileName']]
        if (!grepl("\\.html$", fn, ignore.case = TRUE)) fn <- paste0(fn, ".html")
        fn
      },
      content = function(file) {
        net <- crosscor_network() # This is an htmlwidget
        htmlwidgets::saveWidget(net, file, selfcontained = TRUE)
      }
    )
    
  })
  
}
