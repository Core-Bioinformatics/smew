#' Differential Expression Summary Panel
#' @description Visual summaries of results from the Differential intensity analysis module.
#' Provides heatmap of selected peaks, volcano/MA plots and interactive selection tables.
#'
#' @details The panel pulls peak sets from the linked DE analysis module. Heatmap
#' rows can be clustered and annotated using sample metadata. Volcano/MA visual
#' summarise differential peaks with convenient download options.
#'
#' @param id Module id
#' @param bulk.metadata A data.frame with sample-level metadata
#' @param show Logical flag whether to include this tab in the UI (default TRUE)
#' @keywords internal
#' @name de_summary_panel
#' @rdname de_summary_panel
#' @export
BulkDESummaryUI <- function(id, bulk.metadata, show = TRUE){
  ns <- shiny::NS(id)

  if(show){
    shiny::tabPanel(
      'Differential intensity analysis visualisation',
      shiny::tags$h1("Peak heatmap"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_heatmap"), icon = shiny::icon("info-circle"), status = "info"),
              shiny::tags$h3("Heatmap"),
              shiny::tags$ul(
                shiny::tags$li("The heatmap shows intensity variation for selected peaks across samples."),
                shiny::tags$li("Rows can be clustered to group similar patterns."),
                shiny::tags$li("Annotations highlight sample metadata groups and can be reordered according to user choice."),
                shiny::tags$li("Use the settings to tailor clustering options and select additional peaks, beyond those selected in the DA tab.")
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_heatmap"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download Heatmap"),
                shiny::textInput(ns('plotHeatmapFileName'), 'File name for heatmap download', value ='heatmap.png'),
                shiny::numericInput(ns('heatmapPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('heatmapPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadHeatmapPlot'), 'Download heatmap plot')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shinyjqui::orderInput(ns('heatmap.annotations'), label = "Show annotations", items = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)])),
          shiny::checkboxInput(ns("cluster.heatmap"), label = "Cluster heatmap rows", value = TRUE),
          shiny::selectInput(ns("peakName"), "Additional peaks to include:", multiple = TRUE, choices = character(0)),
          shiny::div("If no peaks are selected in the DE panel or here then the top 50 DE peaks are chosen."),

        ),
        shiny::mainPanel(
          plotly::plotlyOutput(ns('heatmap'), height = 800)
        )
      ),
      shiny::tags$h1("Volcano/MA plot"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_volcano"), icon = shiny::icon("info-circle"), status = "info"),
              shiny::tags$h3("Volcano/MA Plot"),
              shiny::tags$ul(
                shiny::tags$li("The volcano plot shows log2 fold change vs. -log10 BH-adjusted p-value for each peak."),
                shiny::tags$li("Significant peaks are shown through colour."),
                shiny::tags$li("The MA plot shows average intensity vs. log2 fold change."),
                shiny::tags$li("Use the settings to select plot type and further customisation options.")
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_volcano"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download Volcano/MA Plot"),
                shiny::textInput(ns('plotFileNameVolcano'), 'File name for plot download', value ='DEPlot.png'),
                shiny::numericInput(ns('volcanoPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('volcanoPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadVolcano'), 'Download Plot')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::selectInput(ns('plotType'), 'Type of plot:', c('Volcano', 'MA')),
          shiny::conditionalPanel(
            id = ns('conditionalVolcanoOption'),
            ns=ns,
            condition = "input[['plotType']] == 'Volcano'",
            shinyWidgets::switchInput(
              inputId = ns("capPVal"),
              label = "Cap log10(pval) at 10?",
              labelWidth = "80px",
              onLabel = 'No',
              offLabel = 'Yes',
              value = FALSE,
              onStatus = FALSE
            ),
            ),

        ),
        shiny::mainPanel(
          plotly::plotlyOutput(ns('volcanoPlot')),
          shiny::tableOutput(ns('volcanoData'))
        )
      )

    )
  }else{
    NULL
  }
}

# ============ SERVER FUNCTIONS ============

#' @rdname de_summary_panel
#' @param bulk.intensity.matrix Reactive matrix of peak intensities (rows = peaks, cols = samples)
#' @param DEresults Reactive object returned by the DE analysis module (provides DEtable, DEtableSubset, selectedPeaks())
#' @param anno Annotation data.frame for peaks (must contain m_z and display_name/name fields)
#' @export
BulkDESummaryServer <- function(id, bulk.intensity.matrix, bulk.metadata, DEresults, anno){

  # check whether inputs (other than id) are reactive or not
  stopifnot({
    shiny::is.reactive(DEresults)
    shiny::is.reactive(bulk.intensity.matrix)
    shiny::is.reactive(bulk.metadata)
    !shiny::is.reactive(anno)
  })

  shiny::moduleServer(id, function(input, output, session){

    #Set up server-side search for peak names
    shiny::updateSelectizeInput(session, "peakName", choices = anno$m_z, server = TRUE)

    shiny::observe({
      items <- colnames(bulk.metadata)
      include.exclude <- apply(bulk.metadata, 2, function(x){
        l <- length(unique(x))
        (l > 1) & (l < length(x))
      })
      if (sum(include.exclude == TRUE) != 0){
        items <- colnames(bulk.metadata)[include.exclude]
        items <- items[c(length(items), seq_len(length(items) - 1))]
      } else {items = colnames(bulk.metadata)[2:ncol(bulk.metadata)]}
      shinyjqui::updateOrderInput(session, "heatmap.annotations", items = items)
    })

    heatmap.prep <- shiny::reactive({
      selectedPeaks = DEresults()$selectedPeaks()
      if(length(selectedPeaks)){
        selectedPeakNames <- selectedPeaks
        peakSet <- c(selectedPeakNames, input[["peakName"]])
      }else{
        peakSet <- input[["peakName"]]
      }
      if (length(peakSet) == 0){
        peakSet <- utils::head(DEresults()$DE()$DEtableSubset$m_z, 50)
      }
      peakIDs <- peakSet
      subsetExpression <- bulk.intensity.matrix[peakIDs, , drop = FALSE]
      rownames(subsetExpression) <- peakSet
      meta <- lapply(bulk.metadata, function(x)if(!is.factor(x)){factor(x, levels = unique(x))}else{x}) |>
        as.data.frame() |>
        dplyr::arrange(dplyr::across(input[['heatmap.annotations']]))
      scaled = t(scale(t(subsetExpression[, as.character(meta[, 1]), drop = FALSE])))
      return(scaled)
    })
    heatmap.plot <- shiny::reactive({
      scaled = heatmap.prep()
      peaks = rownames(scaled)
      peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
      meta <- lapply(bulk.metadata[,colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)])], function(x)if(!is.factor(x)){factor(x, levels = unique(x))}else{x}) |>
        as.data.frame() |>
        dplyr::arrange(dplyr::across(input[['heatmap.annotations']]))
      mat <- scaled
      mat[] <- peaks
      return(heatmaply::heatmaply_cor(
        scaled,
        Colv = FALSE,
        Rowv = input[['cluster.heatmap']],
        limits = c(-max(abs(scaled)),max(abs(scaled))),
        col_side_colors = meta,
        height=800,
        custom_hovertext = mat
      ))
    })

      # myplot <- expression_heatmap_met(
      #   intensity.matrix.subset = subsetExpression[, as.character(meta[, 1]), drop = FALSE],
      #   top.annotation.ids = match(input[['heatmap.annotations']], colnames(meta)),
      #   metadata = meta,
      #   type = input[["heatmap.processing"]],
      #   show.column.names = (nrow(meta) <= 20),
      #   cluster.peaks = input[['cluster.heatmap']]
      # )
#      return(myplot)

    output[['heatmap']] <- plotly::renderPlotly(heatmap.plot())

    shiny::updateSelectizeInput(session, "peakNameVolcano", choices = anno$m_z, server = TRUE)

    DEplot <- shiny::reactive({
      results = DEresults()$DE()
      selectedPeakNames = DEresults()$selectedPeaks()
      if(input[['plotType']] == 'Volcano'){
        myplot <- volcano_plot(
          peaks.de.results = results$DEtable,
          pval.threshold = results$pvalThreshold,
          lfc.threshold = results$lfcThreshold,
          raster = TRUE,
          add.labels.auto = FALSE,
          n.labels.auto = c(5, 5, 5),
          add.labels.custom = FALSE,
          peaks.to.label = c(),
          log10pval.cap = !(input[['capPVal']])
        )
      }
      if (input[['plotType']] == 'MA'){
        myplot <- ma_plot(
          peaks.de.results = results$DEtable,
          pval.threshold = results$pvalThreshold,
          lfc.threshold = results$lfcThreshold,
          raster = TRUE,
          add.labels.auto = FALSE,
          n.labels.auto = c(5, 5, 5),
          add.labels.custom = FALSE,
          peaks.to.label = c()
        )
      }
      myplot
    })

    #Output MA/volcano plot
    output[['volcanoPlot']] <- plotly::renderPlotly(plotly::ggplotly(DEplot()))

    #Define output table when you click on peak with all peaks or only DE
    output[['volcanoData']] <- shiny::renderTable({
      shiny::req(input[['plot_click']])
      results = DEresults()$DE()
      if (input[['allPeaks']]){
        data <- results$DEtable
      }else{
        data <- results$DEtableSubset
      }
      data <- data |> dplyr::mutate(`-log10pval` = -log10(.data$pvalAdj))
      shiny::nearPoints(df = data, coordinfo = input[['plot_click']], threshold = 20, maxpoints = 10)
    }, digits = 4)

    output[['downloadHeatmapPlot']] <- shiny::downloadHandler(
      filename = function() { input[['plotHeatmapFileName']] },
      content = function(file) {
        scaled = heatmap.prep()
        peaks = rownames(scaled)
        peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
        meta <- lapply(bulk.metadata, function(x)if(!is.factor(x)){factor(x, levels = unique(x))}else{x}) |>
          as.data.frame() |>
          dplyr::arrange(dplyr::across(input[['heatmap.annotations']]))
        mat <- scaled
        mat[] <- peaks
        my.plot = heatmaply::heatmaply_cor(
          scaled,
          Colv = FALSE,
          Rowv = input[['cluster.heatmap']],
          limits = c(-max(abs(scaled)),max(abs(scaled))),
          col_side_colors = meta[,2:ncol(meta)],
          custom_hovertext = mat,
          file = file,
          height = input[['heatmapPlotHeight']]*300,
          width = input[['heatmapPlotWidth']]*300
        )
        rm(my.plot)
      }
    )

    output[['downloadVolcano']] <- shiny::downloadHandler(
      filename = function() { input[['plotFileNameVolcano']] },
      content = function(file) {
        ggplot2::ggsave(file, plot = DEplot(), width=input[['volcanoPlotWidth']],height=input[['volcanoPlotHeight']],units = 'in')
      }
    )


  })
}

#DEsummaryPanelApp <- function(){
#  shinyApp(
#    ui = navbarPage("DE", tabPanel("", tabsetPanel(DEpanelUI('RNA'), DEsummaryPanelUI('RNA')))),
#    server = function(input, output, session){
#      DEresults <- DEpanelServer('RNA')
#      DEsummaryPanelServer('RNA', DEresults)
#    }
#  )
#}
