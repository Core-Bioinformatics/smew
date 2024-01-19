#' @rdname BulkDESummaryPanel
#' @export
BulkDESummaryPanelUI <- function(id, bulk.metadata, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Differential intensity analysis visualisation',
      tags$h1("Peak heatmap"),
      shinyWidgets::dropdownButton(
        tags$h3("Peak heatmap"),
        tags$ul(
          tags$li("Heatmap showing variation of intensity across samples for peaks selected in Differential analysis tab."),
          tags$li("Raw intensities, log2 intensities of Z-score intensities (default) can be shown for each peak."),
          tags$li("Peaks can be clustered to group by similar patterns using complete linkage hierarchical clustering based on Euclidean distance."),
          tags$li("Extra peaks can be added to the heatmap using the search box below."),
        ),
        br(),

        radioButtons(ns('heatmap.processing'), label = "Heatmap values",
                     choices = c('Expression','Log2 Expression','Z-score'),
                     selected = 'Z-score'),
        shinyjqui::orderInput(ns('heatmap.annotations'), label = "Show annotations", items = colnames(bulk.metadata)),
        checkboxInput(ns("cluster.heatmap"), label = "Cluster heatmap rows", value = TRUE),
        selectInput(ns("peakName"), "Additional peaks to include:", multiple = TRUE, choices = character(0)),
        div("\nIf no peaks are selected in the DE panel or here then the top 50 DE peaks are chosen.\n"),
        div(style="margin-bottom:10px"),
        textInput(ns('plotHeatmapFileName'), 'File name for heatmap plot download', value ='HeatmapPlot.png'),
        downloadButton(ns('downloadHeatmapPlot'), 'Download Heatmap Plot'),

        status = "info",
        icon = icon("gear", verify_fa = FALSE),
        tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
      ),
      plotOutput(ns('heatmap'), height = 800),
      tags$h1("Volcano/MA plot"),
      shinyWidgets::dropdownButton(
        tags$h3("Volcano/MA plots"),
        tags$ul(
          tags$li("Volcano/MA plots showing the log2FC, log10 BH-adjusted p-value and average intensity for each peak, colouring peaks showing significant changes."),
          tags$li("Peaks selected in Differential intensity analysis tab are highlighted."),
          tags$li("Extra peaks can be highlighted using the search box below."),
          tags$li("The y-axis scale for volcano plots can be capped at log10(p-value) > 10."),
        ),
        br(),
        selectInput(ns('plotType'), 'Type of plot:', c('Volcano', 'MA')),
        shinyWidgets::switchInput(
          inputId = ns('autoLabel'),
          label = "Auto labels",
          labelWidth = "80px",
          onLabel = 'On',
          offLabel = 'Off',
          value = FALSE,
          onStatus = FALSE
        ),
        shinyWidgets::switchInput(
          inputId = ns("highlightSelected"),
          label = "Highlight selected DE peaks?",
          labelWidth = "80px",
          onLabel = 'No',
          offLabel = 'Yes',
          value = FALSE,
          onStatus = FALSE
        ),
        shinyWidgets::switchInput(
          inputId = ns('allPeaks'),
          label = "Showing on click:",
          labelWidth = "80px",
          onLabel = 'All peaks',
          offLabel = 'Only DE peaks',
          value = FALSE,
          onStatus = FALSE
        ),
        conditionalPanel(
          id = ns('conditionalVolcanoOption'),
          ns=ns,
          condition = "input[['plotType']] == 'Volcano'",
          shinyWidgets::switchInput(
            inputId = ns("capPVal"),
            label = "Cap log10(pval)?",
            labelWidth = "80px",
            onLabel = 'No',
            offLabel = 'Yes',
            value = FALSE,
            onStatus = FALSE
          ),
        ),
        selectInput(ns("peakNameVolcano"), "Other peaks to highlight:", multiple = TRUE, choices = character(0)),
        textInput(ns('plotFileNameVolcano'), 'File name for plot download', value ='DEPlot.png'),
        downloadButton(ns('downloadVolcano'), 'Download Plot'),

        status = "info",
        icon = icon("gear", verify_fa = FALSE),
        tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
      ),
      plotOutput(ns('volcanoPlot'), click = ns('plot_click')),
      tableOutput(ns('volcanoData'))

    )
  }else{
    NULL
  }
}

#' @rdname BulkDESummaryPanel
#' @export
BulkDESummaryPanelServer <- function(id, bulk.intensity.matrix, bulk.metadata, DEresults, anno){

  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(DEresults)
    is.reactive(bulk.intensity.matrix)
    is.reactive(bulk.metadata)
    !is.reactive(anno)
  })

  moduleServer(id, function(input, output, session){

    #Set up server-side search for peak names
    updateSelectizeInput(session, "peakName", choices = anno$m_z, server = TRUE)

    observe({
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
    heatmap.plot <- reactive({
      selectedPeaks = DEresults()$selectedPeaks()
      if(length(selectedPeaks)){
        selectedPeakNames <- selectedPeaks
        peakSet <- c(selectedPeakNames, input[["peakName"]])
      }else{
        peakSet <- input[["peakName"]]
      }
      if (length(peakSet) == 0){
        peakSet <- head(DEresults()$DE()$DEtableSubset$m_z, 50)
      }
      peakIDs <- peakSet
      subsetExpression <- bulk.intensity.matrix[peakIDs, , drop = FALSE]
      rownames(subsetExpression) <- peakSet
      meta <- lapply(bulk.metadata, function(x)if(!is.factor(x)){factor(x, levels = unique(x))}else{x}) |>
        as.data.frame() |>
        dplyr::arrange(dplyr::across(input[['heatmap.annotations']]))
      myplot <- expression_heatmap_met(
        intensity.matrix.subset = subsetExpression[, as.character(meta[, 1]), drop = FALSE],
        top.annotation.ids = match(input[['heatmap.annotations']], colnames(meta)),
        metadata = meta,
        type = input[["heatmap.processing"]],
        show.column.names = (nrow(meta) <= 20),
        cluster.peaks = input[['cluster.heatmap']]
      )
      return(myplot)
    })
    output[['heatmap']] <- renderPlot(heatmap.plot(), height = 800)

    updateSelectizeInput(session, "peakNameVolcano", choices = anno$m_z, server = TRUE)

    DEplot <- reactive({
      results = DEresults()$DE()
      selectedPeaks = DEresults()$selectedPeaks()
      if(!(input[["highlightSelected"]]) & length(selectedPeaks)){
        selectedPeakNames <- selectedPeaks
        highlightPeaks <- c(selectedPeakNames, input[["peakNameVolcano"]])
      }
      else{
        highlightPeaks <- input[["peakNameVolcano"]]
      }

      if(input[['plotType']] == 'Volcano'){
        myplot <- volcano_plot(
          peaks.de.results = results$DEtable,
          pval.threshold = results$pvalThreshold,
          lfc.threshold = results$lfcThreshold,
          raster = TRUE,
          add.labels.auto = input[["autoLabel"]],
          n.labels.auto = c(5, 5, 5),
          add.labels.custom = length(highlightPeaks) > 0,
          peaks.to.label = highlightPeaks,
          log10pval.cap = !(input[['capPVal']])
        )
      }
      if (input[['plotType']] == 'MA'){
        myplot <- ma_plot(
          peaks.de.results = results$DEtable,
          pval.threshold = results$pvalThreshold,
          lfc.threshold = results$lfcThreshold,
          raster = TRUE,
          add.labels.auto = input[["autoLabel"]],
          n.labels.auto = c(5, 5, 5),
          add.labels.custom = length(highlightPeaks) > 0,
          peaks.to.label = highlightPeaks
        )
      }
      myplot
    })

    #Output MA/volcano plot
    output[['volcanoPlot']] <- renderPlot(DEplot())

    #Define output table when you click on peak with all peaks or only DE
    output[['volcanoData']] <- renderTable({
      req(input[['plot_click']])
      results = DEresults()$DE()
      if (input[['allPeaks']]){
        data <- results$DEtable
      }else{
        data <- results$DEtableSubset
      }
      data <- data |> dplyr::mutate(`-log10pval` = -log10(.data$pvalAdj))
      nearPoints(df = data, coordinfo = input[['plot_click']], threshold = 20, maxpoints = 10)
    }, digits = 4)

    output[['downloadHeatmapPlot']] <- downloadHandler(
      filename = function() { input[['plotHeatmapFileName']] },
      content = function(file) {
        if (base::strsplit(input[['plotHeatmapFileName']], split="\\.")[[1]][-1] == 'pdf'){
          grDevices::pdf(file, width = 10, height = 20, pointsize = 20)
          print(heatmap.plot())
          grDevices::dev.off()
        } else if (base::strsplit(input[['plotHeatmapFileName']], split="\\.")[[1]][-1] == 'svg'){
          grDevices::svg(file, width = 10, height = 20, pointsize = 20)
          print(heatmap.plot())
          grDevices::dev.off()
        } else {
          grDevices::png(file, width = 480, height = 1000, units = "px",
                         pointsize = 12, bg = "white", res = NA)
          print(heatmap.plot())
          grDevices::dev.off()
        }
      }
    )

    output[['downloadVolcano']] <- downloadHandler(
      filename = function() { input[['volcanoPlotFileName']] },
      content = function(file) {
        ggsave(file, plot = DEplot(), dpi = 300)
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
