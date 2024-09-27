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
          tags$li("Peaks can be clustered to group by similar patterns using complete linkage hierarchical clustering based on Euclidean distance."),
          tags$li("Extra peaks can be added to the heatmap using the search box below.")
        ),
        br(),

        shinyjqui::orderInput(ns('heatmap.annotations'), label = "Show annotations", items = colnames(bulk.metadata[,sapply(bulk.metadata,n_distinct)!=nrow(bulk.metadata)])),
        checkboxInput(ns("cluster.heatmap"), label = "Cluster heatmap rows", value = TRUE),
        selectInput(ns("peakName"), "Additional peaks to include:", multiple = TRUE, choices = character(0)),
        div("\nIf no peaks are selected in the DE panel or here then the top 50 DE peaks are chosen.\n"),

        status = "success",
        icon = icon("gear", verify_fa = FALSE),
        tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
      ),
      dropMenu(
        circleButton(ns("downloads_heatmap"), icon = icon("download"),status = "success"),
        tags$div(
          tags$h3("Downloads"),
          fluidRow(
            column(10,offset=0,
                   textInput(ns('plotHeatmapFileName'), 'File name for heatmap download', value ='heatmap.png'),
                   numericInput(ns('heatmapPlotWidth'),value = 500,label = 'Width (px)',min = 50,max = 5000,step = 10),
                   numericInput(ns('heatmapPlotHeight'),value =800,label = 'Height (px)',min = 50,max = 5000,step = 10),
                   downloadButton(ns('downloadHeatmapPlot'), 'Download heatmap plot')
            )),
          theme = "light-border",
          placement = "right",
          arrow = FALSE
        )),
      plotly::plotlyOutput(ns('heatmap'), height = 800),
      tags$h1("Volcano/MA plot"),
      shinyWidgets::dropdownButton(
        tags$h3("Volcano/MA plots"),
        tags$ul(
          tags$li("Volcano/MA plots showing the log2FC, log10 BH-adjusted p-value and average intensity for each peak, colouring peaks showing significant changes."),
          tags$li("The y-axis scale for volcano plots can be capped at log10(p-value) > 10.")
        ),
        br(),
        selectInput(ns('plotType'), 'Type of plot:', c('Volcano', 'MA')),
        # shinyWidgets::switchInput(
        #   inputId = ns('autoLabel'),
        #   label = "Auto labels",
        #   labelWidth = "80px",
        #   onLabel = 'On',
        #   offLabel = 'Off',
        #   value = FALSE,
        #   onStatus = FALSE
        # ),
        # shinyWidgets::switchInput(
        #   inputId = ns("highlightSelected"),
        #   label = "Highlight selected DE peaks?",
        #   labelWidth = "80px",
        #   onLabel = 'No',
        #   offLabel = 'Yes',
        #   value = FALSE,
        #   onStatus = FALSE
        # ),
        # shinyWidgets::switchInput(
        #   inputId = ns('allPeaks'),
        #   label = "Showing on click:",
        #   labelWidth = "80px",
        #   onLabel = 'All peaks',
        #   offLabel = 'Only DE peaks',
        #   value = FALSE,
        #   onStatus = FALSE
        # ),
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
#        selectInput(ns("peakNameVolcano"), "Other peaks to highlight:", multiple = TRUE, choices = character(0)),

        status = "success",
        icon = icon("gear", verify_fa = FALSE),
        tooltip = shinyWidgets::tooltipOptions(title = "Click to see inputs!")
      ),
      dropMenu(
        circleButton(ns("downloads_volcano"), icon = icon("download"),status = "success"),
        tags$div(
          tags$h3("Downloads"),
          fluidRow(
            column(10,offset=0,
                   textInput(ns('plotFileNameVolcano'), 'File name for plot download', value ='DEPlot.png'),
                   numericInput(ns('volcanoPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                   numericInput(ns('volcanoPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                   downloadButton(ns('downloadVolcano'), 'Download Plot'),
            )),
          theme = "light-border",
          placement = "right",
          arrow = FALSE
        )),
      plotly::plotlyOutput(ns('volcanoPlot')),
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

    heatmap.prep <- reactive({
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
      scaled = t(scale(t(subsetExpression[, as.character(meta[, 1]), drop = FALSE])))
      return(scaled)
    })
    heatmap.plot <- reactive({
      print(input[['heatmap.annotations']])
      print(head(bulk.metadata[,colnames(bulk.metadata[,sapply(bulk.metadata,n_distinct)!=nrow(bulk.metadata)])]))
      scaled = heatmap.prep()
      peaks = rownames(scaled)
      peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
      meta <- lapply(bulk.metadata[,colnames(bulk.metadata[,sapply(bulk.metadata,n_distinct)!=nrow(bulk.metadata)])], function(x)if(!is.factor(x)){factor(x, levels = unique(x))}else{x}) |>
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

    updateSelectizeInput(session, "peakNameVolcano", choices = anno$m_z, server = TRUE)

    DEplot <- reactive({
      results = DEresults()$DE()
      selectedPeakNames = DEresults()$selectedPeaks()
      # if(!(input[["highlightSelected"]]) & length(selectedPeaks)){
      #   selectedPeakNames <- selectedPeaks
      #   highlightPeaks <- c(selectedPeakNames, input[["peakNameVolcano"]])
      # }
      # else{
      #   highlightPeaks <- input[["peakNameVolcano"]]
      # }

      if(input[['plotType']] == 'Volcano'){
        myplot <- volcano_plot(
          peaks.de.results = results$DEtable,
          pval.threshold = results$pvalThreshold,
          lfc.threshold = results$lfcThreshold,
          raster = TRUE,
          add.labels.auto = F,
          n.labels.auto = c(5, 5, 5),
          add.labels.custom = F,
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
          add.labels.auto = F,
          n.labels.auto = c(5, 5, 5),
          add.labels.custom = F,
          peaks.to.label = c()
        )
      }
      myplot
    })

    #Output MA/volcano plot
    output[['volcanoPlot']] <- plotly::renderPlotly(plotly::ggplotly(DEplot()))

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
      # for this phantomjs has to be available
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
          height = input[['heatmapPlotHeight']],
          width = input[['heatmapPlotWidth']]
        )
        rm(my.plot)
      }
    )

    output[['downloadVolcano']] <- downloadHandler(
      filename = function() { input[['plotFileNameVolcano']] },
      content = function(file) {
        ggsave(file, plot = DEplot(), width=input[['volcanoPlotWidth']],height=input[['volcanoPlotHeight']],units = 'in')
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
