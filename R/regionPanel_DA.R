#' @rdname RegionDEPanel
#' @export
RegionDEpanelUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Differential intensity analysis',
      shinyjs::useShinyjs(),
      sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(
          dropMenu(
            circleButton(ns("info_differential_analysis"), icon = icon("info"),status = "success"),
            tags$div(
              tags$h3("Differential intensity analysis"),
              tags$ul(
                tags$li("First select a segmentation to pseudobulk on, which is either the clusters identified in the previous tab or a column from the metadata table which is not sample-wide. For each region, a sample will be created for each original sample assuming there are sufficient observations (the threshold can be altered)."),
                tags$li("A sample-wide metadata column should also be selected "),
                tags$li("Comparisons are only performed once the button has been pressed."),
                tags$li("The most recent comparison is passed onto other tabs for visualisation, pathway analysis etc."),
                tags$li("The table of peaks showing significant changes can also be downloaded as a csv."),
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),

          selectInput(
            inputId = ns("regionToGroupOn"),
            label = "Select region to group on",
            choices = c(colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y','Sample',colnames(bulk.metadata)))]),
            selected = c(colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y','Sample',colnames(bulk.metadata)))])[1]
          ),
          numericInput(inputId = ns('minimumPixelsPerSample'),
                       label = 'Minimum number of pixels per region per sample to be considered',
                       value = 10,
                       min = 1,
                       max = 1000),
          selectInput(ns('condition'), 'Metadata column to use:', colnames(bulk.metadata)[-1],
                      selected = colnames(bulk.metadata)[ncol(bulk.metadata)]),
          actionButton(ns('pseudoBulk'), label = 'Start Pseudobulking'),

          # Input: Selector variables to compare
          selectInput(ns('variable1'), 'Condition 1:', c()),
          selectInput(ns('variable2'), 'Condition 2:', c(),
                      selected = c()),

          selectInput(ns('pipeline'), 'DE pipeline:', c("t-test", "Wilcox rank sum")),

          #DE thresholds

          sliderInput(ns('pvalThreshold'), label = 'Adjusted p-value threshold',
                      min = 0, value = 0.05, max = 1, step = 0.005),
          sliderInput(ns('lfcThreshold'), label = 'log2 fold change threshold',
                      min = 0, value = 1, max = 5, step = 0.1),

          #Only start DE when button is pressed
          actionButton(ns('goDE'), label = 'Start differential intensity analysis'),

          #download file name and button
          textInput(ns('fileName'),'File name for download', value ='DIAset.csv', placeholder = 'DIAset.csv'),
          downloadButton(ns('download'), 'Download Table'),
          hr(),
          tags$b("Peak selection"),
          div("\nSelect peaks of interest by clicking on the corresponds rows in the table\n"),
          div(style="margin-bottom:10px"),
          actionButton(ns('resetSelection'), label = "Reset row selection"),
          div(style="margin-bottom:10px"),
          actionButton(ns('selectTop50'), label = "Select top 50 peaks")

        ),

        #Main panel for displaying table of DE peaks
        mainPanel(
          textOutput(ns('print_pseudobulk')),
          DT::DTOutput(ns('data'))
        )
      )
    )
  }else{
    NULL
  }
}

#' @rdname RegionDEPanel
#' @export
RegionDEpanelServer <- function(id, full.intensity.matrix, full.metadata, bulk.intensity.matrix, bulk.metadata, region.clusters, anno){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(bulk.intensity.matrix)
    is.reactive(bulk.metadata)
    !is.reactive(anno)
  })

  moduleServer(id, function(input, output, session){
    observe(updateSelectInput(
      session,
      'regionToGroupOn',
      choices = colnames(region.clusters())[!(colnames(region.clusters())%in%c('spot_id','x','y','Sample',colnames(bulk.metadata)))],
      selected=colnames(region.clusters())[!(colnames(region.clusters())%in%c('spot_id','x','y','Sample',colnames(bulk.metadata)))][1])
    )

    pseudobulk_samples <- reactive({
      pseudobulked <- create_bulk_exp_regions(as.data.frame(full.intensity.matrix),
                              region.clusters(),
                              bulk.metadata,
                              colnames(bulk.metadata)[1],
                              region.clusters()[,input[['regionToGroupOn']]],
                              minimum.pixels = input[['minimumPixelsPerSample']])
#      if (is.null(region.clusters())){
        # choices =
        # selected = colnames(full.metadata)[colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata))][1]
#      } else {
#        choices = c('cluster',colnames(full.metadata)[colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata))])
#        selected = 'cluster'
#      }
      updateSelectInput(session, 'variable1', choices = unique(pseudobulked$metadata[[input[["condition"]]]]))
      updateSelectInput(session, 'variable2', choices = unique(pseudobulked$metadata[[input[["condition"]]]]),
                        selected = unique(pseudobulked$metadata[[input[["condition"]]]])[2])

      return(pseudobulked)
    }) %>% bindEvent(input[["pseudoBulk"]])

    output[['print_pseudobulk']]<-renderText(paste0('Conditions extracted: ',paste(unique(pseudobulk_samples()$metadata[[input[["condition"]]]]),collapse = ', ')))

    observe({
      condition.indices <- region.clusters()[[input[["condition"]]]] %in% c(input[['variable1']], input[['variable2']])
      choices <- c("t-test", "Wilcox rank sum")
      updateSelectInput(session, 'pipeline', choices = choices)
    })


    DEresults <- reactive({
      shinyjs::disable("goDE")
      pseudobulk.metadata = pseudobulk_samples()$metadata
      pseudobulk.intensity.matrix = pseudobulk_samples()$intensity.matrix
      pseudobulk.intensity.matrix = pseudobulk.intensity.matrix[,pseudobulk.metadata[,1]]
      condition.indices <- pseudobulk.metadata[[input[["condition"]]]] %in% c(input[['variable1']], input[['variable2']])
      # Need to add error if any group has <2 samples
      DEtable <- DEanalysis(
        intensity.matrix = pseudobulk.intensity.matrix[, condition.indices],
        condition = pseudobulk.metadata[[input[["condition"]]]][condition.indices],
        var1 = input[['variable1']],
        var2 = input[['variable2']],
        test = input[["pipeline"]],
        anno = anno
      )
      DEtableSubset <- DEtable %>%
        dplyr::filter(.data$pvalAdj < input[["pvalThreshold"]] & abs(.data$lfc) > input[['lfcThreshold']]) %>%
        dplyr::arrange(.data$pvalAdj)

      #the thresholds are returned here so that MA/volcano and table display
      #don't use new thresholds without the button being used
      shinyjs::enable("goDE")
      return(list('DEtable' = DEtable,
                  "DEtableSubset" = DEtableSubset,
                  'pvalThreshold' = input[["pvalThreshold"]],
                  'lfcThreshold' = input[['lfcThreshold']]))
    }) %>%
      bindCache(region.clusters, input[["condition"]],
                input[['variable1']], input[['variable2']], input[["pipeline"]],
                input[["pvalThreshold"]],input[['lfcThreshold']]) %>%
      bindEvent(input[["goDE"]])

    #Define output table (only DE peaks)
    dataTable <- reactive({
      DEresults()$DEtableSubset %>%
        DT::datatable() %>%
        DT::formatSignif(columns = c('pval', 'pvalAdj','lfc','log2_intensity'), digits = 3)
    })

    output[['data']] <- DT::renderDataTable(dataTable())

    #DE data download
    output[['download']] <- downloadHandler(
      filename = function() {
        paste(input[['fileName']])
      },
      content = function(file) {
        utils::write.csv(x = DEresults()$DEtableSubset, file = file, row.names = FALSE)
      }
    )

    #Output selected peaks
    selectedPeaks <- reactive({
      DEresults()$DEtableSubset$m_z[input$data_rows_selected]
    })

    proxy = DT::dataTableProxy('data')

    observe({proxy %>% DT::selectRows(NULL)}) %>%
      bindEvent(input[['resetSelection']])

    observe({proxy %>% DT::selectRows(selected = 1:50)}) %>%
      bindEvent(input[['selectTop50']])

    return(reactive(list('DE' = DEresults,
                         'selectedPeaks' = reactive(selectedPeaks())
    )))

  })
}
