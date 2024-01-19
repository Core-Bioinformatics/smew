# need to check that features aren't getting combined!!

#' @rdname BulkDEPanel
#' @export
BulkDEpanelUI <- function(id, bulk.metadata, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Differential intensity analysis',
      tags$h1("Differential intensity analysis"),
      shinyjs::useShinyjs(),
      sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(
          dropMenu(
            circleButton(ns("info_differential_analysis"), icon = icon("info"),status = "success"),
            tags$div(
              tags$h3("Differential intensity analysis"),
              tags$ul(
                tags$li("Select a sample-wide metadata column and 2 condition groups to compare using t-tests or Wilcox rank sum tests"),
                tags$li("log2FC and BH-adjusted p-values can be capped using the slides."),
                tags$li("Comparisons are only performed once the button has been pressed."),
                tags$li("The most recent comparison is passed onto other tabs for visualisation, pathway analysis etc."),
                tags$li("The table of peaks showing significant changes can also be downloaded as a csv."),
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),

          selectInput(ns('condition'), 'Metadata column to use:', colnames(bulk.metadata)[-1],
                      selected = colnames(bulk.metadata)[ncol(bulk.metadata)]),

          # Input: Selector variables to compare
          selectInput(ns('variable1'), 'Condition 1:', unique(bulk.metadata[[ncol(bulk.metadata)]])),
          selectInput(ns('variable2'), 'Condition 2:', unique(bulk.metadata[[ncol(bulk.metadata)]]),
                      selected = unique(bulk.metadata[[ncol(bulk.metadata)]])[2]),

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
          DT::DTOutput(ns('data'))
        )
      )
    )
  }else{
    NULL
  }
}

#' @rdname BulkDEPanel
#' @export
BulkDEpanelServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(bulk.intensity.matrix)
    is.reactive(bulk.metadata)
    !is.reactive(anno)
  })

  moduleServer(id, function(input, output, session){
    r <- reactiveValues(go_DE=0)
    observeEvent(input$goDE,{r$go_DE <- 1})

    observe({
      updateSelectInput(session, 'variable1', choices = unique(bulk.metadata[[input[["condition"]]]]))
      updateSelectInput(session, 'variable2', choices = unique(bulk.metadata[[input[["condition"]]]]),
                        selected = unique(bulk.metadata[[input[["condition"]]]])[2])
    })

    observe({
      condition.indices <- bulk.metadata[[input[["condition"]]]] %in% c(input[['variable1']], input[['variable2']])
      choices <- c("t-test", "Wilcox rank sum")
      updateSelectInput(session, 'pipeline', choices = choices)
    })

    DEresults <- reactive({
      shinyjs::disable("goDE")
      condition.indices <- bulk.metadata[[input[["condition"]]]] %in% c(input[['variable1']], input[['variable2']])
      # Need to add error if any group has <2 samples
      DEtable <- DEanalysis(
        intensity.matrix = bulk.intensity.matrix[, condition.indices],
        condition = bulk.metadata[[input[["condition"]]]][condition.indices],
        var1 = input[['variable1']],
        var2 = input[['variable2']],
        test = input[["pipeline"]],
        anno = anno
      )
      DEtableSubset <- DEtable %>%
        dplyr::filter(.data$pvalAdj < input[["pvalThreshold"]] & abs(.data$lfc) > input[['lfcThreshold']]) |>
        dplyr::arrange(.data$pvalAdj)

      #the thresholds are returned here so that MA/volcano and table display
      #don't use new thresholds without the button being used
      shinyjs::enable("goDE")
      return(list('DEtable' = DEtable,
                  "DEtableSubset" = DEtableSubset,
                  'pvalThreshold' = input[["pvalThreshold"]],
                  'lfcThreshold' = input[['lfcThreshold']]))
    }) %>%
      bindCache(utils::head(bulk.intensity.matrix), bulk.metadata, input[["condition"]],
                input[['variable1']], input[['variable2']], input[["pipeline"]],
                input[["pvalThreshold"]],input[['lfcThreshold']]) |>
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
                         'selectedPeaks' = reactive(selectedPeaks()),
                         'runDE'=r$go_DE
    )))


  })
}
