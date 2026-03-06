#' Performs differential intensity analysis
#'
#' This module provides UI and server logic for performing differential intensity
#' analysis on pseudobulk data. It allows users to compare peak intensities
#' between two condition groups using either t-tests or Wilcox rank-sum tests.
#' @details
#' \itemize{
#'  \item{Interactive selection of metadata columns and condition groups to compare}
#'  \item{Choice of statistical tests (parametric and non-parametric)}
#'  \item{Adjustable p-value and log2 fold-change thresholds}
#'  \item{Interactive table with row selection for downstream analyses}
#'  \item{Export results to CSV}
#' }
#' @name BulkPanel_DATab
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame with bulk sample metadata.
#' @param show Logical; whether to render the panel (default: TRUE).
#'
#' @return A [shiny::tabPanel] containing the UI elements.
#'
#' @rdname BulkPanel_DATab
#' @export
BulkPanel_DATabUI <- function(id, bulk.metadata, show = TRUE) {
  ns <- shiny::NS(id)

  if (!show) {
    return(NULL)
  }

  shiny::tabPanel(
    'Differential Analysis',
    shiny::tags$h1("Differential intensity analysis"),
    shinyjs::useShinyjs(),
    shiny::sidebarLayout(
      # ========================================================================
      # SIDEBAR PANEL - USER INPUTS
      # ========================================================================
      shiny::sidebarPanel(
        # Information menu
        shinyWidgets::dropMenu(
          shinyWidgets::circleButton(ns("info_differential_analysis"), icon = shiny::icon("info-circle"), status = "info"),
          shiny::tags$h3("Differential Intensity Analysis"),
          shiny::tags$ul(
            shiny::tags$li("This tab compares peak intensities between two user-selected condition groups using t-test or Wilcox rank sum test."),
            shiny::tags$li("Results table shows significant peaks, using Benjamini-Hochberg multiple-testing correction and can be downloaded as a CSV."),
            shiny::tags$li("Settings allow selection of metadata, conditions, test type, and filtering thresholds."),
            shiny::tags$li("These differential analysis results are passed on to DA Summary and ORA tabs.")
          ),
          theme = "light-border",
          placement = "right",
          arrow = FALSE
        ),

        # Metadata column selection
        shiny::selectInput(
          ns('condition'), 'Metadata column to use:',
          choices = colnames(bulk.metadata[, sapply(bulk.metadata, dplyr::n_distinct) != nrow(bulk.metadata)]),
          selected = colnames(bulk.metadata[, sapply(bulk.metadata, dplyr::n_distinct) != nrow(bulk.metadata)])[1]
        ),

        # Condition group selection
        shiny::selectInput(ns('variable1'), 'Condition 1:',
                    unique(bulk.metadata[[ncol(bulk.metadata)]])),
        shiny::selectInput(ns('variable2'), 'Condition 2:',
                    unique(bulk.metadata[[ncol(bulk.metadata)]]),
                    selected = unique(bulk.metadata[[ncol(bulk.metadata)]])[2]),

        # Statistical test selection
        shiny::selectInput(ns('pipeline'), 'Statistical test:',
                    c("t-test", "Wilcox rank sum")),

        shiny::hr(),
        shiny::tags$b("Filtering Thresholds"),

        # p-value threshold
        shiny::sliderInput(ns('pvalThreshold'),
                    label = 'Adjusted p-value threshold',
                    min = 0, value = 0.05, max = 1, step = 0.005),

        # log2 fold-change threshold
        shiny::sliderInput(ns('lfcThreshold'),
                    label = 'log2 fold change threshold',
                    min = 0, value = 1, max = 5, step = 0.1),

        shiny::hr(),

        # Run analysis button
        shiny::actionButton(ns('goDE'),
                     label = 'Start differential intensity analysis',
                     class = "btn-primary"),

        shiny::hr(),
        shiny::tags$b("Download Results"),
        shinyWidgets::dropMenu(
          shinyWidgets::circleButton(ns("downloads_table"), icon = shiny::icon("download"), status = "success"),
          shiny::tags$div(
            shiny::tags$h3("Download Results Table"),
            shiny::fluidRow(
              shiny::column(12,
                shiny::textInput(ns('fileName'), 'File name for download', value = 'DIAset.csv', placeholder = 'DIAset.csv'),
                shiny::downloadButton(ns('download'), 'Download Table')
              )
            )
          ),
          theme = "light-border",
          placement = "right",
          arrow = FALSE
        ),

        shiny::hr(),
        shiny::tags$b("Peak Selection"),

        # Peak selection controls
        shiny::div("Select peaks by clicking rows in the table below."),
        shiny::div(style = "margin-bottom: 10px"),
        shiny::actionButton(ns('resetSelection'), label = "Reset selection"),
        shiny::div(style = "margin-bottom: 10px"),
        shiny::actionButton(ns('selectTop50'), label = "Select top 50 peaks")
      ),

      # ========================================================================
      # MAIN PANEL - RESULTS TABLE
      # ========================================================================
      shiny::mainPanel(
        DT::DTOutput(ns('data'))
      )
    )
  )
}

#' @rdname BulkPanel_DATab
#' @param bulk.intensity.matrix Reactive expression returning numeric matrix of
#'   intensity (peaks × samples).
#' @param anno Data frame with peak annotation.
#'
#' @return Reactive list containing:
#'   - `DE`: Function returning list of differential intensity results
#'   - `selectedPeaks`: Reactive vector of selected peak m/z values
#'   - `runDE`: Numeric flag indicating when DE button was pressed
#' @export
BulkPanel_DATabServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno) {

  shiny::moduleServer(id, function(input, output, session) {
    # Track when DE button is pressed
    r <- shiny::reactiveValues(go_DE = 0)
    shiny::observeEvent(input$goDE, {r$go_DE <- 1})

    # Update condition groups when metadata column changes
    shiny::observe({
      condition_choices <- unique(bulk.metadata[[input[["condition"]]]])
      shiny::updateSelectInput(session, 'variable1', choices = condition_choices)
      shiny::updateSelectInput(session, 'variable2', choices = condition_choices,
                        selected = condition_choices[2])
    })

    # Update available statistical tests
    shiny::observe({
      condition_indices <- bulk.metadata[[input[["condition"]]]] %in%
        c(input[['variable1']], input[['variable2']])
      choices <- c("t-test", "Wilcox rank sum")
      shiny::updateSelectInput(session, 'pipeline', choices = choices)
    })

    # ========================================================================
    # REACTIVE: DIFFERENTIAL EXPRESSION COMPUTATION
    # ========================================================================
    DEresults <- shiny::reactive({
      shinyjs::disable("goDE")
      condition_indices <- bulk.metadata[[input[["condition"]]]] %in%
        c(input[['variable1']], input[['variable2']])

      # Perform differential analysis
      DEtable <- bulk_utils_DA(
        intensity_matrix = bulk.intensity.matrix[, condition_indices],
        condition = bulk.metadata[[input[["condition"]]]][condition_indices],
        var1 = input[['variable1']],
        var2 = input[['variable2']],
        test = input[["pipeline"]],
        anno = anno
      )

      # Apply thresholds
      DEtableSubset <- DEtable |>
        dplyr::filter(.data$pvalAdj < input[["pvalThreshold"]] &
                        abs(.data$lfc) > input[['lfcThreshold']]) |>
        dplyr::arrange(.data$pvalAdj)

      shinyjs::enable("goDE")

      return(list(
        'DEtable' = DEtable,
        "DEtableSubset" = DEtableSubset,
        'pvalThreshold' = input[["pvalThreshold"]],
        'lfcThreshold' = input[['lfcThreshold']]
      ))
    }) |>
      shiny::bindCache(utils::head(bulk.intensity.matrix), bulk.metadata,
                input[["condition"]], input[['variable1']], input[['variable2']],
                input[["pipeline"]], input[["pvalThreshold"]],
                input[['lfcThreshold']]) |>
      shiny::bindEvent(input[["goDE"]])

    # ========================================================================
    # OUTPUT: RESULTS TABLE
    # ========================================================================
    dataTable <- shiny::reactive({
      DEresults()$DEtableSubset |>
        DT::datatable(selection = 'multiple') |>
        DT::formatSignif(columns = c('pval', 'pvalAdj', 'lfc', 'log2_intensity'),
                         digits = 3)
    })

    output[['data']] <- DT::renderDataTable(dataTable())

    # ========================================================================
    # OUTPUT: DATA DOWNLOAD
    # ========================================================================
    output[['download']] <- shiny::downloadHandler(
      filename = function() {
        paste(input[['fileName']])
      },
      content = function(file) {
        utils::write.csv(x = DEresults()$DEtableSubset, file = file,
                         row.names = FALSE)
      }
    )

    # ========================================================================
    # REACTIVE: SELECTED PEAKS
    # ========================================================================
    selectedPeaks <- shiny::reactive({
      DEresults()$DEtableSubset$m_z[input$data_rows_selected]
    })

    # Table proxy for row selection management
    proxy <- DT::dataTableProxy('data')

    # Reset row selection
    shiny::observe({
      proxy |> DT::selectRows(NULL)
    }) |>
      shiny::bindEvent(input[['resetSelection']])

    # Select top 50 peaks
    shiny::observe({
      proxy |> DT::selectRows(selected = 1:min(50,
                                                 nrow(DEresults()$DEtableSubset)))
    }) |>
      shiny::bindEvent(input[['selectTop50']])

    # Return results for use by other modules
    return(shiny::reactive(list(
      'DE' = DEresults,
      'selectedPeaks' = shiny::reactive(selectedPeaks()),
      'runDE' = r$go_DE
    )))
  })
}
