#' Infers and visualises region-based covariation networks
#'
#' @description UI and server logic for the Region-Based Covariation Network Inference panel, enabling users to infer, visualise, and compare covariation networks across user-defined regions or clusters in spatial omics data. Supports network inference, visualisation, and download of results for selected regions.
#'
#' @details
#' \itemize{
#'   \item{Infer covariation networks for selected regions or clusters and compare across multiple regions.}
#'   \item{Visualise the resulting network structure and the similarity between networks based on different regions.}
#'   \item{Target metabolite peaks are used to build the network with a user-specified number of edges.}
#'   \item{All network plots and results are available for download.}
#' }
#'
#' @param id Shiny module id (for both UI and server)
#' @param metadata Data frame of sample or region metadata
#' @param show Logical; whether to show the panel (default TRUE)
#' @param intensity.matrix Matrix of intensities (features x samples)
#' @param anno Data frame of peak annotations
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param shared_data Reactive or shared data object
#' @name RegionPanel_GRNTab
#' @rdname RegionPanel_GRNTab
#' @return UI: A shiny tabPanel object for the GRN tab. Server: None (side effects in Shiny module).
#' @export
RegionPanel_GRNTabUI <- function(id, metadata, show = TRUE){
  ns <- shiny::NS(id)
  if(show){
    shiny::tabPanel(
      'Region-Based Covariation Network Inference',
          bslib::accordion(
            bslib::accordion_panel(
              title = "Information",
              icon = bsicons::bs_icon("arrow-right-circle"),
              shiny::tags$ul(
                shiny::tags$li("Infer covariation networks for selected regions or clusters and compare across multiple regions."),
                shiny::tags$li("Visualize the resulting network structure and the similarity between networks based on different regions."),
                shiny::tags$li("Target metabolite peaks are used to build the network with a user-specific number of edges.")
              )
            ),
            bslib::accordion_panel(
              title = "Downloads",
              icon = bsicons::bs_icon("download"),
              shiny::tags$div(
                shiny::tags$h4("Download output plots"),
                shiny::fluidRow(
                  shiny::column(6,
                    shiny::tags$strong("GRN plot"),
                    shiny::textInput(ns('plotFileName'), 'File name for plot download', value ='GRNplot.html')),
                    shiny::column(6,
                    shiny::selectInput(ns('plotId'), 'Select plot to download:', 1:4),
                    shiny::downloadButton(ns('download'), 'Download Plot')
                  )
                )
              )
            ),
            id = ns("acc"),
            open = FALSE
          ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput(ns('n_networks'), 'Number of networks:', 1:4),
          shiny::selectInput(ns('condition'), 'Metadata column to use:', choices = NULL),
          shiny::selectInput(ns('samples1'), 'Samples for GRN #1:', unique(metadata[[ncol(metadata)]]),
            selected = unique(metadata[[ncol(metadata)]]), multiple = TRUE),
          shiny::conditionalPanel(
            condition = "input.n_networks >= 2",
            ns=ns,
            shiny::selectInput(ns('samples2'), 'Samples for GRN #2:', unique(metadata[[ncol(metadata)]]),
              selected = unique(metadata[[ncol(metadata)]]), multiple = TRUE)
          ),
          shiny::conditionalPanel(
            condition = "input.n_networks >= 3",
            ns=ns,
            shiny::selectInput(ns('samples3'), 'Samples for GRN #3:', unique(metadata[[ncol(metadata)]]),
              selected = unique(metadata[[ncol(metadata)]]), multiple = TRUE)
          ),
          shiny::conditionalPanel(
            condition = "input.n_networks >= 4",
            ns=ns,
            shiny::selectInput(ns('samples4'), 'Samples for GRN #4:', unique(metadata[[ncol(metadata)]]),
              selected = unique(metadata[[ncol(metadata)]]), multiple = TRUE)
          ),
          shiny::selectizeInput(ns("targets"), "Target genes:", multiple = TRUE, choices = NULL),
          shinyjs::disabled(shiny::actionButton(ns('goGRN'), label = 'Start GRN inference')),
          shiny::numericInput(ns("plotConnections"), "Connections to plot:", 5, 0, 100)
        ),
        shiny::mainPanel(
          shiny::fluidRow(
            shiny::column(6, visNetwork::visNetworkOutput(ns('plot1'))),
            shiny::column(
              6, 
              shiny::conditionalPanel(
                condition = "input.n_networks >= 2",
                ns = ns,
                visNetwork::visNetworkOutput(ns('plot2'))
              )
            )
          ),
          shiny::conditionalPanel(
            condition = "input.n_networks >= 3",
            ns = ns,
            shiny::fluidRow(
              shiny::column(6, visNetwork::visNetworkOutput(ns('plot3'))),
              shiny::conditionalPanel(
                condition = "input.n_networks > 3",
                ns = ns,
                shiny::column(6, visNetwork::visNetworkOutput(ns('plot4')))
              )
            )
          ),
          shiny::conditionalPanel(
            condition = "input.n_networks > 1",
            ns = ns,
            shiny::plotOutput(ns('plotUpset'))
          )
        )
      )
    )
  } else {
    NULL
  }
}

#' @rdname RegionPanel_GRNTab
#' @export
RegionPanel_GRNTabServer <- function(id, intensity.matrix, metadata, anno, bulk.metadata, shared_data){
  
  shiny::moduleServer(id, function(input, output, session){
    
    shiny::observe({
      categorical_cols <- colnames(shared_data$updated.metadata)[sapply(shared_data$updated.metadata, function(col) {
        n_unique <- length(unique(stats::na.omit(col)))
        n_unique >= 1 && n_unique < 20  # exclude constants & numeric-like columns
      })]
      sub_sample_cols = categorical_cols[!(categorical_cols %in% colnames(bulk.metadata))]
      shiny::updateSelectInput(session, "condition", choices = sub_sample_cols)
      
      shiny::updateSelectizeInput(
        session, "targets", server = TRUE,
        choices = anno$display_name,
        selected = utils::head(anno$display_name)
      )
    })
    
    shiny::observe({
      shiny::req(input$condition)
      shiny::updateSelectInput(session, 'samples1', 
                        choices = unique(shared_data$updated.metadata[[input[["condition"]]]]), 
                        selected = unique(shared_data$updated.metadata[[input[["condition"]]]]))
      shiny::updateSelectInput(session, 'samples2', 
                        choices = unique(shared_data$updated.metadata[[input[["condition"]]]]), 
                        selected = unique(shared_data$updated.metadata[[input[["condition"]]]]))
      shiny::updateSelectInput(session, 'samples3', 
                        choices = unique(shared_data$updated.metadata[[input[["condition"]]]]), 
                        selected = unique(shared_data$updated.metadata[[input[["condition"]]]]))
      shiny::updateSelectInput(session, 'samples4', 
                        choices = unique(shared_data$updated.metadata[[input[["condition"]]]]), 
                        selected = unique(shared_data$updated.metadata[[input[["condition"]]]]))
      
    })
    
    intensity.matrix.sub <- shiny::reactive(# Remove genes of constant expression
      t(intensity.matrix[, matrixStats::colMins(intensity.matrix) !=
                            matrixStats::colMaxs(intensity.matrix) ])
    )

    shiny::observe({
      shiny::req(input$condition)
      shiny::req(input$samples1)
      enable_condition <- length(input[["targets"]]) >= 1 &
        (input[["n_networks"]] < 1 | length(input[["samples1"]]) > 0) &
        (input[["n_networks"]] < 2 | length(input[["samples2"]]) > 0) &
        (input[["n_networks"]] < 3 | length(input[["samples3"]]) > 0) &
        (input[["n_networks"]] < 4 | length(input[["samples4"]]) > 0) &
        (input[["n_networks"]] < 1 | (sum((shared_data$updated.metadata[,input[["condition"]]]%in%input[["samples1"]])) > 1)) &
        (input[["n_networks"]] < 2 | (sum((shared_data$updated.metadata[,input[["condition"]]]%in%input[["samples2"]])) > 1)) &
        (input[["n_networks"]] < 3 | (sum((shared_data$updated.metadata[,input[["condition"]]]%in%input[["samples3"]])) > 1)) &
        (input[["n_networks"]] < 4 | (sum((shared_data$updated.metadata[,input[["condition"]]]%in%input[["samples4"]])) > 1))
      if(enable_condition){
        shinyjs::enable("goGRN")
      }else{
        shinyjs::disable("goGRN")
      }
    }) |>
      shiny::bindEvent(input[["condition"]], input[["targets"]], input[["n_networks"]], input[["samples1"]],
                input[["samples2"]], input[["samples3"]], input[["samples4"]])

    n_networks <- shiny::reactive(input[["n_networks"]]) |> shiny::bindEvent(input[["goGRN"]])
    shiny::observe(shiny::updateSelectInput(session, "plotId", choices = seq_len(n_networks())))
    GRNresults1 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 1){
        weightMat <- NULL
        shiny::withProgress(message = "Inferring GRN for network 1...", value = 0, {
          weightMat <- bulk_utils_infer_GRN(
            intensity_matrix = intensity.matrix.sub(),
            metadata = shared_data$updated.metadata,
            anno = anno,
            targets = input[["targets"]],
            condition = input[["condition"]],
            samples = input[["samples1"]],
            inference_method = "GENIE3"
          )
        })
        shinyjs::enable("goGRN")
      }else{
        weightMat <- NULL
      }
      weightMat
    }) |>
      shiny::bindCache(utils::head(intensity.matrix), shared_data$updated.metadata, input[["condition"]],
                input[['samples1']],input[['targets']]) |>
      shiny::bindEvent(input[["goGRN"]])
    GRNresults2 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 2){
        weightMat <- NULL
        shiny::withProgress(message = "Inferring GRN for network 2...", value = 0, {
          weightMat <- bulk_utils_infer_GRN(
            intensity_matrix = intensity.matrix.sub(),
            metadata = shared_data$updated.metadata,
            anno = anno,
            targets = input[["targets"]],
            condition = input[["condition"]],
            samples = input[["samples2"]],
            inference_method = "GENIE3"
          )
        })
        shinyjs::enable("goGRN")
      }else{
        weightMat <- NULL
      }
      weightMat
    }) |>
      shiny::bindCache(utils::head(intensity.matrix), shared_data$updated.metadata, input[["condition"]],
                input[['samples2']],input[['targets']]) |>
      shiny::bindEvent(input[["goGRN"]])
    GRNresults3 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 3){
        weightMat <- NULL
        shiny::withProgress(message = "Inferring GRN for network 3...", value = 0, {
          weightMat <- bulk_utils_infer_GRN(
            intensity_matrix = intensity.matrix.sub(),
            metadata = shared_data$updated.metadata,
            anno = anno,
            targets = input[["targets"]],
            condition = input[["condition"]],
            samples = input[["samples3"]],
            inference_method = "GENIE3"
          )
        })
        shinyjs::enable("goGRN")
      }else{
        weightMat <- NULL
      }
      weightMat
    }) |>
      shiny::bindCache(utils::head(intensity.matrix), shared_data$updated.metadata, input[["condition"]],
                input[['samples3']],input[['targets']]) |>
      shiny::bindEvent(input[["goGRN"]])
    GRNresults4 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 4){
        weightMat <- NULL
        shiny::withProgress(message = "Inferring GRN for network 4...", value = 0, {
          weightMat <- bulk_utils_infer_GRN(
            intensity_matrix = intensity.matrix.sub(),
            metadata = shared_data$updated.metadata,
            anno = anno,
            targets = input[["targets"]],
            condition = input[["condition"]],
            samples = input[["samples4"]],
            inference_method = "GENIE3"
          )
        })
        shinyjs::enable("goGRN")
      }else{
        weightMat <- NULL
      }
      weightMat
    }) |>
      shiny::bindCache(utils::head(intensity.matrix), shared_data$updated.metadata, input[["condition"]],
                input[['samples4']],input[['targets']]) |>
      shiny::bindEvent(input[["goGRN"]])

    weightMatList <- shiny::reactive({
      weightMatList <- list()
      if(n_networks() >= 1) {
        weightMatList[[1]] <- GRNresults1()
        if(n_networks() >= 2) {
          weightMatList[[2]] <- GRNresults2()
          if(n_networks() >= 3) {
            weightMatList[[3]] <- GRNresults3()
            if(n_networks() >= 4) {
              weightMatList[[4]] <- GRNresults4()
            }
          }
        }
      }
      weightMatList
    })

    recurring_regulators <- shiny::reactive({
      bulk_utils_find_regulators_with_recurring_edges(weightMatList(), input[["plotConnections"]])
    })

    GRNplot1 <- shiny::reactive(bulk_utils_plot_GRN(
      weightMat = GRNresults1(),
      anno = anno,
      plotConnections = input[["plotConnections"]],
      plot_position_grid = 1,
      n_networks = n_networks(),
      recurring_regulators = recurring_regulators()
    ))
    GRNplot2 <- shiny::reactive(bulk_utils_plot_GRN(
      weightMat = GRNresults2(),
      anno = anno,
      plotConnections = input[["plotConnections"]],
      plot_position_grid = 2,
      n_networks = n_networks(),
      recurring_regulators = recurring_regulators()
    ))
    GRNplot3 <- shiny::reactive(bulk_utils_plot_GRN(
      weightMat = GRNresults3(),
      anno = anno,
      plotConnections = input[["plotConnections"]],
      plot_position_grid = 3,
      n_networks = n_networks(),
      recurring_regulators = recurring_regulators()
    ))
    GRNplot4 <- shiny::reactive(bulk_utils_plot_GRN(
      weightMat = GRNresults4(),
      anno = anno,
      plotConnections = input[["plotConnections"]],
      plot_position_grid = 4,
      n_networks = n_networks(),
      recurring_regulators = recurring_regulators()
    ))

    upsetPlot <- shiny::reactive(bulk_utils_plot_upset(weightMatList(), input[["plotConnections"]]))

    output[['plot1']] <- visNetwork::renderVisNetwork(GRNplot1())
    output[['plot2']] <- visNetwork::renderVisNetwork(GRNplot2())
    output[['plot3']] <- visNetwork::renderVisNetwork(GRNplot3())
    output[['plot4']] <- visNetwork::renderVisNetwork(GRNplot4())
    output[['plotUpset']] <- shiny::renderPlot(upsetPlot())

    output[['download']] <- shiny::downloadHandler(
      filename = function() {input[['plotFileName']]},
      content = function(file) {
        GRNplot <- get(paste0("GRNplot", input[["plotId"]]))()
        GRNplot |> visNetwork::visSave(file)
      }
    )

  })
}