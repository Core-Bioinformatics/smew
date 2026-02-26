#' Metabolite Peak Regulatory Network (GRN) UI and Server Module
#'
#' This file contains the UI and server logic for the GRN panel. Helper functions for GRN inference and network processing are now located in bulkPanel_GRNFuns.R.
#'
#' @keywords internal
#' @name GRNpanel
NULL

#' @rdname GRNpanel
#' @export
BulkGRNUI <- function(id, bulk.metadata, show = TRUE){
  ns <- shiny::NS(id)
  
  if(show){
    shiny::tabPanel(
      'Covariation network inference',
      shiny::sidebarLayout(
        shiny::sidebarPanel(
            shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-top: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_grn"), icon = shiny::icon("info-circle"), status = "info"),
              shiny::tags$div(
                shiny::tags$h3("GRN Inference Information"),
                shiny::tags$ul(
                  shiny::tags$li("This panel allows inference of metabolite peak covariation networks from pseudobulked sample-level data for up to 4 sample groups."),
                  shiny::tags$li("Select samples for each network, choose target metabolite peaks, and set the number of connections to plot."),
                  shiny::tags$li("After running inference, interactive network plots are shown for each group (in the order top left, top right, bottom left, bottom right)."),
                  shiny::tags$li("Nodes are coloured as follows: blue for selected targets, grey for predicted regulators, green for overlapping regulators between networks."),
                  shiny::tags$li("Top connections and edge weights are determined by the power of each peak to predict the intensity of the target peaks using GENIE3's random forest approach."),
                  shiny::tags$li("The UpSet plot visualizes overlap of edges between networks."),
                  shiny::tags$li("All plots are downloadable via the download button."),
                  shiny::tags$li("Settings for plot downloads are available in the download menu.")
                )
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_grn"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Downloads"),
                shiny::tags$h4("GRN network plot"),
                shiny::textInput(ns('plotFileName'), 'File name for plot download', value ='GRNplot.html'),
                shiny::selectInput(ns('plotId'), 'Select plot to download:', 1:4),
                shiny::numericInput(ns('plotWidth'), value = 8, label = 'Width (in)', min = 1, max = 50, step = 1),
                shiny::numericInput(ns('plotHeight'), value = 6, label = 'Height (in)', min = 1, max = 50, step = 1),
                shiny::downloadButton(ns('download'), 'Download Plot')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::selectInput(ns('n_networks'), 'Number of networks:', 1:4),
          shiny::selectInput(ns('samples1'), 'Samples for GRN #1:', unique(bulk.metadata$Sample),
                      selected = unique(bulk.metadata$Sample), multiple = TRUE),
          shiny::conditionalPanel(
            id = ns('samples2'),
            ns=ns,
            condition = "input.n_networks >= 2",
            shiny::selectInput(ns('samples2'), 'Samples for GRN #2:', unique(bulk.metadata$Sample),
                        selected = unique(bulk.metadata$Sample), multiple = TRUE),
          ),
          shiny::conditionalPanel(
            id = ns('samples3'),
            ns=ns,
            condition = "input.n_networks >= 3",
            shiny::selectInput(ns('samples3'), 'Samples for GRN #3:', unique(bulk.metadata$Sample),
                        selected = unique(bulk.metadata$Sample), multiple = TRUE),
          ),
          shiny::conditionalPanel(
            id = ns('samples4'),
            ns=ns,
            condition = "input.n_networks >= 4",
            shiny::selectInput(ns('samples4'), 'Samples for GRN #4:', unique(bulk.metadata$Sample),
                        selected = unique(bulk.metadata$Sample), multiple = TRUE),
          ),
          shiny::selectizeInput(ns("targets"), "Target metabolite peaks:", multiple = TRUE, choices = NULL),
          shinyjs::disabled(shiny::actionButton(ns('goGRN'), label = 'Start GRN inference')),
          shiny::numericInput(ns("plotConnections"), "Connections to plot:", 5, 0, 100),
        ),
        
        shiny::mainPanel(
          shiny::fluidRow(
            shiny::column(6, visNetwork::visNetworkOutput(ns('plot1'))),
            shiny::column(
              6, 
              shiny::conditionalPanel(
                id = ns('plot2col'),
                ns = ns,
                condition = "input.n_networks >= 2",
                visNetwork::visNetworkOutput(ns('plot2'))
              )
            )
          ),
          shiny::conditionalPanel(
            id = ns('plotrow'),
            ns = ns,
            condition = "input.n_networks >= 3",
            shiny::fluidRow(
              shiny::column(6, visNetwork::visNetworkOutput(ns('plot3'))),
              shiny::conditionalPanel(
                id = ns('plot2col2'),
                ns = ns,
                condition = "input.n_networks > 3",
                
                shiny::column(
                  6,
                  visNetwork::visNetworkOutput(ns('plot4'))
                )
              )
            ),
            shiny::conditionalPanel(
              id = ns('includeUpset'),
              ns = ns,
              condition = "input.n_networks > 1",
              shiny::plotOutput(ns('plotUpset'))
            )
          )
        )
      )
    )
    
  }else{
    NULL
  }
}

#' @rdname GRNpanel
#' @export
BulkGRNServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno){
  
  # stopifnot({
  #   is.reactive(intensity.matrix)
  #   is.reactive(metadata)
  #   !is.reactive(anno)
  # })
  # 
  shiny::moduleServer(id, function(input, output, session){
    
    shiny::observe({
      shiny::updateSelectizeInput(
        session, "targets", server = TRUE,
        choices = anno$display_name,
        selected = utils::head(anno$display_name)
      )
      
      shiny::updateSelectInput(session, 'samples1', 
                        choices = unique(bulk.metadata[["Sample"]]), 
                        selected = unique(bulk.metadata[["Sample"]]))
      shiny::updateSelectInput(session, 'samples2', 
                        choices = unique(bulk.metadata[["Sample"]]), 
                        selected = unique(bulk.metadata[["Sample"]]))
      shiny::updateSelectInput(session, 'samples3', 
                        choices = unique(bulk.metadata[["Sample"]]),
                        selected = unique(bulk.metadata[["Sample"]]))
      shiny::updateSelectInput(session, 'samples4',
                        choices = unique(bulk.metadata[["Sample"]]),
                        selected = unique(bulk.metadata[["Sample"]]))

      bulk.intensity.matrix = as.matrix(bulk.intensity.matrix)

    })
    
    bulk.intensity.matrix.sub <- shiny::reactive(
      # Remove peaks of constant intensity
      as.matrix(bulk.intensity.matrix[matrixStats::rowMins(as.matrix(bulk.intensity.matrix)) !=
                         matrixStats::rowMaxs(as.matrix(bulk.intensity.matrix)), ])
    )
    
    shiny::observe({
      shiny::req(input$samples1)
      enable_condition <- length(input[["targets"]]) >= 1 &
        (input[["n_networks"]] < 1 | length(input[["samples1"]]) > 0) &
        (input[["n_networks"]] < 2 | length(input[["samples2"]]) > 0) &
        (input[["n_networks"]] < 3 | length(input[["samples3"]]) > 0) &
        (input[["n_networks"]] < 4 | length(input[["samples4"]]) > 0) &
        (input[["n_networks"]] < 1 | (sum((bulk.metadata[,"Sample"]%in%input[["samples1"]])) > 1)) &
        (input[["n_networks"]] < 2 | (sum((bulk.metadata[,"Sample"]%in%input[["samples2"]])) > 1)) &
        (input[["n_networks"]] < 3 | (sum((bulk.metadata[,"Sample"]%in%input[["samples3"]])) > 1)) &
        (input[["n_networks"]] < 4 | (sum((bulk.metadata[,"Sample"]%in%input[["samples4"]])) > 1))
      if(enable_condition){
        shinyjs::enable("goGRN")
      }else{
        shinyjs::disable("goGRN")
      }
    }) |> 
      shiny::bindEvent(input[["targets"]], input[["n_networks"]], input[["samples1"]],
                input[["samples2"]], input[["samples3"]], input[["samples4"]])
    
    n_networks <- shiny::reactive(input[["n_networks"]]) |> shiny::bindEvent(input[["goGRN"]])
    shiny::observe(shiny::updateSelectInput(session, "plotId", choices = seq_len(n_networks())))
    
    GRNresults1 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 1){
        weightMat <- infer_GRN(
          intensity_matrix = bulk.intensity.matrix.sub(),
          metadata = bulk.metadata,
          anno = anno,
          targets = anno$m_z[match(input[["targets"]], anno$display_name)],
          condition = "Sample",
          samples = input[["samples1"]],
          inference_method = "GENIE3"
        )
        shinyjs::enable("goGRN")
      }else{
        weightMat <- NULL
      }
      weightMat
    }) |> 
      shiny::bindCache(utils::head(bulk.intensity.matrix), bulk.metadata, "Sample",
                input[['samples1']],input[['targets']]) |> 
      shiny::bindEvent(input[["goGRN"]])
    GRNresults2 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 2){
        weightMat <- infer_GRN(
          intensity_matrix = bulk.intensity.matrix.sub(),
          metadata = bulk.metadata,
          anno = anno,
          targets = anno$m_z[match(input[["targets"]], anno$display_name)],
          condition = "Sample",
          samples = input[["samples2"]],
          inference_method = "GENIE3"
        )
        shinyjs::enable("goGRN")
      }else{
        weightMat <- NULL
      }
      weightMat
    }) |> 
      shiny::bindCache(utils::head(bulk.intensity.matrix), bulk.metadata, "Sample",
                input[['samples2']],input[['targets']]) |> 
      shiny::bindEvent(input[["goGRN"]])
    GRNresults3 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 3){
        weightMat <- infer_GRN(
          intensity_matrix = bulk.intensity.matrix.sub(),
          metadata = bulk.metadata,
          anno = anno,
          targets = anno$m_z[match(input[["targets"]], anno$display_name)],
          condition = "Sample",
          samples = input[["samples3"]],
          inference_method = "GENIE3"
        )
        shinyjs::enable("goGRN")
      }else{
        weightMat <- NULL
      }
      weightMat
    }) |> 
      shiny::bindCache(utils::head(bulk.intensity.matrix), bulk.metadata, "Sample",
                input[['samples3']],input[['targets']]) |> 
      shiny::bindEvent(input[["goGRN"]])
    GRNresults4 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 4){
        weightMat <- infer_GRN(
          intensity_matrix = bulk.intensity.matrix.sub(),
          metadata = bulk.metadata,
          anno = anno,
          targets = anno$m_z[match(input[["targets"]], anno$display_name)],
          condition = "Sample",
          samples = input[["samples4"]],
          inference_method = "GENIE3"
        )
        shinyjs::enable("goGRN")
      }else{
        weightMat <- NULL
      }
      weightMat
    }) |> 
      shiny::bindCache(utils::head(bulk.intensity.matrix), bulk.metadata, "Sample",
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
      find_regulators_with_recurring_edges(weightMatList(), input[["plotConnections"]])
    })
    
    GRNplot1 <- shiny::reactive(plot_GRN(
      weightMat = GRNresults1(),
      anno = anno,
      plotConnections = input[["plotConnections"]],
      plot_position_grid = 1,
      n_networks = n_networks(),
      recurring_regulators = recurring_regulators()
    ))
    GRNplot2 <- shiny::reactive(plot_GRN(
      weightMat = GRNresults2(),
      anno = anno,
      plotConnections = input[["plotConnections"]],
      plot_position_grid = 2,
      n_networks = n_networks(),
      recurring_regulators = recurring_regulators()
    ))
    GRNplot3 <- shiny::reactive(plot_GRN(
      weightMat = GRNresults3(),
      anno = anno,
      plotConnections = input[["plotConnections"]],
      plot_position_grid = 3,
      n_networks = n_networks(),
      recurring_regulators = recurring_regulators()
    ))
    GRNplot4 <- shiny::reactive(plot_GRN(
      weightMat = GRNresults4(),
      anno = anno,
      plotConnections = input[["plotConnections"]],
      plot_position_grid = 4,
      n_networks = n_networks(),
      recurring_regulators = recurring_regulators()
    ))
    
    upsetPlot <- shiny::reactive(plot_upset(weightMatList(), input[["plotConnections"]]))
    
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