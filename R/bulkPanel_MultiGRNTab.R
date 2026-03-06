#' Infers multi-omics covariation networks
#'
#' @description UI and server logic for comparative GRN inference across multiple sample groups using combined bulk and multi-modal data. Enables side-by-side comparison of metabolite peak regulatory relationships across different conditions or sample sets.
#'
#' @details
#' \itemize{
#'   \item{Multi-network inference (1-4 networks) from independent sample subsets}
#'   \item{GENIE3-based GRN construction from combined bulk/multi-modal data}
#'   \item{Network visualization with recurring regulators highlighted}
#'   \item{Interactive upset plots showing regulator overlap across networks}
#'   \item{Customizable connection filtering and network topology visualization}
#'   \item{Accepts bulk intensity matrix and multi-modal data, combining them for joint GRN inference. Removes constant features before analysis. Supports flexible sample selection per network and identifies regulatory patterns that are conserved or divergent across conditions.}
#' }
#' @name BulkPanel_MultiGRNTab
#' @rdname BulkPanel_MultiGRNTab
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame with bulk sample metadata
#' @param show Logical; whether to render the panel (default: TRUE)
#' @return A shiny::tabPanel containing the UI elements for the multi-modal GRN panel
#' @export
BulkPanel_MultiGRNTabUI <- function(id, bulk.metadata, show = TRUE){
  ns <- shiny::NS(id)
  
  if(show){
    shiny::tabPanel(
      'Multi-Modal Covariation Network Inference',
      shiny::sidebarLayout(
        shiny::sidebarPanel(
            shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-top: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_grn"), icon = shiny::icon("info-circle"), status = "info"),
              shiny::tags$div(
                shiny::tags$h3("Multi-GRN Inference Information"),
                shiny::tags$ul(
                  shiny::tags$li("This panel allows inference of combined multi-modal covariation network for up to 4 sample groups."),
                  shiny::tags$li("Using multi-modal data provided in app generation at the sample level, the variation of entries across multiple modalities can be compared and visualised."),
                  shiny::tags$li("Select samples for each network, choose target metabolite peaks or entries from the other modality provided, and set the number of connections to plot."),
                  shiny::tags$li("After running inference, interactive network plots are shown for each group (in the order top left, top right, bottom left, bottom right)."),
                  shiny::tags$li("Nodes are coloured as follows: blue for selected metabolite peak targets, grey for predicted metabolite peak regulators, green for overlapping metabolite peak regulators between networks, 
                          red for selected other modality targets, light pink for predicted other modality regulators and dark purple for overlapping other modality regulators between networks."),
                  shiny::tags$li("Top connections and edge weights are determined by the power of each peak/other modality entry to predict the intensity of the targets using GENIE3's random forest approach."),
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
                shiny::textInput(ns('plotFileName'), 'File name for download', value ='GRNplot.html', placeholder = 'GRNplot.html'),
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
          shiny::selectInput(ns('samples1'), 'Samples for GRN #1:', unique(data.frame('Sample'=colnames(bulk.metadata),'Sample2'=colnames(bulk.metadata))$Sample),
                      selected = unique(data.frame('Sample'=colnames(bulk.metadata),'Sample2'=colnames(bulk.metadata))$Sample), multiple = TRUE),
          shiny::conditionalPanel(
            id = ns('samples2'),
            ns=ns,
            condition = "input.n_networks >= 2",
            shiny::selectInput(ns('samples2'), 'Samples for GRN #2:', unique(data.frame('Sample'=colnames(bulk.metadata),'Sample2'=colnames(bulk.metadata))$Sample),
                        selected = unique(data.frame('Sample'=colnames(bulk.metadata),'Sample2'=colnames(bulk.metadata))$Sample), multiple = TRUE),
          ),
          shiny::conditionalPanel(
            id = ns('samples3'),
            ns=ns,
            condition = "input.n_networks >= 3",
            shiny::selectInput(ns('samples3'), 'Samples for GRN #3:', unique(data.frame('Sample'=colnames(bulk.metadata),'Sample2'=colnames(bulk.metadata))$Sample),
                        selected = unique(data.frame('Sample'=colnames(bulk.metadata),'Sample2'=colnames(bulk.metadata))$Sample), multiple = TRUE),
          ),
          shiny::conditionalPanel(
            id = ns('samples4'),
            ns=ns,
            condition = "input.n_networks >= 4",
            shiny::selectInput(ns('samples4'), 'Samples for GRN #4:', unique(data.frame('Sample'=colnames(bulk.metadata),'Sample2'=colnames(bulk.metadata))$Sample),
                        selected = unique(data.frame('Sample'=colnames(bulk.metadata),'Sample2'=colnames(bulk.metadata))$Sample), multiple = TRUE),
          ),
          shiny::selectizeInput(ns("targets"), "Target metabolite peaks or from other modality:", multiple = TRUE, choices = NULL),
          shinyjs::disabled(shiny::actionButton(ns('goGRN'), label = 'Start GRN inference')),
          shiny::numericInput(ns("plotConnections"), "Connections to plot:", 5, 0, 100),
          shiny::textInput(ns('plotFileName'), 'File name for plot download', value ='GRNplot.html'),
          shiny::selectInput(ns('plotId'), 'Select plot to download:', 1:4)
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
            )),
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
              shiny::plotOutput(ns('plotUpset')),
              shinyWidgets::dropMenu(
                shinyWidgets::circleButton(ns("downloads_upset"), icon = shiny::icon("download"), status = "success"),
                shiny::tags$div(
                  shiny::tags$h3("Downloads"),
                  shiny::tags$h4("UpSet overlap plot"),
                  shiny::textInput(ns('upsetFileName'), 'File name (overlap)', value ='GRN_OverlapUpset.png'),
                  shiny::numericInput(ns('upsetWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                  shiny::numericInput(ns('upsetHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                  shiny::downloadButton(ns('downloadUpsetPlot'), 'Download overlap')
                ),
                theme = "light-border",
                placement = "right",
                arrow = FALSE
              )
            )
          )
        )
      )
    )
    
  }else{
    NULL
  }
}


#' @rdname BulkPanel_MultiGRNTab
#' @param bulk.intensity.matrix Numeric matrix of bulk sample intensities (features × samples)
#' @param anno Data frame with peak annotation (must include 'display_name' and 'm_z')
#' @param multi_modal Numeric matrix of multi-modal data (features × samples)
#' @export
BulkPanel_MultiGRNTabServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno, multi_modal){
  
  shiny::moduleServer(id, function(input, output, session){
    
    combined_data <- shiny::reactive({
      common_cols <- intersect(
        colnames(bulk.intensity.matrix),
        colnames(multi_modal)
      )
      
      bulk <- bulk.intensity.matrix[, common_cols, drop = FALSE]
      multi <- multi_modal[, common_cols, drop = FALSE]
      combined = as.matrix(rbind(bulk,multi))
      combined = combined[matrixStats::rowMins(combined) !=
                            matrixStats::rowMaxs(combined),]
      
      return(combined)
    })
    
    shiny::observe({
      shiny::req(combined_data())
      shiny::updateSelectizeInput(
        session, "targets", server = TRUE,
        choices = rownames(combined_data()),
        selected = utils::head(rownames(combined_data()))
      )
      
      shiny::updateSelectInput(session, 'samples1', 
                        choices = unique(colnames(combined_data())), 
                        selected = unique(colnames(combined_data())))
      shiny::updateSelectInput(session, 'samples2', 
                        choices = unique(colnames(combined_data())), 
                        selected = unique(colnames(combined_data())))
      shiny::updateSelectInput(session, 'samples3', 
                        choices = unique(colnames(combined_data())), 
                        selected = unique(colnames(combined_data())))
      shiny::updateSelectInput(session, 'samples4', 
                        choices = unique(colnames(combined_data())), 
                        selected = unique(colnames(combined_data())))
      
    })
    
    
    shiny::observe({
      shiny::req(combined_data())
      enable_condition <- length(input[["targets"]]) >= 1 &
        (input[["n_networks"]] < 1 | length(input[["samples1"]]) > 0) &
        (input[["n_networks"]] < 2 | length(input[["samples2"]]) > 0) &
        (input[["n_networks"]] < 3 | length(input[["samples3"]]) > 0) &
        (input[["n_networks"]] < 4 | length(input[["samples4"]]) > 0) &
        (input[["n_networks"]] < 1 | (sum((colnames(combined_data())%in%input[["samples1"]])) > 1)) &
        (input[["n_networks"]] < 2 | (sum((colnames(combined_data())%in%input[["samples2"]])) > 1)) &
        (input[["n_networks"]] < 3 | (sum((colnames(combined_data())%in%input[["samples3"]])) > 1)) &
        (input[["n_networks"]] < 4 | (sum((colnames(combined_data())%in%input[["samples4"]])) > 1))
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
        weightMat <- bulk_utils_infer_GRN(
          intensity_matrix = combined_data(),
          metadata = data.frame('Sample'=colnames(combined_data()),'Sample2'=colnames(combined_data())),
          anno = anno,
          targets = input[["targets"]],
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
      shiny::bindCache(utils::head(combined_data()), data.frame('Sample'=colnames(combined_data()),'Sample2'=colnames(combined_data())), "Sample",
                input[['samples1']],input[['targets']]) |> 
      shiny::bindEvent(input[["goGRN"]])
    GRNresults2 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 2){
        weightMat <- bulk_utils_infer_GRN(
          intensity_matrix = combined_data(),
          metadata = data.frame('Sample'=colnames(combined_data()),'Sample2'=colnames(combined_data())),
          anno = anno,
          targets = input[["targets"]],
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
      shiny::bindCache(utils::head(combined_data()), data.frame('Sample'=colnames(combined_data()),'Sample2'=colnames(combined_data())), "Sample",
                input[['samples2']],input[['targets']]) |> 
      shiny::bindEvent(input[["goGRN"]])
    GRNresults3 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 3){
        weightMat <- bulk_utils_infer_GRN(
          intensity_matrix = combined_data(),
          metadata = data.frame('Sample'=colnames(combined_data()),'Sample2'=colnames(combined_data())),
          anno = anno,
          targets = input[["targets"]],
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
      shiny::bindCache(utils::head(combined_data()), data.frame('Sample'=colnames(combined_data()),'Sample2'=colnames(combined_data())), "Sample",
                input[['samples3']],input[['targets']]) |> 
      shiny::bindEvent(input[["goGRN"]])
    GRNresults4 <- shiny::reactive({
      shinyjs::disable("goGRN")
      if(n_networks() >= 4){
        weightMat <- bulk_utils_infer_GRN(
          intensity_matrix = combined_data(),
          metadata = data.frame('Sample'=colnames(combined_data()),'Sample2'=colnames(combined_data())),
          anno = anno,
          targets = input[["targets"]],
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
      shiny::bindCache(utils::head(combined_data()), data.frame('Sample'=colnames(combined_data()),'Sample2'=colnames(combined_data())), "Sample",
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
    output[['downloadUpsetPlot']] <- utils_create_download_plot_handler(
      plot_func = upsetPlot,
      filename_func = function() input[['upsetFileName']],
      width_func = function() input[['upsetWidth']],
      height_func = function() input[['upsetHeight']]
    )
    
    output[['download']] <- shiny::downloadHandler(
      filename = function() {input[['plotFileName']]},
      content = function(file) {
        GRNplot <- get(paste0("GRNplot", input[["plotId"]]))()
        GRNplot |> visNetwork::visSave(file)
      }
    )
    
  })
}