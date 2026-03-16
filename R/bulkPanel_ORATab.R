#' Performs pathway over-representation analysis
#'
#' @description UI and server logic for over-representation analysis (ORA) of KEGG pathway enrichment. Tests whether differentially expressed metabolites are enriched in specific pathways.
#'
#' @details
#' \itemize{
#'   \item{Pathway filtering by size and hit counts}
#'   \item{P-value based significance thresholding}
#'   \item{Interactive volcano plot for pathway visualization}
#'   \item{KEGG hierarchy categorization of pathways}
#'   \item{Pathway co-occurrence network showing metabolite sharing}
#'   \item{Accepts DE results from the differential intensity analysis module and maps significant metabolites to KEGG pathways. Supports visualization of pathway relationships through network graphs where edge weights represent Jaccard similarity of pathway metabolite composition.}
#' }
#' @name BulkPanel_ORATab
#' @rdname BulkPanel_ORATab
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame with bulk sample metadata
#' @param show Logical; whether to render the panel (default: TRUE)
#' @return A shiny::tabPanel containing the UI elements for the ORA panel
#' @export
BulkPanel_ORATabUI <- function(id, bulk.metadata, show = TRUE){
  ns <- shiny::NS(id)

  if(show){
    shiny::tabPanel(
      'Pathway ORA',
      shiny::sidebarLayout(

        # Sidebar panel for inputs ----
        shiny::sidebarPanel(
          shinyWidgets::dropMenu(
            shinyWidgets::circleButton(ns("info_ora"), icon = shiny::icon("info-circle"), status = "info"),
            shiny::tags$h3("ORA"),
            shiny::tags$ul(
              shiny::tags$li("This tab performs over-representation analysis (ORA) on KEGG pathways using significant metabolites from the differential intensity analysis (DA tab)."),
              shiny::tags$li("Outputs include: a table of significant pathways, volcano plot, pathway category plot, and a network graph of pathway relationships."),
              shiny::tags$li("Settings allow filtering by pathway size, hit count, and p-value threshold.")
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
          shiny::tags$h3("ORA Settings"),
          shiny::sliderInput(
            inputId = ns("min_pathway_size"),
            label = "Discard pathways with overall size of less than x number of metabolites:",
            min = 1,
            max = 10,
            value = 3,
            step = 1
          ),
          shiny::sliderInput(
            inputId = ns("min_pathway_hits"),
            label = "Minimum number of hits in a pathway",
            min = 1,
            max = 10,
            value = 2,
            step = 1
          ),
          shiny::sliderInput(
            inputId = ns("ora_pvalue_cutoff"),
            label = "P-value cutoff for pathway filtering",
            min = 0.01,
            max = 1,
            value = 0.05,
            step = 0.01
          ),
          shiny::actionButton(
            inputId = ns("submit_from_ora"),
            label = "Start ORA",
            class = "btn-primary"),
          shiny::div(style = "margin-top:20px"),
          shinyWidgets::dropMenu(
            shinyWidgets::circleButton(ns("downloads_ora"), icon = shiny::icon("download"), status = "success"),
            shiny::tags$div(
              shiny::tags$h3("Download ORA Results"),
              shiny::tags$h4("ORA Table"),
              shiny::textInput(ns('tableFileName'),'File name for download', value ='ORATable.csv', placeholder = 'ORATable.csv'),
              shiny::downloadButton(ns('downloadTable'), 'Download ORA table'),
              shiny::tags$hr(),
              shiny::tags$h4("Volcano Plot"),
              shiny::textInput(ns('volcanoFileName'),'File name for download', value ='volcano.png', placeholder = 'volcano.png'),
              shiny::numericInput(ns('volcanoWidth'),value = 8,label = 'Width of downloaded figure (in inches)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('volcanoHeight'),value = 6,label = 'Height of downloaded figure (in inches)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadVolcano'), 'Download volcano plot'),
              shiny::tags$hr(),
              shiny::tags$h4("Pathway Network"),
              shiny::textInput(ns('networkFileName'),'File name for download', value ='network.html', placeholder = 'network.html'),
              shiny::downloadButton(ns('downloadNetwork'), 'Download network')
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          )
        ),

        #Main panel for displaying table of enriched pathways
        shiny::mainPanel(
          shiny::tags$h2("Significant Pathways Table"),
          DT::dataTableOutput(ns('data')),
          shiny::tags$h2("Pathway Enrichment Volcano Plot"),
          plotly::plotlyOutput(ns('oraVolcano')),
          shiny::tableOutput(ns('oraVolcanoData')),

        )
      ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::uiOutput(ns('classificationUI')),
        ),
        shiny::mainPanel(
          shiny::uiOutput(ns('categoryPlotUI')),
          shiny::tags$h2("Pathway Co-occurrence Network"),
          visNetwork::visNetworkOutput(ns('ORAnetwork'),height="600")
        )
      )
    )
  }else{
    NULL
  }
}

#' @rdname BulkPanel_ORATab
#' @param bulk.intensity.matrix Numeric matrix of bulk sample intensities (features × samples)
#' @param de.results Reactive expression or data frame of differential expression results
#' @param anno Data frame with peak annotation (must include 'display_name' and 'm_z')
#' @param pathway_table Optional pathway table with columns PathwayID, PathwayName, and MetaboliteIDs.
#' @param pathway_classification Optional pathway classification table. If NULL, category-specific UI/plots are hidden.
#' @export
BulkPanel_ORATabServer <- function(id, bulk.intensity.matrix, bulk.metadata, de.results, anno, pathway_table = NULL, pathway_classification = NULL){

  # check whether inputs (other than id) are reactive or not
  stopifnot({
    shiny::is.reactive(de.results)
  })

  shiny::moduleServer(id, function(input, output, session){
    ns <- session$ns

    has_pathway_table <- !is.null(pathway_table) &&
      is.data.frame(pathway_table) &&
      all(c('PathwayID', 'PathwayName', 'MetaboliteIDs') %in% colnames(pathway_table))

    parse_ids <- function(x) {
      vals <- unlist(strsplit(as.character(x), "[,;|]"))
      vals <- trimws(vals)
      vals[vals != "" & !is.na(vals)]
    }

    has_classification <- !is.null(pathway_classification) &&
      is.data.frame(pathway_classification) &&
      all(c('PathwayName', 'PathwayID', 'Category1', 'Category2') %in% colnames(pathway_classification))

    output$classificationUI <- shiny::renderUI({
      if (has_classification) {
        shiny::tagList(
          shiny::selectInput(ns('selectedClassification'),
                             label = 'Select pathway categories to show in network',
                             choices = c('None','Top Level','More Granular'),
                             selected = 'Top Level',
                             multiple = FALSE),
          shiny::tags$hr(),
          shiny::tags$h4("Pathway Categories"),
          shiny::textInput(ns('categoriesFileName'),'File name for download', value ='pathwayCategories.png', placeholder = 'pathwayCategories.png'),
          shiny::numericInput(ns('categoriesWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
          shiny::numericInput(ns('categoriesHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
          shiny::downloadButton(ns('downloadCategories'), 'Download category figure')
        )
      } else {
        shiny::tags$p("No pathway classification table supplied; category overlays are disabled.")
      }
    })

    output$categoryPlotUI <- shiny::renderUI({
      if (has_classification) {
        shiny::tagList(
          shiny::tags$h2("Pathway Category Plot"),
          shiny::plotOutput(ns('pathwayCategories'), height = "600px")
        )
      } else {
        NULL
      }
    })

    get_ORA <- shiny::reactive({
      shiny::validate(shiny::need(has_pathway_table, "Pathway enrichment is unavailable because pathway_table was not supplied."))
      bulk_utils_execute_ora(
        de_peaks = de.results()$DE()$DEtableSubset,
        anno = anno,
        pathway_db = pathway_table,
        min_path_size = input[['min_pathway_size']],
        min_pathway_hits = input[['min_pathway_hits']]
      )
    }) |> shiny::bindEvent(input[["submit_from_ora"]])

    dataTable <- shiny::reactive({
      get_ORA() |>
        dplyr::filter(.data$FDR<input[['ora_pvalue_cutoff']]) |>
        dplyr::select(-.data$metabolites) |>
        DT::datatable() |>
        DT::formatSignif(columns = c('Raw.p', 'Holm.p','FDR'), digits = 3)
    }) |> shiny::bindEvent(input[["submit_from_ora"]])

    #Output selected pathways
    selectedPathways <- shiny::reactive({
      ora_table <- get_ORA() |>
        dplyr::filter(.data$FDR<input[['ora_pvalue_cutoff']])
      ora_table$pathway[input$data_rows_selected]
    })

    proxy = DT::dataTableProxy('data')

    # observe({proxy |> DT::selectRows(NULL)}) |>
    #   bindEvent(input[['resetSelection']])
    #
    # observe({proxy |> DT::selectRows(selected = 1:10)}) |>
    #   bindEvent(input[['selectTop10']])
    #

    oraNetwork <- shiny::reactive({
      ora <- get_ORA() |>
        dplyr::filter(.data$FDR<input[['ora_pvalue_cutoff']]) |>
        dplyr::mutate(`-log10pval` = -log10(.data$FDR),
                      lfc = ifelse(.data$direction=='up',log2(.data$hits/.data$expected),-log2(.data$hits/.data$expected)))
      pathway_names <- ifelse(is.na(pathway_table$PathwayName) | pathway_table$PathwayName == "",
                              as.character(pathway_table$PathwayID),
                              as.character(pathway_table$PathwayName))
      pathway.list <- split(pathway_table$MetaboliteIDs, pathway_names)
      pathway.list <- lapply(pathway.list, function(vals) unique(parse_ids(paste(vals, collapse = ','))))
      pathway.list <- pathway.list[unique(ora$pathway)]
      # might be issues here if a pathway comes up as up and down
      nodes = data.frame(label=ora$pathway,id=ora$pathway,shape='circle',color=ora$lfc,font.color='white')
      nodes$label = stringr::str_wrap(nodes$label,10)
      if (length(pathway.list) >= 2) {
        edges = data.frame(t(utils::combn(names(pathway.list), 2)))
        edges = edges[edges$X1!=edges$X2,]
        colnames(edges)=c('from','to')
      } else {
        edges <- data.frame(from = character(), to = character())
      }
      weight.vector = c()
      if (nrow(edges) > 0) {
        for (i in seq_len(nrow(edges))){
          from.met = pathway.list[[edges[i,'from']]]
          to.met = pathway.list[[edges[i,'to']]]
          weight.vector = c(weight.vector,(length(intersect(from.met,to.met))/length(union(from.met,to.met))))
        }
      }
      edges$value = weight.vector
      edges = edges[edges$value!=0,]
      nodes$color = scales::col_numeric('RdYlBu',-ceiling(max(abs(nodes$color))):ceiling(max(abs(nodes$color))))(nodes$color)
      return(list('nodes'=nodes,'edges'=edges,'ora'=ora))
    }) |> shiny::bindEvent(input[["submit_from_ora"]])

    oraNetworkCategories <- shiny::reactive({
      network = oraNetwork()
      ora <- network$ora |>
        dplyr::filter(.data$FDR<input[['ora_pvalue_cutoff']]) |>
        dplyr::mutate(`-log10pval` = -log10(.data$FDR),
                      lfc = ifelse(.data$direction=='up',log2(.data$hits/.data$expected),-log2(.data$hits/.data$expected)))

      nodes = network$nodes
      edges = network$edges
      if (has_classification && !is.null(input[['selectedClassification']]) && input[['selectedClassification']]!='None'){
        pathway_classification_sub = pathway_classification[pathway_classification$PathwayName %in% ora$pathway,]
        if (input[['selectedClassification']]=='Top Level'){
          pathway_classification_sub$category = pathway_classification_sub$Category1
        } else {
          pathway_classification_sub$category = pathway_classification_sub$Category2
        }
        min_edge <- if (nrow(edges) > 0) min(edges$value) else 1
        classification_edges = data.frame('from'=pathway_classification_sub$category,'to'=pathway_classification_sub$PathwayName,value=min_edge)
        edges = rbind(edges,classification_edges)
        classification_nodes = data.frame('label'=unique(pathway_classification_sub$category),'id'=unique(pathway_classification_sub$category),
                                          'shape'='box','color'='grey',font.color='white')
        classification_nodes = classification_nodes[!(classification_nodes$id %in% nodes$id),]
        nodes = rbind(nodes,classification_nodes)

      }
      return(visNetwork::visNetwork(nodes,edges) |>
        visNetwork::visPhysics(solver = "forceAtlas2Based",
                               forceAtlas2Based = list(gravitationalConstant = -100)))

    })

    # if there's more than a certain number, just show the dots and allow to hover!
    pathwayCategories <- shiny::reactive({
      shiny::req(has_classification)
      significant.pathways = get_ORA()[get_ORA()$FDR<input[['ora_pvalue_cutoff']],]
      pathway_classification_sub = pathway_classification
      pathway_classification_sub$pathway = pathway_classification_sub$PathwayName
      pathway_classification_sub$category1 = factor(pathway_classification_sub$Category1)
      pathway_classification_sub$category2 = factor(pathway_classification_sub$Category2)
      pathway_classification_sub$pathway_id = paste0('map',pathway_classification_sub$PathwayID)
      significant.pathways = merge(significant.pathways,pathway_classification_sub)
      return(ggplot2::ggplot(significant.pathways,ggplot2::aes(x=.data$category2,y=-log10(.data$FDR),color=.data$category2))+
               ggplot2::geom_point()+
               ggplot2::facet_wrap(~.data$category1,scales = 'free_x',ncol=2)+
               ggplot2::theme_classic()+
               ggplot2::xlab('Pathway sub-category')+
               ggplot2::scale_size_binned(range=c(0,max(-log10(significant.pathways$FDR))),n.breaks=10)+
               ggrepel::geom_label_repel(ggplot2::aes(label = .data$pathway))+
               ggplot2::theme(legend.position="none")
      )

    })

    volcano <- shiny::reactive({
      return(bulk_utils_ora_volcano_plot(get_ORA(),input[['ora_pvalue_cutoff']],selectedPathways()))

    }) |> shiny::bindEvent(input[["submit_from_ora"]])

    output[['data']] <- DT::renderDataTable(dataTable())

    output[['downloadTable']] <- shiny::downloadHandler(
      filename = function() {
        paste(input[['tableFileName']])
      },
      content = function(file) {
        utils::write.csv(x = get_ORA() |>
                           dplyr::filter(.data$FDR<input[['ora_pvalue_cutoff']]), file = file, row.names = FALSE)
      }
    )


    output[['oraVolcano']] <-
      plotly::renderPlotly(volcano()$volcano)

    output[['downloadVolcano']] <- shiny::downloadHandler(
      filename = function() { input[['volcanoFileName']] },
      content = function(file) {
        ggplot2::ggsave(file, plot = bulk_utils_ora_volcano_plot(get_ORA(),input[['ora_pvalue_cutoff']],selectedPathways())$volcano,
                        dpi = 300, width=input[['volcanoWidth']],height=input[['volcanoHeight']])
      }
    )

    # output[['oraVolcanoData']] <- renderTable({
    #   req(input[['plot_click']])

#      data = nearPoints(get_ORA(),coordinfo = input$plot_click,maxpoints=1) |> dplyr::mutate(`-log10pval` = -log10(.data$FDR),
#                                        lfc = ifelse(.data$direction=='up',log2(.data$hits/.data$expected),-log2(.data$hits/.data$expected)))
    # nearPoints(volcano()$data,coordinfo = input$plot_click, maxpoints=1)
    # }, digits = 4)

    output[['pathwayCategories']] <- shiny::renderPlot({
      shiny::req(has_classification)
      pathwayCategories()
    })

    output[['downloadCategories']] <- shiny::downloadHandler(
      filename = function() { input[['categoriesFileName']] },
      content = function(file) {
        shiny::req(has_classification)
        ggplot2::ggsave(file, plot = pathwayCategories(), dpi = 300,
               width=input[['categoriesWidth']],height=input[['categoriesHeight']])
      }
    )

    output[['ORAnetwork']] <- visNetwork::renderVisNetwork(oraNetworkCategories())

    output[['downloadNetwork']] <- shiny::downloadHandler(
      filename = function() {input[['networkFileName']]},
      content = function(file) {
        oraNetworkCategories() |> visNetwork::visSave(file)
      }
    )


  })
}
