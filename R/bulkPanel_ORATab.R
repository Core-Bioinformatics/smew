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
              shiny::tags$h4("Pathway Categories"),
              shiny::textInput(ns('categoriesFileName'),'File name for download', value ='pathwayCategories.png', placeholder = 'pathwayCategories.png'),
              shiny::numericInput(ns('categoriesWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('categoriesHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadCategories'), 'Download category figure'),
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
          shiny::tags$h2("Pathway Category Plot"),
          shiny::plotOutput(ns('pathwayCategories')),

        )
      ),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput(ns('selectedClassification'),label = 'Select pathway categories to show in network',choices = c('None','Top Level','More Granular'),selected = 'Top Level',multiple = F),
        ),
        shiny::mainPanel(
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
#' @param organism Character string specifying the organism (for KEGG mapping)
#' @export
BulkPanel_ORATabServer <- function(id, bulk.intensity.matrix, bulk.metadata, de.results, anno, organism){

  # check whether inputs (other than id) are reactive or not
  stopifnot({
    shiny::is.reactive(de.results)
  })

  shiny::moduleServer(id, function(input, output, session){

    get_ORA <- shiny::reactive({
      bulk_utils_execute_ora(de_peaks = de.results()$DE()$DEtableSubset,
                  path_dict = NULL,
                  background = input[['background_selector']],
                  min_path_size = input[['min_pathway_size']],
                  ora_pvalue_cutoff = input[['ora_pvalue_cutoff']],
                  min_pathway_hits = input[['min_pathway_hits']],
                  organism=organism,
                  anno=anno,
                  kegg_db = kegg_db)
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
      pathway.list = list()
      # could consider switching to overlap in the DE metabolites rather than just overall
      for (pathway in unique(kegg_db$pathway_name)){
        pathway.list[[pathway]]=unique(kegg_db[kegg_db$pathway_name==pathway,]$compound_id)
      }
      # might be issues here if a pathway comes up as up and down
      nodes = data.frame(label=ora$pathway,id=ora$pathway,shape='circle',color=ora$lfc,font.color='white')
      nodes$label = stringr::str_wrap(nodes$label,10)
      edges = data.frame(t(utils::combn(names(pathway.list), 2)))
      edges = edges[edges$X1!=edges$X2,]
      colnames(edges)=c('from','to')
      weight.vector = c()
      for (i in 1:nrow(edges)){
        from.met = pathway.list[[edges[i,'from']]]
        to.met = pathway.list[[edges[i,'to']]]
        weight.vector = c(weight.vector,(length(intersect(from.met,to.met))/length(union(from.met,to.met))))
        #  weight.vector = c(weight.vector,(length(intersect(from.met,to.met))))
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
      if (input[['selectedClassification']]!='None'){
        kegg_classification_sub = kegg_classification[kegg_classification$pathway_name%in%ora$pathway,]
        if (input[['selectedClassification']]=='Top Level'){
          kegg_classification_sub$category = kegg_classification_sub$category1
        } else {
          kegg_classification_sub$category = kegg_classification_sub$category2
        }
        classification_edges = data.frame('from'=kegg_classification_sub$category,'to'=kegg_classification_sub$pathway_name,value=min(edges$value))
        edges = rbind(edges,classification_edges)
        classification_nodes = data.frame('label'=unique(kegg_classification_sub$category),'id'=unique(kegg_classification_sub$category),
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
      significant.pathways = get_ORA()[get_ORA()$FDR<input[['ora_pvalue_cutoff']],]
      kegg_classification_sub = kegg_classification
      kegg_classification_sub$pathway = kegg_classification_sub$pathway_name
      kegg_classification_sub$category1 = factor(kegg_classification_sub$category1)
      kegg_classification_sub$category2 = factor(kegg_classification_sub$category2)
      kegg_classification_sub$pathway_id = paste0('map',kegg_classification_sub$pathway_id)
      significant.pathways = merge(significant.pathways,kegg_classification_sub)
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
      pathwayCategories()
    })

    output[['downloadCategories']] <- shiny::downloadHandler(
      filename = function() { input[['categoriesFileName']] },
      content = function(file) {
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
