#' @rdname BulkORAPanel
#' @export
BulkORAPanelUI <- function(id, bulk.metadata, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'ORA',
      sidebarLayout(

        # Sidebar panel for inputs ----
        sidebarPanel(
          dropMenu(
            circleButton(ns("info_differential_analysis"), icon = icon("info"),status = "success"),
            tags$div(
              tags$h3("Over-representation analysis"),
              tags$ul(
                tags$li("This tab performs over-representation analysis on KEGG pathways."),
                tags$li("The differential metabolites are taken from the last time the Differential intensity analysis tab was run (i.e. the Start differential intensity analysis button was pressed)."),
                tags$li("Pathways with small numbers of present metabolites or significant can be discarded using the options below."),
                tags$li("The p-value threshold for calling a pathway significant can also be adjusted."),
                tags$li("A table of significant pathways is shown, followed by a volcano plot and a separation of pathways into KEGG hierarchy categories and sub-categories."),
                tags$li("Finally, a network is displayed showing the significant pathways with connections between them proportional to the number of shared metabolites between them."),
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
          # add slider to select the minimum number of pathways to be considered before peak mapping
          sliderInput(
            inputId = ns("min_pathway_size"),
            label = "Discard pathways with overall size of less than x number of metabolites:",
            min = 1,
            max = 10,
            value = 3,
            step = 1
          ),
          # slider input for the number of hits in ORA to retain the pathway
          sliderInput(
            inputId = ns("min_pathway_hits"),
            label = "Minimum number of hits in a pathway",
            min = 1,
            max = 10,
            value = 2,
            step = 1
          ),
          # slider input for p-value cutoff for filtering out pathways
          sliderInput(
            inputId = ns("ora_pvalue_cutoff"),
            label = "P-value cutoff for pathway filtering",
            min = 0.01,
            max = 1,
            value = 0.05,
            step = 0.01
          ),
          # button to start ORA
          actionButton(
            inputId = ns("submit_from_ora"),
            label = "Submit",
            icon = icon("play")
          ),
          div(style = "margin-top:10px"),
          tags$b("Pathway selection"),
          div("\nSelect pathways of interest by clicking on the corresponds rows in the table\n"),
          div(style="margin-bottom:10px"),
          actionButton(ns('resetSelection'), label = "Reset row selection"),
          div(style="margin-bottom:10px"),
          actionButton(ns('selectTop10'), label = "Select top 10 pathways"),
          div(style = "margin-top:10px"),
          tags$b("Pathway network"),
          selectInput(ns('selectedClassification'),label = 'Select pathway categories to show in network',choices = c('None','Top Level','More Granular'),selected = 'Top Level',multiple = F),
          div(style = "margin-top:10px"),
          dropMenu(
            circleButton(ns("downloads"), icon = icon("download"),status = "success"),
            tags$div(
              tags$h3("Downloads"),
              fluidRow(
                column(5,offset=0,
                       tags$h4("ORA Table"),
                       textInput(ns('tableFileName'),'File name for download', value ='ORATable.csv', placeholder = 'ORATable.csv'),
                       downloadButton(ns('downloadTable'), 'Download ORA table'),
                       ),
                column(5,offset=1,
                       tags$h4("Volcano Plot"),
                       textInput(ns('volcanoFileName'),'File name for download', value ='volcano.png', placeholder = 'volcano.png'),
                       numericInput(ns('volcanoWidth'),value = 8,label = 'Width of downloaded figure (in inches)',min = 1,max = 50,step = 1),
                       numericInput(ns('volcanoHeight'),value = 6,label = 'Height of downloaded figure (in inches)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadVolcano'), 'Download volcano plot'),
                ),

            ),
            fluidRow(
              column(5,offset=0,
                     tags$h4("Pathway Categories"),
                     textInput(ns('categoriesFileName'),'File name for download', value ='pathwayCategories.png', placeholder = 'pathwayCategories.png'),
                     numericInput(ns('categoriesWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                     numericInput(ns('categoriesHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                     downloadButton(ns('downloadCategories'), 'Download category figure')
              ),
              column(5,offset=1,
                     tags$h4("Pathway Network"),
                     textInput(ns('networkFileName'),'File name for download', value ='network.html', placeholder = 'network.html'),
                     downloadButton(ns('downloadNetwork'), 'Download network')

              ),

            )),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
        ),

        #Main panel for displaying table of enriched pathways
        mainPanel(
          DT::dataTableOutput(ns('data')),
          plotOutput(ns('oraVolcano'),click=ns('plot_click')),
          tableOutput(ns('oraVolcanoData')),
          plotOutput(ns('pathwayCategories')),
          visNetwork::visNetworkOutput(ns('ORAnetwork'),height="600")
        )
      )
    )
  }else{
    NULL
  }
}

#' @rdname BulkORAPanel
#' @export
BulkORAPanelServer <- function(id, bulk.intensity.matrix, bulk.metadata, DEresults, anno){

  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(DEresults)
    is.reactive(bulk.intensity.matrix)
    is.reactive(bulk.metadata)
    !is.reactive(anno)
  })

  moduleServer(id, function(input, output, session){

    get_ORA <- reactive({
      execute_ora(de_peaks = DEresults()$DE()$DEtableSubset,
                  path_dict = NULL,
                  background = input[['background_selector']],
                  min_path_size = input[['min_pathway_size']],
                  ora_pvalue_cutoff = input[['ora_pvalue_cutoff']],
                  min_pathway_hits = input[['min_pathway_hits']])
    }) %>% bindEvent(input[["submit_from_ora"]])

    dataTable <- reactive({
      get_ORA() |>
        dplyr::filter(FDR<input[['ora_pvalue_cutoff']]) |>
        dplyr::select(-metabolites) |>
        DT::datatable() %>%
        DT::formatSignif(columns = c('Raw.p', 'Holm.p','FDR'), digits = 3)
    })

    #Output selected pathways
    selectedPathways <- reactive({
      ora_table <- get_ORA() |>
        dplyr::filter(FDR<input[['ora_pvalue_cutoff']])
      ora_table$pathway[input$data_rows_selected]
    })

    proxy = DT::dataTableProxy('data')

    observe({proxy %>% DT::selectRows(NULL)}) %>%
      bindEvent(input[['resetSelection']])

    observe({proxy %>% DT::selectRows(selected = 1:10)}) %>%
      bindEvent(input[['selectTop10']])


    oraNetwork <- reactive({
      ora <- get_ORA() |>
        dplyr::filter(FDR<input[['ora_pvalue_cutoff']]) |>
        dplyr::mutate(`-log10pval` = -log10(.data$FDR),
                      lfc = ifelse(.data$direction=='up',log2(.data$hits/.data$expected),-log2(.data$hits/.data$expected)))
      pathway.list = list()
      for (pathway in unique(kegg_db$pathway_name)){
        pathway.list[[pathway]]=unique(kegg_db[kegg_db$pathway_name==pathway,]$compound_id)
      }
      nodes = data.frame(label=ora$pathway,id=ora$pathway,shape='circle',color=ora$lfc,font.color='white')
      nodes$label = stringr::str_wrap(nodes$label,10)
      edges = data.frame(t(combn(names(pathway.list), 2)))
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
      return(list('nodes'=nodes,'edges'=edges))
    })
    oraNetworkCategories <- reactive({
      ora <- get_ORA() |>
        dplyr::filter(FDR<input[['ora_pvalue_cutoff']]) |>
        dplyr::mutate(`-log10pval` = -log10(.data$FDR),
                      lfc = ifelse(.data$direction=='up',log2(.data$hits/.data$expected),-log2(.data$hits/.data$expected)))
      network = oraNetwork()
      nodes = network$nodes
      edges = network$edges
      if (input[['selectedClassification']]!='None'){
        kegg_classification = kegg_classification[kegg_classification$pathway_name%in%ora$pathway,]
        if (input[['selectedClassification']]=='Top Level'){
          kegg_classification$category = kegg_classification$category1
        } else {
          kegg_classification$category = kegg_classification$category2
        }
        classification_edges = data.frame('from'=kegg_classification$category,'to'=kegg_classification$pathway_name,value=min(edges$value))
        edges = rbind(edges,classification_edges)
        classification_nodes = data.frame('label'=unique(kegg_classification$category),'id'=unique(kegg_classification$category),
                                          'shape'='box','color'='grey',font.color='white')
        classification_nodes = classification_nodes[!(classification_nodes$id %in% nodes$id),]
        nodes = rbind(nodes,classification_nodes)
        print(tail(nodes))

      }
      return(visNetwork::visNetwork(nodes,edges) |>
        visNetwork::visPhysics(solver = "forceAtlas2Based",
                               forceAtlas2Based = list(gravitationalConstant = -100)))

    })

    pathwayCategories <- reactive({
      significant.pathways = get_ORA()[get_ORA()$FDR<input[['ora_pvalue_cutoff']],]
      kegg_classification$pathway = kegg_classification$pathway_name
      kegg_classification$category1 = factor(kegg_classification$category1)
      kegg_classification$category2 = factor(kegg_classification$category2)
      kegg_classification$pathway_id = paste0('map',kegg_classification$pathway_id)
      significant.pathways = merge(significant.pathways,kegg_classification)
      return(ggplot2::ggplot(significant.pathways,aes(x=category2,y=-log10(FDR),color=category2))+
               ggplot2::geom_point()+
               ggplot2::facet_wrap(~significant.pathways$category1,scales = 'free_x',ncol=2)+
               ggplot2::theme_classic()+
               ggplot2::xlab('Pathway sub-category')+
               ggplot2::scale_size_binned(range=c(0,max(-log10(significant.pathways$FDR))),n.breaks=10)+
               ggrepel::geom_label_repel(aes(label = pathway))+
               ggplot2::theme(legend.position="none")
      )

    })

    output[['data']] <- DT::renderDataTable(dataTable())

    output[['downloadTable']] <- downloadHandler(
      filename = function() {
        paste(input[['tableFileName']])
      },
      content = function(file) {
        utils::write.csv(x = get_ORA() |>
                           dplyr::filter(FDR<input[['ora_pvalue_cutoff']]), file = file, row.names = FALSE)
      }
    )

    output[['oraVolcano']] <- renderPlot({
      ora_volcano_plot(get_ORA(),input[['ora_pvalue_cutoff']],selectedPathways())
    })

    output[['downloadVolcano']] <- downloadHandler(
      filename = function() { input[['volcanoFileName']] },
      content = function(file) {
        ggsave(file, plot = ora_volcano_plot(get_ORA(),input[['ora_pvalue_cutoff']]), dpi = 300,
               width=input[['volcanoWidth']],height=input[['volcanoHeight']])
      }
    )

    output[['oraVolcanoData']] <- renderTable({
      req(input[['plot_click']])
      data = get_ORA() |> dplyr::mutate(`-log10pval` = -log10(.data$FDR),
                                        lfc = ifelse(.data$direction=='up',log2(.data$hits/.data$expected),-log2(.data$hits/.data$expected)))
    }, digits = 4)

    output[['pathwayCategories']] <- renderPlot({
      pathwayCategories()
    })

    output[['downloadCategories']] <- downloadHandler(
      filename = function() { input[['categoriesFileName']] },
      content = function(file) {
        ggsave(file, plot = pathwayCategories(), dpi = 300,
               width=input[['categoriesWidth']],height=input[['categoriesHeight']])
      }
    )

    output[['ORAnetwork']] <- visNetwork::renderVisNetwork(oraNetworkCategories())

    output[['downloadNetwork']] <- downloadHandler(
      filename = function() {input[['networkFileName']]},
      content = function(file) {
        oraNetworkCategories() %>% visNetwork::visSave(file)
      }
    )


  })
}
