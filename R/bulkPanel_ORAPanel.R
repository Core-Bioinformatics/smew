BulkORAPanelUI <- function(id, bulk.metadata, show = TRUE){
  ns <- NS(id)
  
  if(show){
    tabPanel(
      'ORA',
      sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
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
          )
        ),
        
        #Main panel for displaying table of enriched pathways
        mainPanel(
          DT::dataTableOutput(ns('data')),
          plotOutput(ns('oraVolcano'),click=ns('plot_click')),
          tableOutput(ns('oraVolcanoData')),
          plotOutput(ns('pathwayCategories'))
        )
      )
    )
  }else{
    NULL
  }
}

#' @rdname DEsummaryPanel
#' @export
BulkORAPanelServer <- function(id, bulk.expression.matrix, bulk.metadata, DEresults, anno){
  
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(DEresults)
    is.reactive(bulk.expression.matrix)
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
        dplyr::filter(FDR<0.05) |>
        dplyr::select(-metabolites) |>
        DT::datatable() %>%
        DT::formatSignif(columns = c('Raw.p', 'Holm.p','FDR'), digits = 3)
    })
    
    output[['data']] <- DT::renderDataTable(dataTable())
    
    output[['oraVolcano']] <- renderPlot({
      ora_volcano_plot(get_ORA(),input[['ora_pvalue_cutoff']])
    })
    
    output[['oraVolcanoData']] <- renderTable({
      req(input[['plot_click']])
      data = get_ORA() |> dplyr::mutate(`-log10pval` = -log10(.data$FDR),
                                        lfc = ifelse(.data$direction=='up',log2(.data$hits/.data$expected),-log2(.data$hits/.data$expected)))
      nearPoints(df = data, coordinfo = input[['plot_click']], threshold = 20, maxpoints = 10)
    }, digits = 4)
    
    output[['pathwayCategories']] <- renderPlot({
      significant.pathways = get_ORA()[get_ORA()$FDR<0.05,]
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
    
  })
}
