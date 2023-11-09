RegionDEpanelUI <- function(id, bulk.metadata, full.metadata, show = TRUE){
  ns <- NS(id)
  
  if(show){
    tabPanel(
      'Differential expression',
      shinyjs::useShinyjs(),
      sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
          
          selectInput(
            inputId = ns("regionToGroupOn"),
            label = "Select region to group on",
            choices = c(colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata)))]),
            selected = c(colnames(full.metadata)[!(colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata)))])[1]
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
          actionButton(ns('goDE'), label = 'Start DE'),
          
          #download file name and button
          textInput(ns('fileName'),'File name for download', value ='DEset.csv', placeholder = 'DEset.csv'),
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

#' @rdname DEpanel
#' @export
RegionDEpanelServer <- function(id, full.expression.matrix, full.metadata, bulk.metadata, region.clusters, anno){
  # check whether inputs (other than id) are reactive or not
  stopifnot({
    is.reactive(bulk.expression.matrix)
    is.reactive(bulk.metadata)
    !is.reactive(anno)
  })
  
  moduleServer(id, function(input, output, session){
    
    pseudobulk_samples <- reactive({
#      if (!is.null(region.clusters())){full.metadata = merge(full.metadata,region.clusters,all.x=T)}
      print("I'm about to create the bulk expression matrix")
      print(input[['regionToGroupOn']])
      pseudobulked <- create_bulk_exp_regions(full.expression.matrix,
                              full.metadata,
                              bulk.metadata,
                              colnames(bulk.metadata)[1],
                              full.metadata[,input[['regionToGroupOn']]],
                              minimum.pixels = input[['minimumPixelsPerSample']])
#      if (is.null(region.clusters())){
        choices = colnames(full.metadata)[colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata))]
        selected = colnames(full.metadata)[colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata))][1]
#      } else {
#        choices = c('cluster',colnames(full.metadata)[colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata))])
#        selected = 'cluster'
#      }
      print(choices)
      updateSelectInput(session, 'regionToGroupOn', choices = choices,
                        selected=selected)
      updateSelectInput(session, 'variable1', choices = unique(pseudobulked$metadata[[input[["condition"]]]]))
      updateSelectInput(session, 'variable2', choices = unique(pseudobulked$metadata[[input[["condition"]]]]),
                        selected = unique(pseudobulked$metadata[[input[["condition"]]]])[2])
      
      return(pseudobulked)
    }) %>% bindEvent(input[["pseudoBulk"]])

    observe({
    # updateSelectInput(session, 'condition', choices = ifelse(is.null(region.clusters),
    #                                                          colnames(full.metadata)[!(colnames(full.metadata) %in% c('spot_id','x','y'))],
    #                                                          c('cluster',colnames(full.metadata)[!(colnames(full.metadata) %in% c('spot_id','x','y'))])),
    #                   selected=ifelse(is.null(region.clusters),colnames(full.metadata)[ncol(full.metadata)],'cluster'))
    # updateSelectInput(session, 'variable1', choices = ifelse(input[["condition"]]=='cluster',unique(region.clusters$cluster),unique(full.metadata[[input[["condition"]]]])))
    # updateSelectInput(session, 'variable2', choices = ifelse(input[["condition"]]=='cluster',unique(region.clusters$cluster),unique(full.metadata[[input[["condition"]]]])),
    #                   selected = ifelse(input[["condition"]]=='cluster',unique(region.clusters$cluster)[2],unique(full.metadata[[input[["condition"]]]])[2]))
    # if (is.null(region.clusters())){
    #   choices = colnames(full.metadata)[colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata))]
    #   selected = colnames(full.metadata)[colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata))][1]
    # } else {
    #   choices = c('cluster',colnames(full.metadata)[colnames(full.metadata)%in%c('spot_id','x','y',colnames(bulk.metadata))])
    #   selected = 'cluster'
    # }
    # print(choices)
    # updateSelectInput(session, 'regionToGroupOn', choices = choices,
    #                   selected=selected)
    # updateSelectInput(session, 'variable1', choices = unique(pseudobulk_samples()$metadata[[input[["condition"]]]]))
    # updateSelectInput(session, 'variable2', choices = unique(pseudobulk_samples()$metadata[[input[["condition"]]]]),
    #                   selected = unique(pseudobulk_samples()$metadata[[input[["condition"]]]])[2])

    })

    observe({
      condition.indices <- full.metadata[[input[["condition"]]]] %in% c(input[['variable1']], input[['variable2']])
      choices <- c("t-test", "Wilcox rank sum")
      updateSelectInput(session, 'pipeline', choices = choices)
    })


    DEresults <- reactive({
      print('hi')
      shinyjs::disable("goDE")
      pseudobulk.metadata = pseudobulk_samples()$metadata
      print(pseudobulk.metadata)
      pseudobulk.expression.matrix = pseudobulk_samples()$expression_matrix
      print(pseudobulk.expression.matrix)
      condition.indices <- pseudobulk.metadata[[input[["condition"]]]] %in% c(input[['variable1']], input[['variable2']])
      print(condition.indices)
      # Need to add error if any group has <2 samples
      DEtable <- DEanalysis(
        expression.matrix = pseudobulk.expression.matrix[, condition.indices],
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
      bindCache(utils::head(full.expression.matrix), full.metadata, input[["condition"]],
                input[['variable1']], input[['variable2']], input[["pipeline"]],
                input[["pvalThreshold"]],input[['lfcThreshold']]) %>%
      bindEvent(input[["goDE"]])

    #Define output table (only DE peaks)
    dataTable <- reactive({
      DEresults()$DEtableSubset %>%
        DT::datatable() %>%
        DT::formatSignif(columns = c('pval', 'pvalAdj','lfc','log2exp'), digits = 3)
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
      print(DEresults()$DEtableSubset$m_z[input$data_rows_selected])
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