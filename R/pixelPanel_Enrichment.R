#' @rdname PixelEnrichmentPanel
#' @export
PixelEnrichmentPanelUI <- function(id, bulk.metadata, full.metadata, pixel_enrichment, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Spatial enrichment',

      selectInput(
        inputId = ns("samplesToShow"),
        label = "Select samples to show",
        choices = unique(bulk.metadata[,1]),
        selected = unique(bulk.metadata[,1])[1],
        multiple = TRUE,
      ),
      actionButton(ns("fix_samples"),'Fix samples and sort pathways by significance'),
      sidebarLayout(
        sidebarPanel(
          dropMenu(
            circleButton(ns("info_peak"), icon = icon("info"),status = "success"),
            tags$div(
              tags$h3("Spatial enrichment visualisation"),
              tags$ul(
                tags$li("Select a pathway and visualise its p-value across the samples selected at the top of this panel."),
                tags$li("The adjusted p-value can be shown directly or the -log10 p-value"),
                tags$li("Density plots and barplots showing the distribution of p-values and proportion of significant pixels are also shown"),
                tags$li("Caps can also be applied to the colour scale to aid visualisation.")
              )
            ),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),

          selectInput(ns("pathwayName"), "Pathway to show:", multiple = FALSE, choices = c()),
          selectInput(ns("pvalShown"),"How to show p-value:",choices = c('Adjusted p-value','-log10(adjusted p-value)')),
          numericInput(ns('sigThreshold'), 'Significance threshold', value = 0.05, min = 0.001, max = 1,step = 0.001),
          conditionalPanel(
            id = ns('pvalCap_log'),
            ns=ns,
            condition = "input.pvalShown == '-log10(adjusted p-value)'",
            numericInput(ns('logCap'), 'Upper cap on log scale', value = 10, min = 1, max = 1000,step = 1),
          ),
          radioButtons(ns('barplot_metadata'), label = "Color bar plot by",
                       choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[ncol(bulk.metadata)]),
          div(style = "margin-top:10px"),
          dropMenu(
            circleButton(ns("downloadsMetadata"), icon = icon("download"),status = "success"),
            tags$div(
              tags$h3("Downloads"),
              fluidRow(
                column(5,offset=0,
                       tags$h4("Spatial distribution"),
                       textInput(ns('spatialFileName'),'File name for download', value ='spatial.png', placeholder = 'spatial.png'),
                       numericInput(ns('spatialWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('spatialHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadSpatial'), 'Download spatial figure')),
                column(5,offset=1,
                       tags$h4("Proportion of significant pixels plot"),
                       textInput(ns('barFileName'),'File name for download', value ='bar.png', placeholder = 'bar.png'),
                       numericInput(ns('barWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('barHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadBar'), 'Download bar plot')),

              ),
              fluidRow(
                column(10,offset=0,
                       tags$h4("Top metabolite plot"),
                       textInput(ns('topMetFileName'),'File name for download', value ='topMet.png', placeholder = 'topMet.png'),
                       numericInput(ns('topMetWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                       numericInput(ns('topMetHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                       downloadButton(ns('downloadTopMet'), 'Download top metabolites')))),
            theme = "light-border",
            placement = "right",
            arrow = FALSE
          ),
        ),
        mainPanel(
#          plotOutput(ns('peakDensity')),
          fluidRow(column=10,plotOutput(ns('plotPeak'),click = ns('peak_click'))),
          plotOutput(ns('plotBar')),
          plotOutput(ns('plotTopMetab')))),
      )
  }else{
    NULL
  }
}

#' @rdname PixelEnrichmentPanel
#' @export
PixelEnrichmentPanelServer <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno, pixel_enrichment){

  moduleServer(id, function(input, output, session){

    ordered_pathway_list <- reactive({
      pathway_table = pixel_enrichment |> filter(sample %in% input[['samplesToShow']])
      pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
      pathway_table = pathway_table |> dplyr::group_by(pathway) %>% summarise(non_na_count = sum(!is.na(pathwaySig))) |>
        dplyr::arrange(desc(non_na_count), .by_group = FALSE) |> filter(non_na_count!=0) |>
        dplyr::select(pathway)
      print(head(pathway_table))
      return(pathway_table$pathway)
    })  %>%
      bindEvent(input[['fix_samples']])

    observe({
    updateSelectizeInput(session, "pathwayName", choices = ordered_pathway_list(), server = TRUE, selected = ordered_pathway_list()[1])
    })

    selected_pathway_table <- reactive({
      pathway_table = pixel_enrichment |> filter(sample %in% input[['samplesToShow']])
#      pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
      return(pathway_table)
    }) %>%
      bindEvent(input[['fix_samples']])

    selected.samples.metadata <- reactive({
      print(input[['samplesToShow']])
      print(unique(full.metadata$Group))
      print(nrow(full.metadata |> filter(Group %in% input[['samplesToShow']])))
      return(full.metadata |> filter(Group %in% input[['samplesToShow']]))
    }) %>%
      bindEvent(input[['fix_samples']])

    show_peak <- reactive({
      pathway_table = selected_pathway_table() |> filter(pathway==input[['pathwayName']])
      pathway_table = merge(selected.samples.metadata(),pathway_table,all.x=T)
      # pathway_table = pixel_enrichment |> filter(pathway==input[['pathwayName']]) |> filter(sample %in% input[['samplesToShow']])
      # print(nrow(pathway_table))
      # pathway_table = merge(metadata |> filter(Group %in% input[['samplesToShow']]),pathway_table,all.x=T)
      # print(nrow(pathway_table))
      # print(nrow(metadata))
      pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
      if (input[['pvalShown']]=='-log10(adjusted p-value)'){
        pathway_table$pathwaySig = -log10(pathway_table$pathwaySig)
        pathway_table$pathwaySig = ifelse(pathway_table$pathwaySig>input[['logCap']],input[['logCap']],pathway_table$pathwaySig)
      }
      legend_title=input[['pvalShown']]
      print(head(pathway_table))
      spatial.plot = ggplot2::ggplot(pathway_table,ggplot2::aes(x=x,y=y,color=pathwaySig,fill=pathwaySig))+geom_tile()+
        ggplot2::facet_wrap(~pathway_table$Group, scales = 'free',nrow=max(1,floor(sqrt(length(input[['samplesToShow']])/2))))  +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(legend_title),
                       aspect.ratio = 1)
      if (input[['pvalShown']]=='-log10(adjusted p-value)'){
        spatial.plot <- spatial.plot +
          ggplot2::scale_fill_gradient(name=legend_title,low = "lightgrey", high = "brown",na.value = 'lightgrey')+
          ggplot2::scale_color_gradient(name=legend_title,low = "lightgrey", high = "brown",na.value = 'lightgrey')
      } else {
        spatial.plot <- spatial.plot +
        ggplot2::scale_fill_gradient(name=legend_title,high = "lightgrey", low = "brown",na.value = 'lightgrey')+
        ggplot2::scale_color_gradient(name=legend_title,high = "lightgrey", low = "brown",na.value = 'lightgrey')}

         density.plot = ggplot2::ggplot(pathway_table,ggplot2::aes(x=FDR,color=Group))+geom_density()+theme_classic()
      return(list('spatial'=spatial.plot,'density'=density.plot))
    }) %>%
      bindCache(input[['pathwayName']],input[['pvalShown']],input[['logCap']],input[['sigThreshold']],input[['fix_samples']])


    prop_significant <- reactive({
      if (!is.null(input[['pathwayName']])){
      # pathway_table = pixel_enrichment |> filter(pathway==input[['pathwayName']]) |> filter(sample %in% input[['samplesToShow']])
      # print(nrow(pathway_table))
      # pathway_table = merge(metadata |> filter(Group %in% input[['samplesToShow']]),pathway_table,all.x=T)
      # print(nrow(pathway_table))
      pathway_table = selected_pathway_table() |> filter(pathway==input[['pathwayName']])
    # NEED TO WORK OUT WHY THE NUMBER OF PIXELS DOESNT MATCH!!!!
      pathway_table = merge(selected.samples.metadata(),pathway_table,all.x=T)
      pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
      pathway_table = pathway_table |> dplyr::group_by(dplyr::across(dplyr::all_of(c('Group',input[['barplot_metadata']])))) %>% summarise(non_na_count = sum(!is.na(pathwaySig)),count = n()) |>
        dplyr::arrange(desc(non_na_count), .by_group = FALSE) |> filter(non_na_count!=0)

      ggplot(pathway_table,aes(x=Group,y=non_na_count/count,fill=get(input[['barplot_metadata']])))+
        geom_bar(stat='identity')+
        ylab('Proportion of significant pixels')+
        xlab('Sample')+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        scale_fill_discrete(name=input[['barplot_metadata']])
      }
    }) #%>%
      #bindEvent(input[['fix_samples']])

    show_peak_zoom <- reactive({

      pathway_table = pixel_enrichment |> filter(pathway==input[['pathwayName']]) |> filter(sample == input$peak_click$panelvar1)
      pathway_table = merge(full.metadata |> filter(Group == input$peak_click$panelvar1),pathway_table,all.x=T)

      pathway_table$pathwaySig = ifelse(pathway_table$FDR>input[['sigThreshold']],NA,pathway_table$FDR)
      if (input[['pvalShown']]=='-log10(adjusted p-value)'){
        pathway_table$pathwaySig = -log10(pathway_table$pathwaySig)
        pathway_table$pathwaySig = ifelse(pathway_table$pathwaySig>input[['logCap']],input[['logCap']],pathway_table$pathwaySig)
      }
      legend_title=input[['pvalShown']]
      spatial.plot = ggplot2::ggplot(pathway_table,ggplot2::aes(x=x,y=y,color=pathwaySig,fill=pathwaySig))+geom_tile()+
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks.y=ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(legend_title),
                       aspect.ratio = 1)
      if (input[['pvalShown']]=='-log10(adjusted p-value)'){
        spatial.plot <- spatial.plot +
          ggplot2::scale_fill_gradient(name=legend_title,low = "lightgrey", high = "brown",na.value = 'lightgrey')+
          ggplot2::scale_color_gradient(name=legend_title,low = "lightgrey", high = "brown",na.value = 'lightgrey')
      } else {
        spatial.plot <- spatial.plot +
          ggplot2::scale_fill_gradient(name=legend_title,high = "lightgrey", low = "brown",na.value = 'lightgrey')+
          ggplot2::scale_color_gradient(name=legend_title,high = "lightgrey", low = "brown",na.value = 'lightgrey')}
     return(spatial.plot)
    })

    top_metabolites <- reactive({
      pathway_table = selected_pathway_table() |>
        filter(pathway==input[['pathwayName']]) |>
        filter(FDR<input[['sigThreshold']])
      full.metabolite.list = unlist(strsplit(pathway_table$metabolites,split = '; '))
      print(head(full.metabolite.list))
      full.metabolite.freq = data.frame(table(full.metabolite.list))
      colnames(full.metabolite.freq)=c('kegg','freq')
      print(head(full.metabolite.freq))
      print(nrow(pathway_table))
      full.metabolite.freq$prop = full.metabolite.freq$freq/nrow(pathway_table)
      print(summary(full.metabolite.freq$prop))
      full.metabolite.freq = head(full.metabolite.freq[order(-full.metabolite.freq$prop),],30)
      full.metabolite.freq$kegg = factor(full.metabolite.freq$kegg,levels=rev(full.metabolite.freq$kegg))
      ggplot(full.metabolite.freq,aes(x=prop,y=kegg))+geom_bar(stat='identity')+theme_classic()
    })
    #%>% bindEvent(input[["go_plot_peak"]])

    output[['plotBar']] <- renderPlot({
      prop_significant()})
    #%>% bindEvent(input[["go_plot_metadata"]])
    output[['plotPeak']] <- renderPlot({
      show_peak()$spatial})

    output[['plotTopMetab']] <- renderPlot({
      top_metabolites()})

    output[['plotPeakZoom']] <- renderPlot({
      show_peak_zoom()})

    output[['peakDensity']] <- renderPlot({
      show_peak()$density
    })

    output[['downloadSpatial']] <- downloadHandler(
      filename = function() { input[['spatialFileName']] },
      content = function(file) {
        ggsave(file, plot = show_peak()$spatial, dpi = 300,
               width=input[['spatialWidth']],height=input[['spatialHeight']])
      }
    )

    output[['downloadBar']] <- downloadHandler(
      filename = function() { input[['barFileName']] },
      content = function(file) {
        ggsave(file, plot = prop_significant(), dpi = 300,
               width=input[['barWidth']],height=input[['barHeight']])
      }
    )

    output[['downloadTopMet']] <- downloadHandler(
      filename = function() { input[['topMetFileName']] },
      content = function(file) {
        ggsave(file, plot = top_metabolites(), dpi = 300,
               width=input[['topMetWidth']],height=input[['topMetHeight']])
      }
    )


    observeEvent(input$peak_click, {
      ns <- session$ns
      showModal(
        modalDialog(
          plotOutput(ns('plotPeakZoom')),
          easyClose = TRUE,
          footer = NULL
        )
      )
    })


  })
}
