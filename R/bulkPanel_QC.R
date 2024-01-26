#' @rdname BulkQCPanel
#' @export
BulkQCpanelUI <- function(id, bulk.metadata, show = TRUE){
  ns <- NS(id)

  if(show){
    tabPanel(
      'Quality checks',
      tags$h1("Principal Component Analysis (PCA)"),
      shinyWidgets::dropMenu(
        circleButton(ns("info_pca"), icon = icon("gear"),status = "success"),
        tags$h3("PCA"),
        tags$ul(
          tags$li("Distribution of samples across the first 2 principal components."),
          tags$li("Each sample is coloured by the selected sample-wide metadata information."),
          tags$li("Metadata groups are surrounded by minimal ellipses containing all samples (ggplot2::stat_ellipse) or 95% confidence ellipses as selected in the options below."),
          tags$li("Each sample can also be labelled.")
        ),
        br(),
        radioButtons(ns('pca.annotation'), label = "Group by",
                     choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[ncol(bulk.metadata)]),
        checkboxInput(ns("pca.show.labels"), label = "Show sample labels", value = FALSE),
        checkboxInput(ns('pca.show.ellipses'),label = "Show ellipses around groups",value=TRUE),
        checkboxInput(ns('pca.show.confidence.ellipses'),label = "Show 95% confidence ellipses around groups",value=FALSE),
        theme = "light-border",
        placement = "right",
        arrow = FALSE
      ),
      div(style = "margin-top:10px"),
      dropMenu(
        circleButton(ns("downloads_pca"), icon = icon("download"),status = "success"),
        tags$div(
          tags$h3("Downloads"),
          fluidRow(
            column(10,offset=0,
                   textInput(ns('plotPCAFileName'), 'File name for PCA plot download', value ='PCAPlot.png'),
                   numericInput(ns('PCAPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                   numericInput(ns('PCAPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                   downloadButton(ns('downloadPCAPlot'), 'Download PCA Plot'),
            )),
          theme = "light-border",
          placement = "right",
          arrow = FALSE
        )),
      plotOutput(ns('pca')),

      tags$h1("Partial Least Squares Discriminant Analysis (PLS-DA)"),
      shinyWidgets::dropMenu(
        circleButton(ns("info_plsda"), icon = icon("gear"),status = "success"),
        tags$h3("PLS-DA"),
        tags$ul(
          tags$li("Distribution of samples across the first 2 PLS-DA components as computed using mixOmics."),
          tags$li("The discriminating metadata (used to compute PLS-DA which maximises separation between these groups) and the metadata to colour samples by can both be selected below."),
          tags$li("As with PCA, metadata groups are surrounded by minimal ellipses containing all samples (ggplot2::stat_ellipse) or 95% confidence ellipses as selected in the options below."),
          tags$li("Each sample can also be labelled."),
          tags$li("The top peaks contributing to the PLS-DA are also shown below.")
        ),
        br(),

        radioButtons(ns('plsda.separator'), label = "Discriminating bulk.metadata",
                     choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[ncol(bulk.metadata)]),
        radioButtons(ns('plsda.annotation'), label = "Group by",
                     choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[ncol(bulk.metadata)]),
        checkboxInput(ns("plsda.show.labels"), label = "Show sample labels", value = FALSE),
        checkboxInput(ns('plsda.show.ellipses'),label = "Show ellipses around groups",value=TRUE),
        checkboxInput(ns('plsda.show.confidence.ellipses'),label = "Show 95% confidence ellipses around groups",value=FALSE),
        numericInput(ns('plsda.comp'),label = "PLS-DA component's contributions to show",min=1,max=2,step = 1,value = 1),
        theme = "light-border",
        placement = "right",
        arrow = FALSE
      ),
      div(style = "margin-top:10px"),
      dropMenu(
        circleButton(ns("downloads_plsda"), icon = icon("download"),status = "success"),
        tags$div(
          tags$h3("Downloads"),
          fluidRow(
            column(10,offset=0,
                   textInput(ns('plotPLSDAFileName'), 'File name for PLS-DA plot download', value ='PLSDAPlot.png'),
                   numericInput(ns('PLSDAPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                   numericInput(ns('PLSDAPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                   downloadButton(ns('downloadPLSDAPlot'), 'Download PLS-DA Plot'),
                   textInput(ns('plotPLSDAContribFileName'), 'File name for PLS-DA contribution plot download', value ='PLSDAContribPlot.png'),
                   numericInput(ns('PLSDAContribPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                   numericInput(ns('PLSDAContribPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                   downloadButton(ns('downloadPLSDAContribPlot'), 'Download PLS-DA contribution plot'),
            )),
          theme = "light-border",
          placement = "right",
          arrow = FALSE
        )),
      plotOutput(ns('plsda')),
      plotOutput(ns('plsda_contrib'),click = ns("plsda_hover")),
      tableOutput(ns("plsda_contrib_data")),

      tags$h1("Individual peak intensity barplots"),
      dropMenu(
        circleButton(ns("info_peak_barplot"), icon = icon("info"),status = "success"),
        tags$div(
          tags$h3("Peak intensity barplots"),
          tags$ul(
            tags$li("Select peaks and visualise their sample-wide average intensity compared to other samples in barplots."),
            tags$li("Choose sample-wide metadata information to colour each bar by."),
          )
        ),
        theme = "light-border",
        placement = "right",
        arrow = FALSE
      ),
      selectInput(ns("barPeakName"), "Peaks to include:", multiple = TRUE, choices = character(0)),
      radioButtons(ns('peak.barplot.colour'), label = "Group by",
                   choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[ncol(bulk.metadata)]),
      shinyWidgets::dropdownButton(
        textInput(ns('plotBarFileName'), 'File name for bar plot download', value ='BarPlot.png'),
        numericInput(ns('barPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
        numericInput(ns('barPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
        downloadButton(ns('downloadBarPlot'), 'Download Bar Plot'),

        status = "success",
        icon = icon("download", verify_fa = FALSE),
        tooltip = shinyWidgets::tooltipOptions(title = "Click to see information and options!")
      ),
      plotOutput(ns('barplot')),
      tags$h1("Individual peak intensity box plots"),
      dropMenu(
        circleButton(ns("info_peak_boxplot"), icon = icon("info"),status = "success"),
        tags$div(
          tags$h3("Peak intensity boxplots"),
          tags$ul(
            tags$li("Select peaks and visualise their sample-wide average intensity compared to other samples in boxplots, grouped by the selected sample-wide metadata information."),
          )
        ),
        theme = "light-border",
        placement = "right",
        arrow = FALSE
      ),

      selectInput(ns("boxPeakName"), "Peaks to include:", multiple = TRUE, choices = character(0)),
      radioButtons(ns('boxplot.metadata'), label = "Group by",
                   choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[ncol(bulk.metadata)]),
      shinyWidgets::dropdownButton(
        textInput(ns('plotBoxFileName'), 'File name for box plot download', value ='BoxPlot.png'),
        numericInput(ns('boxPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
        numericInput(ns('boxPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
        downloadButton(ns('downloadBoxPlot'), 'Download Box Plot'),

        status = "success",
        icon = icon("download", verify_fa = FALSE),
        tooltip = shinyWidgets::tooltipOptions(title = "Click to see information and options!")
      ),
      plotOutput(ns('boxplot'),click = ns('boxplot_click')),
      tableOutput(ns("box_data"))
    )
  }else{
    NULL
  }
}

#' @rdname BulkQCPanel
#' @export
BulkQCpanelServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno){
  ns <- NS(id)
  # check whether inputs (other than id) are reactive or not

  moduleServer(id, function(input, output, session){

    #Set up server-side search for peak names
    updateSelectizeInput(session, "barPeakName", choices = anno$display_name, server = TRUE, selected = anno$display_name[1:2])
    updateSelectizeInput(session, "boxPeakName", choices = anno$display_name, server = TRUE, selected = anno$display_name[1:2])

    #Set up server-side search for peak names
    updateSelectizeInput(session, "peakName", choices = anno$display_name, server = TRUE)

    pca.plot <- reactive({
      myplot <- plot_pca(
        intensity.matrix = bulk.intensity.matrix,
        metadata = bulk.metadata,
        annotation.id = match(input[['pca.annotation']], colnames(bulk.metadata)),
        n.abundant = nrow(bulk.intensity.matrix),
        show.labels = input[['pca.show.labels']],
        show.ellipses = input[['pca.show.ellipses']],
        show.confidence.ellipses = input[['pca.show.confidence.ellipses']]
      )
      myplot
    })
    output[['pca']] <- renderPlot(pca.plot())

    plsda.plot <- reactive({
      myplot <- plot_plsda(
        intensity.matrix = bulk.intensity.matrix,
        metadata = bulk.metadata,
        separator.id = match(input[['plsda.separator']], colnames(bulk.metadata)),
        annotation.id = match(input[['plsda.annotation']], colnames(bulk.metadata)),
        show.labels = input[['plsda.show.labels']],
        show.ellipses = input[['plsda.show.ellipses']],
        show.confidence.ellipses = input[['plsda.show.confidence.ellipses']],
      )
      myplot
    })
    output[['plsda']] <- renderPlot(plsda.plot())

    plsda.contrib <- reactive({
      myplot <- plsda_contrib(intensity.matrix=bulk.intensity.matrix,
                              metadata=bulk.metadata,
                              separator.id = match(input[['plsda.separator']], colnames(bulk.metadata)),
                              comp = input[['plsda.comp']],
                              anno)
      myplot
    })
    output[['plsda_contrib']] <- renderPlot(plsda.contrib()$plot)
    bar.plot <- reactive({
      peak.ids <- anno$m_z[match(input[["barPeakName"]],anno$display_name)]
      if (length(peak.ids)==1){
        sub.intensity.matrix <- data.frame(bulk.intensity.matrix[peak.ids,,drop=F])
      } else {
        sub.intensity.matrix <- data.frame(bulk.intensity.matrix[peak.ids,,drop=F])
      }
      rownames(sub.intensity.matrix) <- input[["barPeakName"]]
      myplot <- peaks_barplot(
        sub.intensity.matrix = sub.intensity.matrix,
        log.transformation = F,
        condition.vector = bulk.metadata[,input[['peak.barplot.colour']]])
      myplot
    })
    output[['barplot']] <- renderPlot(bar.plot())

    box.plot <- reactive({
      peak.ids <- anno$m_z[match(input[["boxPeakName"]],anno$display_name)]
      if (length(peak.ids)==1){
        sub.intensity.matrix <- data.frame(bulk.intensity.matrix[peak.ids,,drop=F])
      } else {
        sub.intensity.matrix <- data.frame(bulk.intensity.matrix[peak.ids,,drop=F])
      }
      rownames(sub.intensity.matrix) <- input[["boxPeakName"]]
      myplot <- peaks_boxplot(
        sub.intensity.matrix = sub.intensity.matrix,
        log.transformation = F,
        metadata = bulk.metadata,
        metadata.column = input[['boxplot.metadata']])
      return(myplot)
    })
    output[['boxplot']] <- renderPlot(box.plot()$plot)

    # output$data <- renderTable({
    #   nearPoints(box.plot()$table, input$boxplot_click)
    # })

    output$plsda_contrib_data <- renderTable({
      if (is.null(input$plsda_hover$y)) return()
      selected.metadata = unique(sort(plsda.contrib()$table$metab))[round(input$plsda_hover$y)]
      keeprows <- selected.metadata == plsda.contrib()$table$metab
      keeprows <- as.data.frame(plsda.contrib()$table[keeprows, ])
      keeprows <- keeprows[,c('PLSDAComp','metab','display_metab')]
      colnames(keeprows)=c('PLSDA Contribution','m_z','Metabolite annotation')
      keeprows
    })

    output$box_data <- renderTable({
      if (is.null(input$boxplot_click$x)) return()
      panel = input$boxplot_click$panelvar1
      selected.metadata = sort(unique(box.plot()$table$metadata))[round(input$boxplot_click$x)]
      keeprows <- selected.metadata == box.plot()$table$metadata & box.plot()$table$peak==panel
      keeprows <- as.data.frame(box.plot()$table[keeprows, ])
      keeprows$peak_name <- anno$name[match(keeprows$peak,anno$m_z)]
      keeprows = keeprows[order(abs(input$boxplot_click$y-keeprows$value)),]
      head(keeprows,5)
    })

    output[['downloadPCAPlot']] <- downloadHandler(
      filename = function() { input[['plotPCAFileName']] },
      content = function(file) {
        ggsave(file, plot = pca.plot(), width=input[['PCAPlotWidth']],height=input[['PCAPlotHeight']],units='in')
      }
    )

    output[['downloadPLSDAPlot']] <- downloadHandler(
      filename = function() { input[['plotPLSDAFileName']] },
      content = function(file) {
        ggsave(file, plot = plsda.plot(), , width=input[['PLSDAPlotWidth']],height=input[['PLSDAPlotHeight']],units='in')
      }
    )

    output[['downloadPLSDAContribPlot']] <- downloadHandler(
      filename = function() { input[['plotPLSDAContribFileName']] },
      content = function(file) {
        ggsave(file, plot = plsda.contrib()$plot, width=input[['PLSDAContribPlotWidth']],height=input[['PLSDAContribPlotHeight']],units='in')
      }
    )

    output[['downloadBarPlot']] <- downloadHandler(
      filename = function() { input[['plotBarFileName']] },
      content = function(file) {
        ggsave(file, plot = bar.plot() , width=input[['barPlotWidth']],height=input[['barPlotHeight']],units='in')
      }
    )

    output[['downloadBoxPlot']] <- downloadHandler(
      filename = function() { input[['plotBoxFileName']] },
      content = function(file) {
        ggsave(file, plot = box.plot()$plot, , width=input[['boxPlotWidth']],height=input[['boxPlotHeight']],units='in')
      }
    )


  })
}

# QCpanelApp <- function(){
#   shinyApp(
#     ui = fluidPage(QCpanelMetabUI('qc', bulk.metadata)),
#     server = function(input, output, session){
#       QCpanelMetabServer('qc', bulk.intensity.matrix[[1]], bulk.metadata[[1]])
#     }
#   )
# }
