#' Visualises bulk-level quality control
#'
#' @description UI and server logic for quality control and exploratory analysis of bulk mass spectrometry data. Provides multiple perspectives on sample relationships and peak distributions.
#'
#' @details
#' \itemize{
#'   \item{Principal Component Analysis (PCA) with metadata annotation}
#'   \item{Partial Least Squares Discriminant Analysis (PLS-DA) for supervised analysis}
#'   \item{Individual peak intensity bar plots and box plots}
#'   \item{Interactive visualization with downloadable plots}
#'   \item{The module performs dimension reduction on bulk peak intensity data, visualizing sample relationships and peak contributions. Supports filtering by metadata variables and includes confidence ellipses for grouped samples. All plots are interactive (plotly) with customizable download options.}
#' }
#' @name BulkPanel_QCTab
#' @rdname BulkPanel_QCTab
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame with bulk sample metadata
#' @param show Logical; whether to render the panel (default: TRUE)
#' @return A shiny::tabPanel containing the UI elements for the QC panel
#' @export
BulkPanel_QCTabUI <- function(id, bulk.metadata, show = TRUE){
  ns <- shiny::NS(id)

  if(show){
    shiny::tabPanel(
      'Quality Checks',
      shiny::tags$h1("Principal Component Analysis (PCA)"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_pca"), icon = shiny::icon("info-circle"),status = "info"),
              shiny::tags$h3("PCA"),
              shiny::tags$ul(
                shiny::tags$li("The PCA plot shows the distribution of samples across the first two principal components, colored by selected metadata."),
                shiny::tags$li("Confidence ellipses (if enabled) highlight group separation."),
                shiny::tags$li("The PCA contributions plot shows the top peaks contributing to each component."),
                shiny::tags$li("Both plots are interactive and downloadable.")
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_pca"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download PCA plots"),
                shiny::textInput(ns('plotPCAFileName'), 'File name for PCA plot download', value ='PCAPlot.png'),
                shiny::numericInput(ns('PCAPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('PCAPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadPCAPlot'), 'Download PCA Plot'),
                shiny::tags$hr(),
                shiny::textInput(ns('plotPCAContribFileName'), 'File name for PCA contribution plot download', value ='PCAContribPlot.png'),
                shiny::numericInput(ns('PCAContribPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('PCAContribPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadPCAContribPlot'), 'Download PCA contribution plot')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
            shiny::radioButtons(ns('pca.annotation'), label = "Group by",
              choices = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)]), selected = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)])[1]),
          shiny::checkboxInput(ns('pca.show.confidence.ellipses'),label = "Show 95% confidence ellipses around groups",value=TRUE),
          shiny::numericInput(ns('pca.comp'),label = "PCA's contributions to show",min=1,max=2,step = 1,value = 1)
        ),
        shiny::mainPanel(
          plotly::plotlyOutput(ns('pca'),width='800px'),
          plotly::plotlyOutput(ns('pca_contrib'),width='800px')
        )
      ),

      shiny::tags$h1("Partial Least Squares Discriminant Analysis (PLS-DA)"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_plsda"), icon = shiny::icon("info-circle"),status = "info"),
              shiny::tags$h3("PLS-DA"),
              shiny::tags$ul(
                shiny::tags$li("The PLS-DA plot shows sample separation based on selected discriminating metadata (first user selection)."),
                shiny::tags$li("Confidence ellipses (if enabled) highlight group separation based on the second user selection."),
                shiny::tags$li("The PLS-DA contributions plot shows the top peaks contributing to each component."),
                shiny::tags$li("Both plots are interactive and downloadable.")
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_plsda"), icon = shiny::icon("download"), status = "success"),
              shiny::div(
                shiny::tags$h3("Download PLS-DA plots"),
                shiny::textInput(ns('plotPLSDAFileName'), 'File name for PLS-DA plot download', value ='PLSDAPlot.png'),
                shiny::numericInput(ns('PLSDAPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('PLSDAPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadPLSDAPlot'), 'Download PLS-DA Plot'),
                shiny::tags$hr(),
                shiny::textInput(ns('plotPLSDAContribFileName'), 'File name for PLS-DA contribution plot download', value ='PLSDAContribPlot.png'),
                shiny::numericInput(ns('PLSDAContribPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('PLSDAContribPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadPLSDAContribPlot'), 'Download PLS-DA contribution plot')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
            shiny::radioButtons(ns('plsda.separator'), label = "Condition to run PLS-DA on",
              choices = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)]), selected = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)])[1]),
            shiny::radioButtons(ns('plsda.annotation'), label = "Condition to color PLS-DA plot by",
                  choices = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)]), selected = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)])[1]),
          shiny::checkboxInput(ns('plsda.show.confidence.ellipses'),label = "Show 95% confidence ellipses around groups",value=TRUE),
          shiny::numericInput(ns('plsda.comp'),label = "PLS-DA component's contributions to show",min=1,max=2,step = 1,value = 1),

        ),
        shiny::mainPanel(
          plotly::plotlyOutput(ns('plsda'),width='800px'),
          plotly::plotlyOutput(ns('plsda_contrib'),width='800px'),
          shiny::tableOutput(ns("plsda_contrib_data"))
        )
      ),

      shiny::tags$h1("Individual peak intensity barplots"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_peak_barplot"), icon = shiny::icon("info-circle"),status = "info"),
              shiny::tags$h3("Barplot Information"),
              shiny::tags$ul(
                shiny::tags$li("The barplot shows average intensity of selected peaks across samples, colored by metadata group."),
                shiny::tags$li("Use the settings to select peaks and grouping variable.")
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_barplot"), icon = shiny::icon("download"), status = "success"),
              shiny::div(
                shiny::tags$h3("Download Bar Plot"),
                shiny::textInput(ns('plotBarFileName'), 'File name for bar plot download', value ='BarPlot.png'),
                shiny::numericInput(ns('barPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('barPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadBarPlot'), 'Download Bar Plot')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::selectInput(ns("barPeakName"), "Peaks to include:", multiple = TRUE, choices = character(0)),
          shiny::radioButtons(ns('peak.barplot.colour'), label = "Group by",
                      choices = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)]), selected = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)])[1]),

        ),
        shiny::mainPanel(
          plotly::plotlyOutput(ns('barplot'))
        )
      ),
      shiny::tags$h1("Individual peak intensity box plots"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("info_peak_boxplot"), icon = shiny::icon("info-circle"),status = "info"),
              shiny::tags$h3("Boxplot Information"),
              shiny::tags$ul(
                shiny::tags$li("The boxplot shows intensity distribution of selected peaks across samples, grouped by metadata."),
                shiny::tags$li("Use the settings to select peaks and grouping variable.")
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            ),
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_boxplot"), icon = shiny::icon("download"), status = "success"),
              shiny::div(
                shiny::tags$h3("Download Box Plot"),
                shiny::textInput(ns('plotBoxFileName'), 'File name for box plot download', value ='BoxPlot.png'),
                shiny::numericInput(ns('boxPlotWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
                shiny::numericInput(ns('boxPlotHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
                shiny::downloadButton(ns('downloadBoxPlot'), 'Download Box Plot')
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          shiny::selectInput(ns("boxPeakName"), "Peaks to include:", multiple = TRUE, choices = character(0)),
          shiny::radioButtons(ns('boxplot.metadata'), label = "Group by",
                      choices = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)]), selected = colnames(bulk.metadata[,sapply(bulk.metadata,dplyr::n_distinct)!=nrow(bulk.metadata)])[1]),

        ),
        shiny::mainPanel(
          plotly::plotlyOutput(ns('boxplot')),
          shiny::tableOutput(ns("box_data"))
        )
      )
    )
  }else{
    NULL
  }
}

#' @rdname BulkPanel_QCTab
#' @param bulk.intensity.matrix Numeric matrix of bulk sample intensities (features × samples)
#' @param anno Data frame with peak annotation (must include 'display_name' and 'm_z')
#' @export
BulkPanel_QCTabServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno){
  ns <- shiny::NS(id)
  # check whether inputs (other than id) are reactive or not

  shiny::moduleServer(id, function(input, output, session){

    # Peak selectize setup via helper
    shiny::req(!is.null(anno), !is.null(anno$display_name))
    utils_update_peak_selectize(session, "barPeakName", choices = anno$display_name, selected_n = 2)
    utils_update_peak_selectize(session, "boxPeakName", choices = anno$display_name, selected_n = 2)
    utils_update_peak_selectize(session, "peakName", choices = anno$display_name, selected_n = 0)

    run_pca <- shiny::reactive({
      expr.PCA.list <- bulk.intensity.matrix |>
        as.data.frame() |>
        t()
      expr.PCA.res <- expr.PCA.list[, apply(expr.PCA.list, 2, function(x) max(x) != min(x))] |>
        stats::prcomp(center = TRUE, scale = TRUE)
      return(expr.PCA.res)
    })
    pca.plot <- shiny::reactive({
      myplot <- bulk_utils_plot_pca(
        run_pca(),
        metadata = bulk.metadata,
        annotation.id = match(input[['pca.annotation']], colnames(bulk.metadata)),
        show.confidence.ellipses = input[['pca.show.confidence.ellipses']]
      )
      myplot$plot
    })
    output[['pca']] <- plotly::renderPlotly(plotly::ggplotly(pca.plot()))

    pca.contrib <- shiny::reactive({
      myplot <- bulk_utils_pca_contrib(
        run_pca(),
        comp = input[['pca.comp']],
        anno = anno
      )
      myplot$plot
    })
    output[['pca_contrib']] <- plotly::renderPlotly(plotly::ggplotly(pca.contrib()))

    plsda.plot <- shiny::reactive({
      myplot <- bulk_utils_plot_plsda(
        intensity.matrix = bulk.intensity.matrix,
        metadata = bulk.metadata,
        separator.id = match(input[['plsda.separator']], colnames(bulk.metadata)),
        annotation.id = match(input[['plsda.annotation']], colnames(bulk.metadata)),
        show.confidence.ellipses = input[['plsda.show.confidence.ellipses']],
      )
      myplot
    })
    output[['plsda']] <- plotly::renderPlotly(plotly::ggplotly(plsda.plot()))

    plsda.contrib <- shiny::reactive({
      myplot <- bulk_utils_plsda_contrib(intensity.matrix=bulk.intensity.matrix,
                              metadata=bulk.metadata,
                              separator.id = match(input[['plsda.separator']], colnames(bulk.metadata)),
                              comp = input[['plsda.comp']],
                              anno)
      myplot
    })
    output[['plsda_contrib']] <- plotly::renderPlotly(plotly::ggplotly(plsda.contrib()$plot))
    bar.plot <- shiny::reactive({
      peak.ids <- anno$m_z[match(input[["barPeakName"]],anno$display_name)]
      if (length(peak.ids)==1){
        sub.intensity.matrix <- data.frame(bulk.intensity.matrix[peak.ids,,drop=F])
      } else {
        sub.intensity.matrix <- data.frame(bulk.intensity.matrix[peak.ids,,drop=F])
      }
      rownames(sub.intensity.matrix) <- input[["barPeakName"]]
      myplot <- bulk_utils_peaks_barplot(
        sub.intensity.matrix = sub.intensity.matrix,
        log.transformation = F,
        condition.vector = bulk.metadata[,input[['peak.barplot.colour']]])
      myplot
    })
    output[['barplot']] <- plotly::renderPlotly(plotly::ggplotly(bar.plot()))

    box.plot <- shiny::reactive({
      peak.ids <- anno$m_z[match(input[["boxPeakName"]],anno$display_name)]
      if (length(peak.ids)==1){
        sub.intensity.matrix <- data.frame(bulk.intensity.matrix[peak.ids,,drop=F])
      } else {
        sub.intensity.matrix <- data.frame(bulk.intensity.matrix[peak.ids,,drop=F])
      }
      rownames(sub.intensity.matrix) <- input[["boxPeakName"]]
      myplot <- bulk_utils_peaks_boxplot(
        sub.intensity.matrix = sub.intensity.matrix,
        log.transformation = F,
        metadata = bulk.metadata,
        metadata.column = input[['boxplot.metadata']])
      return(myplot)
    })
    output[['boxplot']] <- plotly::renderPlotly(plotly::ggplotly(box.plot()$plot))

    output$plsda_contrib_data <- shiny::renderTable({
      if (is.null(input$plsda_hover$y)) return()
      selected.metadata = unique(sort(plsda.contrib()$table$metab))[round(input$plsda_hover$y)]
      keeprows <- selected.metadata == plsda.contrib()$table$metab
      keeprows <- as.data.frame(plsda.contrib()$table[keeprows, ])
      keeprows <- keeprows[,c('PLSDAComp','metab','display_metab')]
      colnames(keeprows)=c('PLSDA Contribution','m_z','Metabolite annotation')
      keeprows
    })

    output$box_data <- shiny::renderTable({
      if (is.null(input$boxplot_click$x)) return()
      panel = input$boxplot_click$panelvar1
      selected.metadata = sort(unique(box.plot()$table$metadata))[round(input$boxplot_click$x)]
      keeprows <- selected.metadata == box.plot()$table$metadata & box.plot()$table$peak==panel
      keeprows <- as.data.frame(box.plot()$table[keeprows, ])
      keeprows$peak_name <- anno$name[match(keeprows$peak,anno$m_z)]
      keeprows = keeprows[order(abs(input$boxplot_click$y-keeprows$value)),]
      utils::head(keeprows,5)
    })

    output[['downloadPCAPlot']] <- utils_create_download_plot_handler(
      plot_func = pca.plot,
      filename_func = function() input[['plotPCAFileName']],
      width_func = function() input[['PCAPlotWidth']],
      height_func = function() input[['PCAPlotHeight']]
    )

    output[['downloadPCAContribPlot']] <- utils_create_download_plot_handler(
      plot_func = pca.contrib,
      filename_func = function() input[['plotPCAContribFileName']],
      width_func = function() input[['PCAContribPlotWidth']],
      height_func = function() input[['PCAContribPlotHeight']]
    )

    output[['downloadPLSDAPlot']] <- utils_create_download_plot_handler(
      plot_func = plsda.plot,
      filename_func = function() input[['plotPLSDAFileName']],
      width_func = function() input[['PLSDAPlotWidth']],
      height_func = function() input[['PLSDAPlotHeight']]
    )

    output[['downloadPLSDAContribPlot']] <- utils_create_download_plot_handler(
      plot_func = function() plsda.contrib()$plot,
      filename_func = function() input[['plotPLSDAContribFileName']],
      width_func = function() input[['PLSDAContribPlotWidth']],
      height_func = function() input[['PLSDAContribPlotHeight']]
    )

    output[['downloadBarPlot']] <- utils_create_download_plot_handler(
      plot_func = bar.plot,
      filename_func = function() input[['plotBarFileName']],
      width_func = function() input[['barPlotWidth']],
      height_func = function() input[['barPlotHeight']]
    )

    output[['downloadBoxPlot']] <- utils_create_download_plot_handler(
      plot_func = function() box.plot()$plot,
      filename_func = function() input[['plotBoxFileName']],
      width_func = function() input[['boxPlotWidth']],
      height_func = function() input[['boxPlotHeight']]
    )


  })
}