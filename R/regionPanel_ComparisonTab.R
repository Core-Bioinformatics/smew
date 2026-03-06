#' Compares multiple regions and/or clusters
#'
#' @description UI and server logic for comparing clustering results between different methods or parameters in the SMEW app. Includes contingency table visualization, similarity indices, proportional composition, and spatial concordance plots.
#'
#' @details
#' \itemize{
#'  \item{Contingency table visualization (heatmap of cluster overlap)}
#'  \item{Jaccard similarity indices between clusters}
#'  \item{Proportional composition analysis across cluster assignments}
#'  \item{Spatial visualization of cluster concordance}
#'  \item{Download handlers for all plots}
#' }
#' @param id Shiny module id
#' @param shared_data Reactive or shared data object (server only)
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param show Logical; whether to render the panel (default: TRUE)
#' @return UI: A shiny::tabPanel object for the comparison panel. Server: None; called for side effects in Shiny module.
#'
#' @name RegionPanel_ComparisonTab
#' @rdname RegionPanel_ComparisonTab
#' @export
RegionPanel_ComparisonTabUI <- function(id, show = TRUE) {
  ns <- shiny::NS(id)
  if (show){
  shiny::tabPanel(
    'Region Comparison',
    bslib::accordion(
      bslib::accordion_panel(
        title = "Information",
        icon = bsicons::bs_icon("arrow-right-circle"),
        shiny::tags$ul(
          shiny::tags$li("Compare clustering results from different methods and/or parameters using contingency tables and similarity metrics (Jaccard similarity index, number of overlapping pixels, element-centric similarity (ECS)."),
          shiny::tags$li("Visualise cluster overlap as heatmaps and proportional barplots."),
          shiny::tags$li("Assess spatial concordance of cluster assignments across samples."),
          shiny::tags$li("Download all plots for reporting or further analysis.")
        )
      ),
      bslib::accordion_panel(
        title = "Downloads",
        icon = bsicons::bs_icon("download"),
        shiny::tags$div(
          shiny::tags$h4("Download output plots"),
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$strong("Overlap plot"),
              shiny::textInput(ns('overlapFileName'), 'File name', value ='OverlapCounts.png'),
              shiny::numericInput(ns('overlapWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('overlapHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadOverlap'), 'Download overlap counts')
            ),
            shiny::column(6,
              shiny::tags$strong("Jaccard heatmap"),
              shiny::textInput(ns('overlapJSIFileName'), 'File name', value ='OverlapJSI.png'),
              shiny::numericInput(ns('overlapJSIWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('overlapJSIHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadOverlapJSI'), 'Download Jaccard heatmap')
            )
          ),
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$strong("Proportions plot"),
              shiny::textInput(ns('proportionFileName'), 'File name', value ='Proportions.png'),
              shiny::numericInput(ns('proportionWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('proportionHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadProportion'), 'Download proportions')
            ),
            shiny::column(6,
              shiny::tags$strong("Spatial similarity plot"),
              shiny::textInput(ns('ecsSpatialFileName'), 'File name', value ='SpatialSimilarity.png'),
              shiny::numericInput(ns('ecsSpatialWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('ecsSpatialHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadECSSpatial'), 'Download spatial similarity')
            )
          ),
          shiny::fluidRow(
            shiny::column(6,
              shiny::tags$strong("Cluster/Metadata Proportion plot"),
              shiny::textInput(ns('clusterMetaProportionFileName'), 'File name', value ='ClusterMetaProportion.png'),
              shiny::numericInput(ns('clusterMetaProportionWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('clusterMetaProportionHeight'),value = 6,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadClusterMetaProportion'), 'Download cluster/meta proportion')
            )
          )
        )
      ),
      open = NULL
    ),
    shiny::tags$h4("Compare region proportions per sample across sample-level metadata"),
    # --- New section: Cluster vs sample-level metadata ---
    bslib::card(
      class = "mb-3",
      style = "width:100%;",

      shiny::fluidRow(
        shiny::column(4,
          shiny::selectInput(ns("clusterCol"), "Select clustering column:", choices = c())
        ),
        shiny::column(4,
          shiny::selectInput(ns("sampleMetaCol"), "Select sample-level metadata column:", choices = c())
        ),
        shiny::column(4,
          shiny::checkboxInput(ns("clusterMetaIncludeNAs"),label = "Include NA values?",value = F)
        )
      )
    ),
      shiny::plotOutput(ns("clusterMetaProportionPlot")),
    shiny::tags$h4("Compare multiple regions"),
    # --- Existing comparison UI ---
    bslib::card(
      class = "mb-3",
      style = "width:100%;",
      shiny::fluidRow(
        shiny::column(4,
          shiny::selectInput(ns("col1"), "Select first condition:", choices = c())
        ),
        shiny::column(4,
          shiny::selectInput(ns("col2"), "Select second condition:", choices = c())
        ),
        shiny::column(4,
          shiny::checkboxInput(ns("includeNAs"), label = "Include NA values?",value = F),
          shiny::checkboxInput(ns("scaleByNumber"), label = "Scale bar plot proportions by size of condition 2?", value = F)
        )
      )),
    shiny::fluidRow(
      shiny::column(4,
        shiny::tags$h4("Cluster Overlap Counts"),
        shiny::plotOutput(ns("overlapPlot"))
      ),
      shiny::column(4,
        shiny::tags$h4("Jaccard Similarity Index (JSI) Heatmap"),
        shiny::plotOutput(ns("overlapJSIPlot"))
      ),
      shiny::column(4,
        shiny::tags$h4("Cluster Proportions Barplot"),
        shiny::plotOutput(ns("proportionPlot"))
      )
    ),
    shiny::tags$h4("Spatial Similarity (ECS) Plot"),
    shiny::plotOutput(ns("ecsSpatial"),height = '600px')
  )
  } else {
    NULL
  }
}

#' @rdname RegionPanel_ComparisonTab
#' @param id Shiny module id
#' @param shared_data Reactive or shared data object
#' @param bulk.metadata Data frame of bulk sample metadata
#' @return None; called for side effects in Shiny module
#' @export
RegionPanel_ComparisonTabServer <- function(id, shared_data, bulk.metadata) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observe({
      categorical_cols <- colnames(shared_data$updated.metadata)[sapply(shared_data$updated.metadata, function(col) {
        n_unique <- length(unique(stats::na.omit(col)))
        n_unique >= 1 && n_unique < 20  # exclude constants & numeric-like columns
      })]
      shiny::updateSelectInput(session, 'col1', choices = categorical_cols)
      shiny::updateSelectInput(session, 'col2', choices = categorical_cols, selected = categorical_cols[2])
      shiny::updateSelectInput(session, 'clusterCol', choices = categorical_cols)
      shiny::updateSelectInput(session, 'sampleMetaCol', choices = colnames(bulk.metadata))
    })

    clusterPropPlot <- shiny::reactive({
      shiny::req(input$clusterCol, input$sampleMetaCol)
      df <- shared_data$updated.metadata
      if (!(input$clusterMetaIncludeNAs)) {
        df <- df[!is.na(df[,input$clusterCol]) & !is.na(df[,input$sampleMetaCol]), ]
      }
      df$sample_meta <- as.factor(df[[input$sampleMetaCol]])
      df$cluster_col <- as.factor(df[[input$clusterCol]])
      df_summary <- df |>
        dplyr::group_by(Sample, .data$sample_meta, .data$cluster_col) |>
        dplyr::summarise(count = dplyr::n()) |>
        dplyr::group_by(Sample, .data$sample_meta) |>
        dplyr::mutate(proportion = .data$count / sum(.data$count))

      # Group by sample and metadata
      ggplot2::ggplot(df_summary, ggplot2::aes(x = .data$sample_meta, y = .data$proportion, fill = .data$cluster_col)) +
        ggplot2::geom_boxplot() +
        ggplot2::theme_minimal()

    })
    output$clusterMetaProportionPlot <- shiny::renderPlot({
      clusterPropPlot()
    })

      output$overlapPlot <- shiny::renderPlot({
        df <- dplyr::select(shared_data$updated.metadata, dplyr::all_of(unique(c(input$col1, input$col2))))
        comp1 <- df[[input$col1]]
        comp2 <- df[[input$col2]]

        na.value <- ifelse(input$includeNAs, 'ifany', 'no')
        tab <- as.data.frame(table(comp1, comp2, useNA = na.value))

        ggplot2::ggplot(tab, ggplot2::aes(x = .data$comp1, y = .data$comp2, fill = .data$Freq)) +
          ggplot2::geom_tile(color = "white") +
          ggplot2::geom_text(ggplot2::aes(label = .data$Freq)) +
          ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
          ggplot2::labs(x = input$col1, y = input$col2, fill = "Count") +
          ggplot2::theme_minimal()
      })
      output[['downloadOverlap']] <- utils_create_download_plot_handler(
        plot_func = function(){
          df <- dplyr::select(shared_data$updated.metadata, dplyr::all_of(unique(c(input$col1, input$col2))))
          comp1 <- df[[input$col1]]
          comp2 <- df[[input$col2]]
          na.value <- ifelse(input$includeNAs, 'ifany', 'no')
          tab <- as.data.frame(table(comp1, comp2, useNA = na.value))
          ggplot2::ggplot(tab, ggplot2::aes(x = .data$comp1, y = .data$comp2, fill = .data$Freq)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::geom_text(ggplot2::aes(label = .data$Freq)) +
            ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
            ggplot2::labs(x = input$col1, y = input$col2, fill = "Count") +
            ggplot2::theme_minimal()
        },
        filename_func = function() input[['overlapFileName']],
        width_func = function() input[['overlapWidth']],
        height_func = function() input[['overlapHeight']]
      )

      output$overlapJSIPlot <- shiny::renderPlot({
        df <- dplyr::select(shared_data$updated.metadata, dplyr::all_of(unique(c(input$col1, input$col2))))
        comp1 <- df[[input$col1]]
        comp2 <- df[[input$col2]]

        na.value = ifelse(input$includeNAs,'ifany','no')
        tab <- as.data.frame(table(comp1, comp2, useNA = na.value))

        jaccard_df <- dplyr::rowwise(tab) |>
          dplyr::mutate(
            intersection = .data$Freq,
            union = sum(df[[input$col1]] == comp1 | df[[input$col2]] == comp2, na.rm = TRUE),
            JSI = ifelse(.data$union > 0, .data$intersection / .data$union, NA)
          )

        ggplot2::ggplot(jaccard_df, ggplot2::aes(x = .data$comp1, y = .data$comp2, fill = .data$JSI)) +
          ggplot2::geom_tile(color = "white") +
          ggplot2::geom_text(ggplot2::aes(label = round(.data$JSI, 2))) +
          ggplot2::scale_fill_gradient(low = "white", high = "darkgreen", na.value = "grey90") +
          ggplot2::labs(x = input$col1, y = input$col2, fill = "Jaccard Index") +
          ggplot2::theme_minimal()

      })
      output[['downloadOverlapJSI']] <- utils_create_download_plot_handler(
        plot_func = function(){
          df <- dplyr::select(shared_data$updated.metadata, dplyr::all_of(unique(c(input$col1, input$col2))))
          comp1 <- df[[input$col1]]
          comp2 <- df[[input$col2]]
          na.value = ifelse(input$includeNAs,'ifany','no')
          tab <- as.data.frame(table(comp1, comp2, useNA = na.value))
          jaccard_df <- dplyr::rowwise(tab) |>
            dplyr::mutate(
              intersection = .data$Freq,
              union = sum(df[[input$col1]] == .data$comp1 | df[[input$col2]] == .data$comp2, na.rm = TRUE),
              JSI = ifelse(.data$union > 0, .data$intersection / .data$union, NA)
            )
          ggplot2::ggplot(jaccard_df, ggplot2::aes(x = .data$comp1, y = .data$comp2, fill = .data$JSI)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::geom_text(ggplot2::aes(label = round(.data$JSI, 2))) +
            ggplot2::scale_fill_gradient(low = "white", high = "darkgreen", na.value = "grey90") +
            ggplot2::labs(x = input$col1, y = input$col2, fill = "Jaccard Index") +
            ggplot2::theme_minimal()
        },
        filename_func = function() input[['overlapJSIFileName']],
        width_func = function() input[['overlapJSIWidth']],
        height_func = function() input[['overlapJSIHeight']]
      )

      propBar <- shiny::reactive({
        if (input$includeNAs == F){
          df = shared_data$updated.metadata[!is.na(shared_data$updated.metadata[,input$col1]) &
                                              !is.na(shared_data$updated.metadata[,input$col2]), ]
        } else {
          df = shared_data$updated.metadata
        }
        df[,input$col2] = factor(df[,input$col2])
        if (input$scaleByNumber) {
          df_summary <- df |>
            dplyr::group_by(.data[[input$col1]], .data[[input$col2]]) |>
            dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
            dplyr::group_by(.data[[input$col2]]) |>
            dplyr::mutate(proportion = .data$count / sum(.data$count))

        } else {
          df_summary <- df |>
            dplyr::group_by(.data[[input$col1]], .data[[input$col2]]) |>
            dplyr::summarise(count = dplyr::n()) |>
            dplyr::ungroup() |>
            dplyr::mutate(proportion = .data$count)
        }
        ggplot2::ggplot(df_summary, ggplot2::aes(x = .data[[input$col1]], fill = .data[[input$col2]], y = .data$proportion)) +
          ggplot2::geom_bar(position='fill', stat = 'identity') +
          ggplot2::theme_minimal()
      })

      output$proportionPlot <- shiny::renderPlot({
        propBar()
      })

      output[['downloadProportion']] <- utils_create_download_plot_handler(
        plot_func = propBar,
        filename_func = function() input[['proportionFileName']],
        width_func = function() input[['proportionWidth']],
        height_func = function() input[['proportionHeight']]
      )

      output$ecsSpatial <- shiny::renderPlot({
        df = shared_data$updated.metadata[!is.na(shared_data$updated.metadata[,input$col1]) &
                                            !is.na(shared_data$updated.metadata[,input$col2]), ]
        df$sim = ClustAssess::element_sim_elscore(df[,input$col1],df[,input$col2])

        ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data$sim, fill = .data$sim)) +
          ggplot2::geom_tile() +
          ggplot2::theme_classic() +
          ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                         axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank(),
                         axis.line.x = ggplot2::element_blank(),
                         axis.title.y=ggplot2::element_blank(),
                         axis.text.y=ggplot2::element_blank(),
                         axis.ticks.y=ggplot2::element_blank(),
                         axis.line.y = ggplot2::element_blank()) + ggplot2::facet_wrap(ggplot2::vars(Sample),scales='free', nrow = max(1,floor(sqrt(length(unique(df$Sample))/1.5)))) +
          viridis::scale_color_viridis() + viridis::scale_fill_viridis() + ggplot2::theme(aspect.ratio = 1)
      })
      output[['downloadECSSpatial']] <- utils_create_download_plot_handler(
        plot_func = function(){
          df = shared_data$updated.metadata[!is.na(shared_data$updated.metadata[,input$col1]) &
                                              !is.na(shared_data$updated.metadata[,input$col2]), ]
          df$sim = ClustAssess::element_sim_elscore(df[,input$col1],df[,input$col2])
          ggplot2::ggplot(df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, color = .data$sim, fill = .data$sim)) +
            ggplot2::geom_tile() +
            ggplot2::theme_classic() +
            ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank(),
                           axis.line.x = ggplot2::element_blank(),
                           axis.title.y=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks.y=ggplot2::element_blank(),
                           axis.line.y = ggplot2::element_blank()) +
                           ggplot2::facet_wrap(ggplot2::vars(Sample), scales='free', nrow = max(1,floor(sqrt(length(unique(df$Sample))/1.5)))) +
            viridis::scale_color_viridis() + viridis::scale_fill_viridis() + ggplot2::theme(aspect.ratio = 1)
        },
        filename_func = function() input[['ecsSpatialFileName']],
        width_func = function() input[['ecsSpatialWidth']],
        height_func = function() input[['ecsSpatialHeight']]
      )

      output[['downloadClusterMetaProportion']] <- utils_create_download_plot_handler(
      plot_func = clusterPropPlot,
      filename_func = function() input[['clusterMetaProportionFileName']],
      width_func = function() input[['clusterMetaProportionWidth']],
      height_func = function() input[['clusterMetaProportionHeight']]
    )
  })
}



