
#' Performs and visualises multi-group comparisons
#'
#' @description UI and server logic for the Bulk Multi-Comparison panel, enabling statistical comparison of multiple groups in the bulk dataset using ANOVA or Kruskal-Wallis tests, with post-hoc analysis, boxplots, cross-plots, heatmaps, and Venn diagrams.
#'
#' @details
#' \itemize{
#'   \item{Allows users to perform multi-group comparisons using ANOVA or Kruskal-Wallis tests on bulk-level data.}
#'   \item{Displays results in an interactive table, with options to view and download ANOVA/Kruskal results, post-hoc boxplots, and significance tiles.}
#'   \item{Cross-plot: Compare two independent group contrasts, visualise log2 fold changes, and overlap of significant peaks via Venn diagram.}
#'   \item{Heatmaps: Visualise scaled intensities and log2 fold changes for significant peaks across groups.}
#'   \item{All tables and plots are downloadable via the download menus, with customisable file names and plot sizes.}
#'   \item{Tip: Click a row in the ANOVA table to view post-hoc plots for that peak.}
#'   \item{Only groups with at least 2 samples are available for selection in comparisons.}
#' }
#'
#' @name BulkPanel_MultiCompareTab
#' @rdname BulkPanel_MultiCompareTab

utils::globalVariables(c("ComparisonGroup", "Peak", "Sample", "LFCGroup", "HeatmapGroup"))

#' @rdname BulkPanel_MultiCompareTab
#' @param id Shiny module id (for both UI and server)
#' @param bulk.metadata Data frame of bulk sample metadata
#' @param bulk.intensity.matrix Matrix of bulk sample intensities
#' @param show Logical; whether to render the panel (default: TRUE)
#' @export
BulkPanel_MultiCompareTabUI <- function(id, bulk.metadata, bulk.intensity.matrix, show = TRUE) {
  ns <- shiny::NS(id)
  if(show){
  shiny::tabPanel(
    "Multiple Comparisons",
    shiny::fluidPage(
      shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
        shinyWidgets::dropMenu(
          shinyWidgets::circleButton(ns("info_multi_compare"), icon = shiny::icon("info-circle"), status = "info"),
          shiny::tags$div(
            shiny::tags$h3("Bulk Multi-Comparison Panel Info"),
            shiny::tags$ul(
              shiny::tags$li("Perform multi-group comparisons using ANOVA or Kruskal-Wallis tests."),
              shiny::tags$li("View and download ANOVA/Kruskal results, post-hoc boxplots, and significance tiles."),
              shiny::tags$li("Cross-plot: Compare two independent group contrasts, visualise log2 fold changes, and overlap of significant peaks via Venn diagram."),
              shiny::tags$li("All tables and plots are downloadable via the download menus."),
              shiny::tags$li("Settings for plot downloads are available in the download menu."),
              shiny::tags$li("Tip: Click a row in the ANOVA table to view post-hoc plots for that peak."),
              shiny::tags$li(shiny::tags$b("Note: Only groups with at least 2 samples are available for selection in comparisons."))
            )
          ),
          theme = "light-border",
          placement = "right",
          arrow = FALSE
        )
      ),
      shiny::h3("ANOVA and Kruskal-Wallis Comparisons"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::selectInput(ns("group_col"), "Metadata column for grouping:", choices = names(bulk.metadata), selected = names(bulk.metadata)[1]),
          shiny::uiOutput(ns("group_choices_ui")),
          shiny::selectInput(ns("test"), "Statistical test:", choices = c("ANOVA", "Kruskal-Wallis")),
          shiny::actionButton(ns("run"), "Run comparison"),
          shiny::hr(),
          shiny::textOutput(ns("posthoc"))
        ),
        shiny::mainPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_anova"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download ANOVA/Kruskal Results"),
                shiny::textInput(ns('anova_table_filename'), 'ANOVA Table file name', value = 'anova_table.csv'),
                shiny::downloadButton(ns("download_anova_table"), "Download ANOVA Table (CSV)"),
                shiny::textInput(ns('selected_row_filename'), 'Selected Row file name', value = 'anova_selected_row.csv'),
                shiny::downloadButton(ns("download_selected_row"), "Download Selected Row (CSV)"),
                shiny::textInput(ns('posthoc_plot_filename'), 'Boxplot/Tile Plot file name', value = 'posthoc_plot.png'),
                shiny::numericInput(ns('posthoc_plot_width'), 'Plot width (in)', value = 10, min = 1, max = 50, step = 1),
                shiny::numericInput(ns('posthoc_plot_height'), 'Plot height (in)', value = 5, min = 1, max = 50, step = 1),
                shiny::downloadButton(ns("download_posthoc_plot"), "Download box/tile plot"),
                shiny::tags$hr(),
                shiny::tags$h3("Download Heatmaps"),
                shiny::textInput(ns('signif_heatmap_filename'), 'Scaled intensity heatmap file name', value = 'signif_heatmap.png'),
                shiny::numericInput(ns('signif_heatmap_width'), 'Scaled intensity heatmap width (in)', value = 10, min = 1, max = 50, step = 1),
                shiny::numericInput(ns('signif_heatmap_height'), 'Scaled intensity heatmap height (in)', value = 6, min = 1, max = 50, step = 1),
                shiny::downloadButton(ns("download_signif_heatmap"), "Download scaled intensity heatmap"),
                shiny::textInput(ns('lfc_heatmap_filename'), 'LFC heatmap file name', value = 'lfc_heatmap.png'),
                shiny::numericInput(ns('lfc_heatmap_width'), 'LFC heatmap width (in)', value = 10, min = 1, max = 50, step = 1),
                shiny::numericInput(ns('lfc_heatmap_height'), 'LFC heatmap height (in)', value = 6, min = 1, max = 50, step = 1),
                shiny::downloadButton(ns("download_lfc_heatmap"), "Download LFC heatmap")
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          DT::DTOutput(ns("anovaTable")),
          shiny::plotOutput(ns("posthocPlot"))
        )),
          shiny::sidebarLayout(
            shiny::sidebarPanel(
              shiny::numericInput(ns("heatmap_npeaks"), "Number of top peaks to show in heatmap:", value = 30, min = 1, max = 200, step = 1),
              shiny::selectInput(ns("heatmap_group_col"), "Metadata column for heatmap aggregation:", choices = names(bulk.metadata), selected = names(bulk.metadata)[1]),
              shiny::actionButton(ns("run_heatmap"), "Run heatmap")
            ),
            shiny::mainPanel(
              plotly::plotlyOutput(ns("signifHeatmap"),height='600px')
            )
          ),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
            shiny::selectizeInput(ns('user_peaks'), 'Peaks to include in heatmaps:', choices = NULL, multiple = TRUE, options = list(placeholder = 'Select peaks or leave blank for top N.')),
            shiny::helpText("If no peaks are selected, the top peaks (by adjusted p-value) will be used for both heatmaps."),

          shiny::fluidRow(
      shiny::column(6,
        shiny::selectInput(ns("lfc1_group_col"), "Grouping column (LFC 1):", choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[1]),
        shiny::selectInput(ns("lfc1_treated"), "Treated group(s) (LFC 1):", choices = c(), multiple = TRUE),
        shiny::selectInput(ns("lfc1_control"), "Control group(s) (LFC 1):", choices = c(), multiple = TRUE)
      ),
      shiny::column(6,
        shiny::selectInput(ns("lfc2_group_col"), "Grouping column (LFC 2):", choices = colnames(bulk.metadata), selected = colnames(bulk.metadata)[1]),
        shiny::selectInput(ns("lfc2_treated"), "Treated group(s) (LFC 2):", choices = c(), multiple = TRUE),
        shiny::selectInput(ns("lfc2_control"), "Control group(s) (LFC 2):", choices = c(), multiple = TRUE)
                ),
      shiny::actionButton(ns("run_lfc_heatmaps"), "Run LFC heatmaps")
              )
            ),
    shiny::mainPanel(
      plotly::plotlyOutput(ns("lfcHeatmap1"), height = '600px')
          )
          ),
      shiny::h3("Cross plot"),
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::tags$h4("Cross-plot Controls"),
          shiny::selectInput(ns("group_col"), "Metadata column for grouping:", choices = names(bulk.metadata), selected = names(bulk.metadata)[1]),
          shiny::uiOutput(ns("crossplot_pair_select")),
          shiny::uiOutput(ns("crossplot_pair_select2")),
          shiny::selectInput(ns("crossplot_test"), "Cross plot test:", choices = c("t-test", "Wilcox rank sum"), selected = "t-test"),
          shiny::sliderInput(ns("crossplot_pval"), "Adjusted p-value threshold", min = 0, max = 1, value = 0.05, step = 0.005),
          shiny::sliderInput(ns("crossplot_lfc"), "log2 fold change threshold", min = 0, max = 5, value = 1, step = 0.1),
          shiny::actionButton(ns("run_crossplot"), "Run cross plot")
        ),
        shiny::mainPanel(
          shiny::div(style = "display: flex; gap: 10px; align-items: center; margin-bottom: 10px;",
            shinyWidgets::dropMenu(
              shinyWidgets::circleButton(ns("downloads_crossplot"), icon = shiny::icon("download"), status = "success"),
              shiny::tags$div(
                shiny::tags$h3("Download cross plot results"),
                shiny::textInput(ns('crossplot_filename'), 'Cross plot file name', value = 'crossplot.png'),
                shiny::numericInput(ns('crossplot_width'), 'Plot width (in)', value = 7, min = 1, max = 50, step = 1),
                shiny::numericInput(ns('crossplot_height'), 'Plot height (in)', value = 7, min = 1, max = 50, step = 1),
                shiny::downloadButton(ns("download_crossplot"), "Download cross plot"),
                shiny::textInput(ns('crossvenn_filename'), 'Cross plot venn diagram file name', value = 'crossvenn.png'),
                shiny::numericInput(ns('crossvenn_width'), 'Venn diagram width (in)', value = 6, min = 1, max = 50, step = 1),
                shiny::numericInput(ns('crossvenn_height'), 'Venn diagram height (in)', value = 5, min = 1, max = 50, step = 1),
                shiny::downloadButton(ns("download_crossvenn"), "Download cross plot venn diagram"),
                shiny::textInput(ns('crossplot_table_filename'), 'Cross plot table file name', value = 'crossplot_table.csv'),
                shiny::downloadButton(ns("download_crossplot_table"), "Download cross plot table")
              ),
              theme = "light-border",
              placement = "right",
              arrow = FALSE
            )
          ),
          plotly::plotlyOutput(ns("crossPlot"), height = "600px"),
          shiny::plotOutput(ns("crossVenn"), height = "300px")
        )
      )
    )
  )
  } else {
    NULL
  }
}

#' @param anno Data frame of peak annotations
#' @rdname BulkPanel_MultiCompareTab
#' @export
BulkPanel_MultiCompareTabServer <- function(id, bulk.intensity.matrix, bulk.metadata, anno) {

  anova_table_last <- shiny::reactiveVal(NULL)
  df_last <- shiny::reactiveVal(NULL)
  crossplot_last <- shiny::reactiveVal(NULL)
  crossvenn_last <- shiny::reactiveVal(NULL)

  # Function to run post-hoc test for selected peak
  run_posthoc_test <- function(df, peak, test_type) {
    subdf <- df[df$Peak == peak, ]
    if (length(unique(subdf$ComparisonGroup)) < 2) return(NULL)
    if (test_type == "ANOVA") {
      fit <- stats::aov(Intensity ~ ComparisonGroup, data = subdf)
      tukey <- stats::TukeyHSD(fit)
      res <- as.data.frame(tukey$ComparisonGroup)
      res$Comparison <- rownames(tukey$ComparisonGroup)
      res <- res[, c("Comparison", "diff", "lwr", "upr", "p adj")]
      colnames(res) <- c("Comparison", "Difference", "Lower", "Upper", "p_adj")
      return(res)
    } else if (test_type == "Kruskal-Wallis") {
      dunn <- dunn.test::dunn.test(subdf$Intensity, subdf$ComparisonGroup, method = "BH")
      res <- data.frame("Comparison" = dunn$comparisons, "Z" = dunn$Z, "p_adj" = dunn$P.adjusted)
      return(res)
    }
    NULL
  }
  selected_peak <- shiny::reactiveVal(NULL)

  shiny::moduleServer(id, function(input, output, session) {

  output$posthoc <- shiny::renderText({
      posthoc_test = if (selected_test() == "ANOVA") "Tukey HSD" else "Dunn"
      paste0("Select a peak from the table to perform ", posthoc_test, " post-hoc test.")
  })

    # UI for group selection based on chosen metadata column
    output$group_choices_ui <- shiny::renderUI({
      shiny::req(input$group_col)
      # Only offer groups with at least 2 samples
      group_counts <- table(bulk.metadata[[input$group_col]])
      valid_choices <- names(group_counts)[group_counts >= 2]
      shiny::selectInput(session$ns("groups"), "Groups to compare:", choices = valid_choices, multiple = TRUE, selected = valid_choices)
    })

    # UI for first comparison (A vs B)
    output$crossplot_pair_select <- shiny::renderUI({
      shiny::req(input$group_col)
      # Only offer groups with at least 2 samples
      group_counts <- table(bulk.metadata[[input$group_col]])
      valid_choices <- names(group_counts)[group_counts >= 2]
      shiny::tagList(
        shiny::selectInput(session$ns("crossplot_group1a"), "Comparison 1: Group A", choices = valid_choices, selected = valid_choices[1]),
        shiny::selectInput(session$ns("crossplot_group1b"), "Comparison 1: Group B", choices = valid_choices, selected = valid_choices[min(2, length(valid_choices))])
      )
    })

    # UI for second comparison (C vs D)
    output$crossplot_pair_select2 <- shiny::renderUI({
      shiny::req(input$group_col)
      # Only offer groups with at least 2 samples
      group_counts <- table(bulk.metadata[[input$group_col]])
      valid_choices <- names(group_counts)[group_counts >= 2]
      shiny::tagList(
        shiny::selectInput(session$ns("crossplot_group2a"), "Comparison 2: Group C", choices = valid_choices, selected = valid_choices[min(3, length(valid_choices))]),
        shiny::selectInput(session$ns("crossplot_group2b"), "Comparison 2: Group D", choices = valid_choices, selected = valid_choices[min(4, length(valid_choices))])
      )
    })

        # Update selectizeInput for user_peaks with all available peaks
    shiny::observe({
      shiny::updateSelectizeInput(
        session,
        'user_peaks',
        choices = rownames(bulk.intensity.matrix),
        server = TRUE
      )
    })

    selected_groups <- shiny::reactive({ input$groups })
    selected_test <- shiny::reactive({ input$test }) |> shiny::bindEvent(input$run)

    shiny::observeEvent(input$run, {
      shiny::req(selected_groups())
      shiny::req(length(selected_groups()) > 1)
      idx <- bulk.metadata[[input$group_col]] %in% selected_groups()
      data_sub <- as.data.frame(bulk.intensity.matrix[, idx, drop = FALSE])
      meta_sub <- bulk.metadata[idx, , drop = FALSE]
      if (is.null(rownames(data_sub))) rownames(data_sub) <- seq_len(nrow(data_sub))
      if (is.null(colnames(data_sub))) colnames(data_sub) <- seq_len(ncol(data_sub))
      data_sub$Peak <- rownames(data_sub)
      df <- tidyr::pivot_longer(data_sub,cols=1:(ncol(data_sub)-1), names_to = "Sample", values_to = "Intensity")
      if (ncol(df) == 3) colnames(df) <- c("Peak", "Sample", "Intensity")
      df$ComparisonGroup <- meta_sub[[input$group_col]][match(df$Sample, meta_sub$Sample)]

      # ANOVA or Kruskal-Wallis
      shiny::withProgress(message = 'Running multi-group comparison...', value = 0, {
      shiny::incProgress(0.3)
      if (selected_test() == "ANOVA") {
        anova_table <- data.frame("Peak" = unique(df$Peak))
        anova_table$pval <- sapply(anova_table$Peak, function(pk) {
          subdf <- df[df$Peak == pk, ]
          fit <- stats::aov(Intensity ~ ComparisonGroup, data = subdf)
          return(summary(fit)[[1]]["Pr(>F)"][[1]][1])
         })
        anova_table$adj_pval <- stats::p.adjust(anova_table$pval, method = "BH")
      } else {
        anova_table <- data.frame("Peak" = unique(df$Peak))
        anova_table$pval <- sapply(anova_table$Peak, function(pk) {
          subdf <- df[df$Peak == pk, ]
          fit <- stats::kruskal.test(Intensity ~ ComparisonGroup, data = subdf)
          return(fit$p.value)
        })
        anova_table$adj_pval <- stats::p.adjust(anova_table$pval, method = "BH")
      }
      anova_table <- merge(anova_table, anno[, c("m_z", "display_name")], by.x = "Peak", by.y = "m_z", all.x = TRUE)
      anova_table <- anova_table[base::order(anova_table$adj_pval), ]
      anova_table_last(anova_table)
      df_last(df)
      shiny::incProgress(0.9)
      output$anovaTable <- DT::renderDT({
        DT::datatable(anova_table, selection = 'single', options = list(pageLength = 10))
      })
      })
      shiny::observeEvent(input$anovaTable_rows_selected, {
        sel <- input$anovaTable_rows_selected
        if (!is.null(sel) && length(sel) == 1) {
          selected_peak(anova_table$Peak[sel])
        }
      })
    })


    # Download ANOVA table
    output$download_anova_table <- shiny::downloadHandler(
      filename = function() input$anova_table_filename,
      content = function(file) {
        atab <- anova_table_last()
        if (!is.null(atab)) {
          utils::write.csv(atab, file, row.names = FALSE)
        } else {
          shiny::showNotification("Please run the comparison first.", type = "error")
        }
      }
    )
    output$download_selected_row <- shiny::downloadHandler(
      filename = function() input$selected_row_filename,
      content = function(file) {
        atab <- anova_table_last()
        sel <- input$anovaTable_rows_selected
        if (!is.null(atab) && !is.null(sel) && length(sel) == 1) {
          utils::write.csv(atab[sel, , drop = FALSE], file, row.names = FALSE)
        } else {
          shiny::showNotification("Please select a row in the ANOVA table first.", type = "error")
        }
      }
    )

    posthoc_plot <- shiny::reactive({
        pk <- selected_peak()
        if (is.null(pk)) return(NULL)
        df = df_last()
        subdf <- df[df$Peak == pk, ]
        plot1 <- ggplot2::ggplot(subdf, ggplot2::aes(x = .data$ComparisonGroup, y = .data$Intensity, fill = .data$ComparisonGroup)) +
          ggplot2::geom_boxplot() +
          ggplot2::ggtitle(paste("Peak:", pk)) +
          ggplot2::theme_classic()
        res <- run_posthoc_test(df, pk, selected_test())
        if (is.null(res)) return(NULL)
        res$star <- cut(res[,'p_adj'], breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c('***', '**', '*', ''))
        res$pval_round <- signif(res[,'p_adj'], 2)
        pairs <- strsplit(as.character(gsub(' ','',res$Comparison)), "-")
        groups <- sort(unique(unlist(pairs)))
        tile_df <- data.frame(
          "Var1" = vapply(pairs, `[`, character(1), 1),
          "Var2" = vapply(pairs, `[`, character(1), 2),
          "pval" = res$pval_round,
          "star" = res$star
        )
        # Fill in symmetric tiles
        tile_df_sym <- rbind(tile_df, stats::setNames(tile_df, c('Var2', 'Var1', 'pval', 'star')))
        plot2 <- ggplot2::ggplot(tile_df_sym, ggplot2::aes(x = .data$Var1, y = .data$Var2, fill = .data$pval)) +
          ggplot2::geom_tile(color = "grey80") +
          ggplot2::geom_text(ggplot2::aes(label = paste0(.data$pval, '\n', .data$star)), size = 5) +
          ggplot2::scale_fill_gradient(high = "white", low = "red") +
          ggplot2::scale_x_discrete(limits = groups) +
          ggplot2::scale_y_discrete(limits = groups) +
          ggplot2::theme_classic() +
          ggplot2::coord_fixed()
        patchwork::wrap_plots(plot1, plot2, ncol = 2) + patchwork::plot_annotation(title = paste("Post-hoc test for peak", pk))
    })

    # Boxplot for selected peak
    output$posthocPlot <- shiny::renderPlot({
      posthoc_plot()
    })

    output[['download_posthoc_plot']] <- utils_create_download_plot_handler(
      plot_func = posthoc_plot,
      filename_func = function() input[['posthoc_plot_filename']],
      width_func = function() input[['posthoc_plot_width']],
      height_func = function() input[['posthoc_plot_height']]
    )

    # Cross-plot analysis (pairwise DA)
    shiny::observeEvent(input$run_crossplot, {

      shiny::withProgress(message = 'Computing cross plot...', value = 0, {
      shiny::req(input$crossplot_group1a, input$crossplot_group1b, input$crossplot_group2a, input$crossplot_group2b)
      test_type <- input$crossplot_test
      # First comparison (A vs B)
      idx1 <- bulk.metadata[[input$group_col]] %in% c(input$crossplot_group1a, input$crossplot_group1b)
      data_sub1 <- as.data.frame(bulk.intensity.matrix[, idx1, drop = FALSE])
      meta_sub1 <- bulk.metadata[idx1, , drop = FALSE]
      if (is.null(rownames(data_sub1))) rownames(data_sub1) <- seq_len(nrow(data_sub1))
      if (is.null(colnames(data_sub1))) colnames(data_sub1) <- seq_len(ncol(data_sub1))
      group_col <- input$group_col
      shiny::incProgress(0.3)
      da1 <- bulk_utils_DA(data_sub1, meta_sub1[[group_col]], input$crossplot_group1a, input$crossplot_group1b, test = test_type, anno = anno)

      # Second comparison (C vs D)
      idx2 <- bulk.metadata[[input$group_col]] %in% c(input$crossplot_group2a, input$crossplot_group2b)
      data_sub2 <- as.data.frame(bulk.intensity.matrix[, idx2, drop = FALSE])
      meta_sub2 <- bulk.metadata[idx2, , drop = FALSE]
      if (is.null(rownames(data_sub2))) rownames(data_sub2) <- seq_len(nrow(data_sub2))
      if (is.null(colnames(data_sub2))) colnames(data_sub2) <- seq_len(ncol(data_sub2))
      shiny::incProgress(0.6)
      da2 <- bulk_utils_DA(data_sub2, meta_sub2[[group_col]], input$crossplot_group2a, input$crossplot_group2b, test = test_type, anno = anno)

      # Merge by m_z and add annotation
      merged <- merge(da1[, c('m_z', 'pvalAdj', 'lfc')], da2[, c('m_z', 'pvalAdj', 'lfc')], by = 'm_z', suffixes = c('_1', '_2'))
      merged <- merge(merged, anno, by = 'm_z', all.x = TRUE)
      # Use user-set thresholds
      pval_thresh <- input$crossplot_pval
      lfc_thresh <- input$crossplot_lfc
      merged$DE1 <- abs(merged$lfc_1) > lfc_thresh & merged$pvalAdj_1 < pval_thresh
      merged$DE2 <- abs(merged$lfc_2) > lfc_thresh & merged$pvalAdj_2 < pval_thresh
      merged$DEgroup <- 'not DE'
      merged$DEgroup[merged$DE1 & merged$DE2] <- 'DE both'
      merged$DEgroup[merged$DE1 & !merged$DE2] <- paste0('DE: ', input$crossplot_group1a, ' vs ', input$crossplot_group1b)
      merged$DEgroup[!merged$DE1 & merged$DE2] <- paste0('DE: ', input$crossplot_group2a, ' vs ', input$crossplot_group2b)
      merged$DEgroup <- factor(merged$DEgroup, levels = c('not DE', 'DE both',
        paste0('DE: ', input$crossplot_group1a, ' vs ', input$crossplot_group1b),
        paste0('DE: ', input$crossplot_group2a, ' vs ', input$crossplot_group2b)))
      # Interactive cross plot with annotation
      shiny::incProgress(0.9)

      output$crossPlot <- plotly::renderPlotly({
        plotly::plot_ly(
          data = merged,
          x = ~.data$lfc_1,
          y = ~.data$lfc_2,
          color = ~.data$DEgroup,
          colors = c('grey', 'purple', 'dodgerblue', 'lightcoral'),
          type = 'scatter',
          mode = 'markers',
          marker = list(size = 8, opacity = 0.7),
          text = ~paste0(
            'm/z: ', .data$m_z,
            if ('display_name' %in% colnames(merged)) paste0('<br>Name: ', .data$display_name) else '',
            '<br>LFC 1: ', signif(.data$lfc_1, 3),
            '<br>LFC 2: ', signif(.data$lfc_2, 3),
            '<br>DE group: ', .data$DEgroup
          ),
          hoverinfo = 'text'
        ) |>
          plotly::layout(
            xaxis = list(title = paste('LFC:', input$crossplot_group1a, 'vs', input$crossplot_group1b)),
            yaxis = list(title = paste('LFC:', input$crossplot_group2a, 'vs', input$crossplot_group2b), scaleanchor = 'x', scaleratio = 1),
            title = 'Cross plot of log2 fold changes',
            legend = list(title = list(text = 'DE group'))
          )
      })

      # Venn diagram for overlap in significant peaks
      output$crossVenn <- shiny::renderPlot({
        sig1 <- merged$m_z[merged$DE1]
        sig2 <- merged$m_z[merged$DE2]
        venn_list <- stats::setNames(list(sig1, sig2),
          c(paste0('DE: ', input$crossplot_group1a, ' vs ', input$crossplot_group1b),
            paste0('DE: ', input$crossplot_group2a, ' vs ', input$crossplot_group2b)))
        ggVennDiagram::ggVennDiagram(venn_list, label_alpha = 0.7) +
          ggplot2::theme_void()
      })

        crossplot_last(list(merged = merged, input = list(
            group1a = input$crossplot_group1a,
            group1b = input$crossplot_group1b,
            group2a = input$crossplot_group2a,
            group2b = input$crossplot_group2b
        )))
        crossvenn_last(list(merged = merged, input = list(
            group1a = input$crossplot_group1a,
            group1b = input$crossplot_group1b,
            group2a = input$crossplot_group2a,
            group2b = input$crossplot_group2b
        )))

        # Download cross plot
        output$download_crossplot <- shiny::downloadHandler(
        filename = function() input$crossplot_filename,
        content = function(file) {
          cp <- crossplot_last()
          if (is.null(cp)) { shiny::showNotification("Please run the cross plot first.", type = "error"); return() }
            merged <- cp$merged
            inputvals <- cp$input
            p <- ggplot2::ggplot(merged, ggplot2::aes(x = .data$lfc_1, y = .data$lfc_2, color = .data$DEgroup)) +
            ggplot2::geom_point(size = 2, alpha = 0.7) +
            ggplot2::scale_color_manual(values = c('grey', 'purple', 'dodgerblue', 'lightcoral')) +
            ggplot2::labs(
                x = paste('LFC:', inputvals$group1a, 'vs', inputvals$group1b),
                y = paste('LFC:', inputvals$group2a, 'vs', inputvals$group2b),
                title = 'Cross-plot of log2 fold changes',
                color = 'DE group'
            ) +
            ggplot2::theme_classic() +
            ggplot2::coord_fixed()
            ggplot2::ggsave(file, plot = p, width = input$crossplot_width, height = input$crossplot_height)
        }
        )
        output$download_crossvenn <- shiny::downloadHandler(
        filename = function() input$crossvenn_filename,
        content = function(file) {
          cv <- crossvenn_last()
          if (is.null(cv)) { shiny::showNotification("Please run the cross plot first.", type = "error"); return() }
            merged <- cv$merged
            inputvals <- cv$input
            sig1 <- merged$m_z[merged$DE1]
            sig2 <- merged$m_z[merged$DE2]
          venn_list <- stats::setNames(list(sig1, sig2),
            c(paste0('DE: ', inputvals$group1a, ' vs ', inputvals$group1b),
                paste0('DE: ', inputvals$group2a, ' vs ', inputvals$group2b)))
            p <- ggVennDiagram::ggVennDiagram(venn_list, label_alpha = 0.7) + ggplot2::theme_void()
            ggplot2::ggsave(file, plot = p, width = input$crossvenn_width, height = input$crossvenn_height)
        }
        )
        output$download_crossplot_table <- shiny::downloadHandler(
        filename = function() input$crossplot_table_filename,
        content = function(file) {
          cp <- crossplot_last()
          if (is.null(cp)) { shiny::showNotification("Please run the cross plot first.", type = "error"); return() }
            merged <- cp$merged
          utils::write.csv(merged, file, row.names = FALSE)
        }
        )
      })
    })

    shiny::observeEvent(input$lfc1_group_col, {
      groups <- unique(bulk.metadata[[input$lfc1_group_col]])
      groups <- groups[!is.na(groups)]
      shiny::updateSelectInput(session, "lfc1_treated", choices = groups)
      shiny::updateSelectInput(session, "lfc1_control", choices = groups)
    })

    shiny::observeEvent(input$lfc2_group_col, {
      groups <- unique(bulk.metadata[[input$lfc2_group_col]])
      groups <- groups[!is.na(groups)]
      shiny::updateSelectInput(session, "lfc2_treated", choices = groups)
      shiny::updateSelectInput(session, "lfc2_control", choices = groups)
    })

    signif_heatmap_preproc <- shiny::reactive({
      # Only run when button pressed
      # 1. Get the latest ANOVA/Kruskal results and long-format data
      atab <- anova_table_last()
      dff <- df_last()
      shiny::req(atab, dff, input$heatmap_npeaks, input$heatmap_group_col)

      # 2. Filter for significant peaks (adj_pval < 0.05)
      sig_atab <- atab[atab$adj_pval < 0.05, ]
      if (nrow(sig_atab) == 0) {
        # No significant peaks found
        stop("No significant peaks found.")
      }

      # 3. Use user-selected peaks or top N by adjusted p-value
      sig_atab <- sig_atab[order(sig_atab$adj_pval), ]
      user_peaks <- input$user_peaks
      if (!is.null(user_peaks) && length(user_peaks) > 0) {
        selected_peaks <- intersect(sig_atab$Peak, user_peaks)
      } else {
        selected_peaks <- utils::head(sig_atab$Peak, input$heatmap_npeaks)
      }
      dff_sig <- dff[dff$Peak %in% selected_peaks, ]
      if (nrow(dff_sig) == 0) {
        # No data for selected peaks
        stop("No data for selected peaks.")
      }

      # 4. Assign group for aggregation using selected metadata column
      group_col <- input$heatmap_group_col
      meta_samples <- bulk.metadata$Sample
      dff_sig$HeatmapGroup <- NA
      if (!is.null(meta_samples)) {
        match_idx <- match(as.character(dff_sig$Sample), as.character(meta_samples))
        dff_sig$HeatmapGroup <- bulk.metadata[[group_col]][match_idx]
      }
      if (all(is.na(dff_sig$HeatmapGroup))) {
        # No matching group assignments
        stop("No matching group assignments for selected metadata column.")
      }

      # 5. Aggregate mean intensity by Peak and HeatmapGroup using dplyr
      #    (instead of base aggregate)
      agg <- dff_sig |>
        dplyr::group_by(.data$Peak, .data$HeatmapGroup) |>
        dplyr::summarise(Intensity = mean(.data$Intensity, na.rm = TRUE), .groups = 'drop')
      if (nrow(agg) == 0) {
        # No aggregated data available
        stop("No aggregated data available.")
      }

      # 6. Cast to wide format: rows=peaks, columns=groups
      heatmap_df <- reshape2::dcast(agg, Peak ~ HeatmapGroup, value.var = "Intensity")
      if (nrow(heatmap_df) == 0 || ncol(heatmap_df) <= 1) {
        # Not enough data for heatmap
        stop("Not enough data for heatmap.")
      }
      rownames(heatmap_df) <- heatmap_df$Peak
      mat <- as.matrix(heatmap_df[, -1, drop=FALSE])
      if (nrow(mat) == 0 || ncol(mat) == 0) {
        # Not enough data for heatmap
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::geom_text(ggplot2::aes(0.5,0.5,label="Not enough data for heatmap."), size=6))
      }

      # 7. Z-score normalization by row (peak)
      mat_z <- t(scale(t(mat)))
      return(list(mat_z = mat_z, group_col = group_col))
    }) |> shiny::bindEvent(input$run_heatmap)

    output$signifHeatmap <- plotly::renderPlotly({
      res <- signif_heatmap_preproc()
      shiny::req(res)
        # Use heatmaply for interactive clustering
        heatmaply::heatmaply(
          res$mat_z,
          Rowv = TRUE,
          Colv = TRUE,
          colors = grDevices::colorRampPalette(c(scales::muted('blue'), 'white', scales::muted('red')))(100),
          xlab = res$group_col,
          ylab = "Peak",
          main = paste0("Top ", input$heatmap_npeaks, " Significant Peaks Heatmap (", res$group_col, ")"),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red'), midpoint = 0, name = "Z-score"),
          dendrogram = "both",
          showticklabels = c(TRUE, TRUE),
          fontsize_row = 8,
          fontsize_col = 10,
          hide_colorbar = FALSE
        )
    })

    # Helper to compute LFC matrix for a set of treated/control groups and a grouping column
    compute_lfc_matrix <- function(group_col, treated, control) {
      if (is.null(group_col) || length(treated) == 0 || length(control) == 0) return(NULL)
      # 1. Long-format data
      data_sub <- as.data.frame(bulk.intensity.matrix)
      if (is.null(rownames(data_sub))) rownames(data_sub) <- seq_len(nrow(data_sub))
      if (is.null(colnames(data_sub))) colnames(data_sub) <- seq_len(ncol(data_sub))
      data_sub$Peak <- rownames(data_sub)
      df <- tidyr::pivot_longer(data_sub, cols = 1:(ncol(data_sub)-1), names_to = "Sample", values_to = "Intensity")
      if (ncol(df) == 3) colnames(df) <- c("Peak", "Sample", "Intensity")
      # 2. Assign group for aggregation using selected metadata column
      meta_samples <- bulk.metadata$Sample
      df$LFCGroup <- NA
      if (!is.null(meta_samples)) {
        match_idx <- match(as.character(df$Sample), as.character(meta_samples))
        df$LFCGroup <- bulk.metadata[[group_col]][match_idx]
      }
      # 3. Aggregate mean intensity by Peak and LFCGroup
      agg <- df |>
        dplyr::filter(.data$LFCGroup %in% c(treated, control)) |>
        dplyr::group_by(.data$Peak, .data$LFCGroup) |>
        dplyr::summarise(Intensity = mean(.data$Intensity, na.rm = TRUE), .groups = 'drop')
      # 4. Cast to wide format: rows=peaks, columns=groups
      heatmap_df <- reshape2::dcast(agg, Peak ~ LFCGroup, value.var = "Intensity")
      if (nrow(heatmap_df) == 0) return(NULL)
      # 5. For each treated group, calculate LFC vs mean of all controls
      control_means <- rowMeans(heatmap_df[, control, drop = FALSE], na.rm = TRUE)
      lfc_mat <- sapply(treated, function(tr) {
        tr_means <- heatmap_df[, tr]
        log2(tr_means / control_means)
      })
      if (is.null(dim(lfc_mat))) lfc_mat <- matrix(lfc_mat, ncol = 1)
      rownames(lfc_mat) <- heatmap_df$Peak
      colnames(lfc_mat) <- treated
      lfc_mat[!is.finite(lfc_mat)] <- NA
      lfc_mat <- lfc_mat[apply(lfc_mat, 1, function(x) any(!is.na(x))), , drop=FALSE]
      lfc_mat
    }

    shiny::observeEvent(input$run_lfc_heatmaps, {
            mat1 <- compute_lfc_matrix(input$lfc1_group_col, input$lfc1_treated, input$lfc1_control)
            mat2 <- compute_lfc_matrix(input$lfc2_group_col, input$lfc2_treated, input$lfc2_control)

            # Find common peaks
            common_peaks <- intersect(rownames(mat1), rownames(mat2))
            if (length(common_peaks) == 0) {
              output$lfcHeatmap1 <- plotly::renderPlotly({
                ggplot2::ggplot() + ggplot2::theme_void() +
                  ggplot2::geom_text(ggplot2::aes(0.5,0.5,label="No common peaks between comparisons."), size=6)
              })
              output$lfcHeatmap2 <- plotly::renderPlotly({NULL})
              return()
            }
            # Combine matrices by common peaks
            mat1_sub <- mat1[common_peaks, , drop=FALSE]
            mat2_sub <- mat2[common_peaks, , drop=FALSE]
            combined_mat <- cbind(mat1_sub, mat2_sub)

            # Use user-selected peaks or top N by max absolute LFC
            user_peaks <- input$user_peaks
            if (!is.null(user_peaks) && length(user_peaks) > 0) {
              selected_peaks <- intersect(rownames(combined_mat), user_peaks)
            } else {
              n_peaks <- input$heatmap_npeaks
              lfc_max <- apply(abs(combined_mat), 1, max, na.rm=TRUE)
              selected_peaks <- names(sort(lfc_max, decreasing=TRUE))[1:min(n_peaks, length(lfc_max))]
            }
            combined_mat_top <- combined_mat[selected_peaks, , drop=FALSE]

            output$lfcHeatmap1 <- plotly::renderPlotly({
              if (is.null(combined_mat_top) || nrow(combined_mat_top) == 0) {
                return(ggplot2::ggplot() + ggplot2::theme_void() +
                         ggplot2::geom_text(ggplot2::aes(0.5,0.5,label="No data for selected groups/peaks."), size=6))
              }
              # Insert a gap column of NA values between the two sets
              n1 <- ncol(mat1_sub)
              n2 <- ncol(mat2_sub)
              gap_col <- matrix(NA, nrow = nrow(combined_mat_top), ncol = 1)
              colnames(gap_col) <- " "
              combined_with_gap <- cbind(
                combined_mat_top[, 1:n1, drop=FALSE],
                gap_col,
                combined_mat_top[, (n1+1):(n1+n2), drop=FALSE]
              )
              # Dynamically set colour scale limits centered on 0
              lfc_max <- max(abs(combined_with_gap), na.rm = TRUE)
              heatmaply::heatmaply(
                combined_with_gap,
                Rowv = TRUE,
                Colv = FALSE,
                colors = grDevices::colorRampPalette(c(scales::muted('blue'), 'white', scales::muted('red')))(100),
                xlab = "Comparison",
                ylab = "Peak",
                main = "Combined LFC Heatmap",
                dendrogram = "row",
                showticklabels = c(TRUE, TRUE),
                fontsize_row = 8,
                fontsize_col = 10,
                hide_colorbar = FALSE,
                limits = c(-lfc_max, lfc_max)
              )
            })

    })

      # --- Download handler for significance heatmap ---
  output$download_signif_heatmap <- shiny::downloadHandler(
    filename = function() input$signif_heatmap_filename,
    content = function(file) {
        res <- signif_heatmap_preproc()
        shiny::req(res)
        heatmaply::heatmaply(
          res$mat_z,
          Rowv = TRUE,
          Colv = TRUE,
          colors = grDevices::colorRampPalette(c(scales::muted('blue'), 'white', scales::muted('red')))(100),
          xlab = res$group_col,
          ylab = "Peak",
          main = paste0("Top ", input$heatmap_npeaks, " Significant Peaks Heatmap (", res$group_col, ")"),
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = scales::muted('blue'), high = scales::muted('red'), midpoint = 0, name = "Z-score"),
          dendrogram = "both",
          showticklabels = c(TRUE, TRUE),
          fontsize_row = 8,
          fontsize_col = 10,
          hide_colorbar = FALSE,
          file = file,
          selfcontained = TRUE,
          width = input$signif_heatmap_width*300,
          height = input$signif_heatmap_height*300
      )
    }
  )

  # --- Download handler for combined LFC heatmap ---
  output$download_lfc_heatmap <- shiny::downloadHandler(
    filename = function() input$lfc_heatmap_filename,
    content = function(file) {
      mat1 <- compute_lfc_matrix(input$lfc1_group_col, input$lfc1_treated, input$lfc1_control)
      mat2 <- compute_lfc_matrix(input$lfc2_group_col, input$lfc2_treated, input$lfc2_control)
      common_peaks <- intersect(rownames(mat1), rownames(mat2))
      if (length(common_peaks) == 0) stop("No common peaks between comparisons.")
      mat1_sub <- mat1[common_peaks, , drop=FALSE]
      mat2_sub <- mat2[common_peaks, , drop=FALSE]
      combined_mat <- cbind(mat1_sub, mat2_sub)
      user_peaks <- input$user_peaks
      if (!is.null(user_peaks) && length(user_peaks) > 0) {
        selected_peaks <- intersect(rownames(combined_mat), user_peaks)
      } else {
        n_peaks <- input$heatmap_npeaks
        lfc_max <- apply(abs(combined_mat), 1, max, na.rm=TRUE)
        selected_peaks <- names(sort(lfc_max, decreasing=TRUE))[1:min(n_peaks, length(lfc_max))]
      }
      combined_mat_top <- combined_mat[selected_peaks, , drop=FALSE]
      n1 <- ncol(mat1_sub)
      n2 <- ncol(mat2_sub)
      gap_col <- matrix(NA, nrow = nrow(combined_mat_top), ncol = 1)
      colnames(gap_col) <- " "
      combined_with_gap <- cbind(
        combined_mat_top[, 1:n1, drop=FALSE],
        gap_col,
        combined_mat_top[, (n1+1):(n1+n2), drop=FALSE]
      )
      # Dynamically set colour scale limits centered on 0
      lfc_max <- max(abs(combined_with_gap), na.rm = TRUE)
      heatmaply::heatmaply(
        combined_with_gap,
        Rowv = TRUE,
        Colv = FALSE,
        colors = grDevices::colorRampPalette(c(scales::muted('blue'), 'white', scales::muted('red')))(100),
        xlab = "Comparison",
        ylab = "Peak",
        main = "Combined LFC Heatmap",
        dendrogram = "row",
        showticklabels = c(TRUE, TRUE),
        fontsize_row = 8,
        fontsize_col = 10,
        hide_colorbar = FALSE,
        file = file,
        selfcontained = TRUE,
        width = input$lfc_heatmap_width*300,
        height = input$lfc_heatmap_height*300,
        limits = c(-lfc_max, lfc_max)
      )
    }
  )

  })

}
