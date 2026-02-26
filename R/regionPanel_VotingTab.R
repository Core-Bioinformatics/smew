##' RegionVotingServer
##'
##' Server logic for the region-level voting scheme tab in the SMEW app.
##'
##' @param id Shiny module id
##' @param bulk.metadata Data frame of bulk sample metadata
##' @param full.metadata Data frame of full sample metadata
##' @param full.intensity.matrix Matrix of intensities (features x samples)
##' @param anno Data frame of peak annotations
##' @param shared_data Reactive or shared data object
##' @return None; called for side effects in Shiny module
##' @export
 RegionVotingServer <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno, shared_data){
   shiny::moduleServer(id, function(input, output, session){
     voting_result <- shiny::reactiveVal(NULL)

     shiny::updateSelectizeInput(session, "voting_peaks", "Select peaks:", choices = anno$display_name, server = TRUE)

     shiny::observeEvent(input$run_voting, {
       shiny::withProgress(message = 'Applying voting scheme...', value = 0, {
         shiny::req(input$voting_peaks, input$voting_percent, input$voting_min_peaks, input$samplesToFactor)
         selected_peaks <- input$voting_peaks
         percent <- input$voting_percent
         min_peaks <- input$voting_min_peaks
         sample_col <- colnames(bulk.metadata)[1]
         sample_mask <- full.metadata[, sample_col] %in% input[['samplesToFactor']]
         current.intensity.matrix <- t(full.intensity.matrix)[, sample_mask, drop=FALSE]
         current.metadata <- full.metadata[sample_mask, , drop=FALSE]
         peak_indices <- anno$m_z[anno$display_name %in% selected_peaks]
         vote_matrix <- matrix(0, nrow = nrow(current.metadata), ncol = length(selected_peaks))
         for(i in seq_along(peak_indices)){
           peak_id <- as.character(peak_indices[i])
           if (!peak_id %in% rownames(current.intensity.matrix)) {
             warning(sprintf("Peak %s not found in intensity matrix, skipping.", peak_id))
             next
           }
           intensities <- current.intensity.matrix[peak_id, ]
           threshold <- stats::quantile(intensities, probs = 1 - percent/100)
           vote_matrix[, i] <- as.integer(intensities >= threshold)
           shiny::incProgress(1/length(peak_indices))
         }
         votes <- rowSums(vote_matrix)
         mask <- votes >= min_peaks
         voting_result(mask)
       })
     })

    voting_plot <- shiny::reactive({
       shiny::req(voting_result(), input$samplesToFactor)
       mask <- ifelse(voting_result(),'Pass','Fail')
       sample_col <- colnames(bulk.metadata)[1]
       plot_df <- full.metadata[full.metadata[, sample_col] %in% input[['samplesToFactor']], , drop=FALSE]
       plot_df$voting_mask <- mask
       ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x_tf, y = .data$y_tf, fill = .data$voting_mask, color = .data$voting_mask)) +
         ggplot2::geom_tile() +
         ggplot2::facet_wrap(~.data$Sample, scales = 'free') +
         ggplot2::scale_fill_manual(values = c('Fail' = 'grey90', 'Pass' = 'red'), name = 'Region') +
         ggplot2::scale_color_manual(values = c('Fail' = 'grey90', 'Pass' = 'red'), name = 'Region') +
         ggplot2::theme_classic() +
         ggplot2::theme(axis.title = ggplot2::element_blank(),
                        axis.text = ggplot2::element_blank(),
                        axis.ticks = ggplot2::element_blank(),
                        axis.line = ggplot2::element_blank(),
                        aspect.ratio = 1)
    })
     output$voting_plot <- shiny::renderPlot({
      voting_plot()
     })

    shiny::observeEvent(input$add_to_shared, {
      shiny::req(voting_result(), input$voting_region_name)
      region_name <- input$voting_region_name
      if (is.null(region_name) || region_name == "") region_name <- "voting_region"
      meta <- shared_data$updated.metadata
      sample_col <- colnames(bulk.metadata)[1]
      sample_mask <- meta[, sample_col] %in% input[['samplesToFactor']]
      # Ensure the column exists
      if (!region_name %in% colnames(meta)) {
        meta[[region_name]] <- NA
      }
      mask <- voting_result()
      meta[[region_name]][sample_mask] <- ifelse(mask, 'Pass', 'Fail')
      shared_data$updated.metadata <- shiny::isolate(meta)
    })

      output[['downloadSpatialVoting']] <- create_download_plot_handler(
        plot_func = voting_plot,
        filename_func = function() input[['spatialVotingFileName']],
        width_func = function() input[['spatialVotingWidth']],
        height_func = function() input[['spatialVotingHeight']]
    )
  })
}

##' RegionVotingUI
##'
##' UI for the region-level voting scheme tab in the SMEW app.
##'
##' @param id Shiny module id
##' @param bulk.metadata Data frame of bulk sample metadata
##' @param full.metadata Data frame of full sample metadata
##' @param full.intensity.matrix Matrix of intensities (features x samples)
##' @param anno Data frame of peak annotations
##' @param show Logical; whether to show the panel (default TRUE)
##' @return A shiny tabPanel object for the voting scheme tab
##' @export
RegionVotingUI <- function(id, bulk.metadata, full.metadata, full.intensity.matrix, anno, show = TRUE){
  ns <- shiny::NS(id)
  if(show){
    shiny::tabPanel(
      'Voting Scheme',
      shiny::br(),
      bslib::accordion(
        bslib::accordion_panel(
          title = "Information",
          icon = bsicons::bs_icon("info-circle"),
          open = FALSE,
          shiny::tags$ul(
            shiny::tags$li("This tab allows you to create regions of interest using a voting scheme across multiple peaks."),
            shiny::tags$li("Select multiple peaks, set the top x% threshold for each, and specify the minimum number of peaks a pixel must pass to be included."),
            shiny::tags$li("The resulting region can be added to shared_data for downstream analysis.")
          )
        ),
        bslib::accordion_panel(
          title = "Sample selection",
          icon = bsicons::bs_icon("gear"),
          shiny::selectInput(
            inputId = ns("samplesToFactor"),
            label = "Select samples to include",
            choices = unique(bulk.metadata[, 1]),
            selected = unique(bulk.metadata[, 1]),
            multiple = TRUE
          )
        ),
        bslib::accordion_panel(
        title = "Downloads",
        icon = bsicons::bs_icon("download"),
        shiny::tags$div(
          shiny::tags$h4("Download output plots and tables"),
          shiny::fluidRow(
              shiny::tags$strong("Spatial distance plot"),
              shiny::textInput(ns('spatialVotingFileName'), 'File name', value ='SpatialVotingScheme.png'),
              shiny::numericInput(ns('spatialVotingWidth'),value = 8,label = 'Width (in)',min = 1,max = 50,step = 1),
              shiny::numericInput(ns('spatialVotingHeight'),value = 4,label = 'Height (in)',min = 1,max = 50,step = 1),
              shiny::downloadButton(ns('downloadSpatialVoting'), 'Download spatial distance')
          )
        )
        ),
        id = ns("acc"),
        open = "Sample selection"
      ),
      shiny::fluidRow(
        shiny::column(4,
          shiny::selectizeInput(ns("voting_peaks"), "Select peaks:", choices = c(), multiple = TRUE),
          shiny::sliderInput(ns("voting_percent"), "Top x% of pixels per peak:", min = 1, max = 50, value = 10),
          shiny::numericInput(ns("voting_min_peaks"), "Minimum number of peaks to pass:", value = 1, min = 1),
          shiny::actionButton(ns("run_voting"), "Apply Voting Scheme"),
          shiny::textInput(ns("voting_region_name"), "Name for region:", value = "VotingRegion"),
          shiny::actionButton(ns("add_to_shared"), "Create region")
        ),
        shiny::column(8,
          shiny::plotOutput(ns("voting_plot"), height = '600px')
        )
      )
    )
  } else {
    NULL
  }
}
