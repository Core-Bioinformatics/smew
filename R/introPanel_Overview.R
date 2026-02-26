##' SMEW App Overview UI
##'
##' UI for the SMEW app overview tab, with placeholders for manuscript and GitHub links.
##'
##' @param id Shiny module id
##' @return A shiny tabPanel object for the overview tab
##' @export
IntroOverviewUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tabPanel("Overview",
    shiny::div(
      style = "display: flex; align-items: center; justify-content: space-between; background: #f5f6fa; padding: 36px 48px 28px 48px; border-radius: 10px; margin-bottom: 32px; min-height: 140px;",
      shiny::imageOutput(ns("overview_logo"), height = "200px", width = "auto")
    ),

    shiny::div(
      style = "display: flex; flex-wrap: wrap; gap: 32px; align-items: flex-start; justify-content: flex-start;",
      # Left column: App Structure
      shiny::div(
        style = "flex: 2 1 400px; min-width: 340px; background: #e3f2fd; border-radius: 10px; padding: 24px 32px; border: 1px solid #90caf9;",
        shiny::h2("SMEW Application Overview"),
        shiny::p("Welcome to SMEW (Spatial Metabolomics Enhanced Workflow), a modular Shiny application for interactive analysis and visualisation of spatial metabolomics data."),
        shiny::h3("App Structure"),
        shiny::div(
          shiny::h4("Introduction"),
          shiny::tags$ul(
            shiny::tags$li(shiny::tags$b("Annotation:"), "Explore metabolite peak annotations and search for masses of interest."),
            shiny::tags$li(shiny::tags$b("Spatial Visualization:"), "Interactive spatial plots of metabolite peak intensities and metadata information.")
          ),
          shiny::h4("Pseudobulk Analysis"),
          shiny::tags$ul(
            shiny::tags$li(shiny::tags$b("Quality Control (QC):"), "Visualize and assess data quality at the whole sample level and compare intensity profiles, including PCA, PLS-DA, boxplots and barplots across experimental conditions. This tab also provides a platform to identify any outliers or batch effects in your data which can be factored into further downstream analysis."),
            shiny::tags$li(shiny::tags$b("Differential Analysis (DA):"), "Perform pairwise comparisons between experimental conditions using parametric and non-parametric tests with multiple testing correction."),
            shiny::tags$li(shiny::tags$b("DA Summary:"), "Summarize and visualise results from DA comparisons using heatmaps and volcano plots."),
            shiny::tags$li(shiny::tags$b("Pathway ORA (if available):"), "Test for enrichment of differentially abundant peaks in known pathways using over-representation analysis based on KEGG pathways."),
            shiny::tags$li(shiny::tags$b("Covariation Network Inference:"), "Infer regulatory networks from sample-level data using GENIE3, with interactive network visualisation and comparison between multiple networks."),
            shiny::tags$li(shiny::tags$b("Multi-modal Covariation Network Inference (if available):"), " Infer and compare regulatory networks across multiple sample groups using multi-modal data."),
            shiny::tags$li(shiny::tags$b("Bulk Multi-Comparison:"), " Compare multiple conditions using ANOVA or Kruskal-Wallis tests with post-hoc analysis and cross plots."),

          ),
          shiny::h4("Region Analysis"),
          shiny::tags$ul(
            shiny::tags$li(shiny::tags$b("Dimensionality Reduction:"), " Identify and visualize spatial patterns across multiple samples using PCA, NMF and/or UMAP. Clusters can be created by thresholding the resulting dimensionality reductions or using individual peak intensities."),
            shiny::tags$li(shiny::tags$b("Clustering:"), " Create spatial regions/clusters based on molecular profiles and visualise cluster assignments."),
            shiny::tags$li(shiny::tags$b("Histology Integration:"), " Overlay molecular data with histological images to draw manual regions of interest."),
            shiny::tags$li(shiny::tags$b("Voting Scheme:"), " Use thresholding on one or multiple metabolic features to create consensus regions or clusters."),
            shiny::tags$li(shiny::tags$b("Spatial Clustering:"), " Perform spatial-informed clustering using BayesSpace."),
            shiny::tags$li(shiny::tags$b("Spatial Smoothing:"), " Apply spatial smoothing to any regions identified through this app or outside to reduce noise and highlight spatially-resolved patterns."),
            shiny::tags$li(shiny::tags$b("Cross-Cluster Comparison:"), " Compare different clustering and region-identification options to understand the overlap between regions."),
            shiny::tags$li(shiny::tags$b("Region-Based Differential Analysis:"), " Compare multiple regions to find metabolite peaks driving regions, optionally taking into account experimental conditions"),
            shiny::tags$li(shiny::tags$b("Radial Distance Analysis:"), " Analyse molecular changes as a function of distance from a reference point or region to find spatially-refined patterns."),
            shiny::tags$li(shiny::tags$b("Region-Based Covariation Network Inference:"), " Infer covariation networks within specific tissue regions and compare regulatory relationships across regions.")
          ),
          shiny::h4("Pixel-level Analysis"),
          shiny::tags$ul(
            shiny::tags$li(shiny::tags$b("SVM Identification (if available):"), " Identify metabolic peaks with distinct spatial patterns and groups of peaks with common spatial patterns auto-correlation and cross-correlation metrics."),
            shiny::tags$li(shiny::tags$b("Pixel Enrichment (if available):"), " Perform spatially-informed pathway enrichment analysis to identify regions with distinct metabolic changes.")
          )
        ),
      ),
      # Right column: Documentation and Manuscript
      shiny::div(
        style = "flex: 1 1 260px; min-width: 260px; display: flex; flex-direction: column; gap: 24px;",
        shiny::div(
          style = "background: #fff3e0; border-radius: 10px; padding: 20px 32px; border: 1px solid #ffb74d;",
          shiny::h3("Links to Documentation and GitHub"),
          shiny::p("Links to the documentation and GitHub repository will be added here when public."),
          shiny::tags$div(id = ns("github-link-placeholder"), style = "margin-top: 10px;", "[GitHub link placeholder]")
        ),
        shiny::div(
          style = "background: #e8f5e9; border-radius: 10px; padding: 20px 32px; border: 1px solid #81c784;",
          shiny::h3("Manuscript"),
          shiny::tags$div(id = ns("manuscript-link-placeholder"), style = "margin-top: 10px;", "[Manuscript link placeholder]")
        )
      )
    )
  )
}

##' SMEW App Overview Server
##'
##' Server logic for the SMEW app overview tab, with placeholders for manuscript and GitHub links.
##'
##' @param id Shiny module id
##' @return A shiny tabPanel object for the overview tab
##' @export
IntroOverviewServer <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    
    output$overview_logo <- shiny::renderImage({
      list(
        src = file.path("figures", "logo banner.png"),
        contentType = "image/png",
        width = NULL,
        height = 200,
        alt = "SMEW Logo"
      )
    }, deleteFile = FALSE)
    
  })
}
