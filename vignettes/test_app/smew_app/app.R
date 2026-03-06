library(shiny)
library(smew)
bulked <- readRDS('processed_bulked.rds')
bulk.metadata = bulked$bulk_metadata
bulk.intensity.matrix = bulked$bulk_intensity
intensity.matrix <- readRDS('processed_data.rds')
metadata <- readRDS('processed_metadata.rds')
anno <- readRDS('processed_anno.rds')

ui <- function(request) {
bslib::page_navbar(
title = tags$img(src='logo.png',width='70px'),
window_title = 'SMEW: Spatial Metabolomics Enhanced Workflow',
theme = bslib::bs_theme(bootswatch = 'flatly'),
shiny::tabPanel('Introduction',
shiny::tabsetPanel(
IntroPanel_OverviewTabUI(id = 'Overview'),
IntroPanel_AnnoTabUI(id = 'Anno', bulk.metadata = bulk.metadata),
IntroPanel_SpatialVisTabUI(id = 'spatialVis', bulk.metadata = bulk.metadata, full.metadata = metadata)
)
),
shiny::tabPanel('Pseudobulk Analysis',
shiny::tabsetPanel(BulkPanel_QCTabUI(id = 'BulkQC', bulk.metadata = bulk.metadata),
BulkPanel_DATabUI(id = 'BulkDA', bulk.metadata = bulk.metadata),
BulkPanel_DASummaryTabUI(id = 'BulkSummaryDA', bulk.metadata = bulk.metadata),
BulkPanel_ORATabUI(id = 'BulkORA', bulk.metadata = bulk.metadata),
BulkPanel_GRNTabUI(id = 'BulkGRN', bulk.metadata = bulk.metadata),
BulkPanel_MultiCompareTabUI(id = 'BulkMultiCompare', bulk.metadata = bulk.metadata, bulk.intensity.matrix = bulk.intensity.matrix))
),
shiny::tabPanel('Region Analysis',
shiny::tabsetPanel(
RegionPanel_DimRedTabUI(id = 'DimRed', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix),
RegionPanel_ClusterTabUI(id = 'Cluster', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix),
RegionPanel_VotingTabUI(id = 'Voting', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix, anno = anno),
RegionPanel_HistologyTabUI(id = 'HistoTab', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix),
RegionPanel_SpatialClustersTabUI(id = 'SpatialClusterPanel', bulk.metadata = bulk.metadata, full.metadata = metadata),
RegionPanel_SpatialSmoothingTabUI(id = 'Smoothing'),
RegionPanel_ComparisonTabUI(id = 'Comparison'),
RegionPanel_DATabUI('ClusterDA', bulk.metadata = bulk.metadata, full.metadata = metadata),
RegionPanel_RadialDistanceTabUI(id = 'Radial', bulk.metadata),
RegionPanel_GRNTabUI(id = 'RegionGRN', bulk.metadata))
),
)
}

server <- function(input, output, session) {
thematic::thematic_shiny(
bg = '#FFFFFF',
fg = '#000000',
accent = '#007BFF',
font = 'sans'
)

# --- Introduction Panel Servers ---
IntroPanel_OverviewTabServer(id = 'Overview')
IntroPanel_AnnoTabServer(id = 'Anno',
bulk.intensity.matrix = bulk.intensity.matrix,
bulk.metadata = bulk.metadata, anno = anno)
IntroPanel_SpatialVisTabServer(id = 'spatialVis',
bulk.metadata = bulk.metadata,
full.intensity.matrix = intensity.matrix,
full.metadata = metadata, anno = anno)

# --- Pseudobulk Analysis Servers ---
BulkPanel_QCTabServer(id = 'BulkQC', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)
de_res <- BulkPanel_DATabServer(id = 'BulkDA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)
BulkPanel_DASummaryTabServer(id = 'BulkSummaryDA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno, de.results = de_res)
BulkPanel_ORATabServer(id = 'BulkORA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, de.results = de_res, anno = anno, organism = organism)
BulkPanel_GRNTabServer('BulkGRN', bulk.intensity.matrix, bulk.metadata, anno)
BulkPanel_MultiCompareTabServer('BulkMultiCompare', bulk.intensity.matrix, bulk.metadata, anno)

# --- Region-level Analysis Servers ---
shared_data <- reactiveValues(updated.metadata = metadata)
RegionPanel_DimRedTabServer(id = 'DimRed', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)
RegionPanel_ClusterTabServer(id = 'Cluster', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)
RegionPanel_VotingTabServer('Voting', bulk.metadata, metadata, intensity.matrix, anno, shared_data)
RegionPanel_HistologyTabServer(id = 'HistoTab', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)
RegionPanel_SpatialClustersTabServer(id = 'SpatialClusterPanel', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)
RegionPanel_SpatialSmoothingTabServer('Smoothing', shared_data, metadata)
RegionPanel_ComparisonTabServer('Comparison', shared_data, bulk.metadata)
RegionPanel_DATabServer('ClusterDA', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, shared_data = shared_data, anno = anno)
RegionPanel_RadialDistanceTabServer(id = 'Radial', intensity.matrix, metadata, bulk.metadata, shared_data, anno)
RegionPanel_GRNTabServer(id = 'RegionGRN', intensity.matrix, metadata, anno, bulk.metadata, shared_data)

# --- Pixel-level Analysis Servers ---

}

# --- Organism parameter for reference ---
organism <- "Mouse"

shiny::shinyApp(ui, server)

