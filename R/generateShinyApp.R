generateShinyApp <- function(shiny.dir='MSIToolKitApp',
                             expression.matrix,
                             metadata,
                             sample.id.column,
                             sample.wide.columns,
                             only.annotated=TRUE,
                             organism='Mouse',
                             adducts = c('M-H [1-]','M-H20-H [1-]','M+Cl [1-]'),
                             mode = 'Negative',
                             ppm = 5){
  # check inputs
  bulked <- create_bulk_exp(expression.matrix,
                            metadata,
                            sample.id.column = sample.id.column,
                            sample.wide.columns = sample.wide.columns,
                            organism = organism,
                            adducts = adducts,
                            mode = mode,
                            ppm = ppm)
  # show the user the assigned identities and ask if they want to include everything or just the annotated ones
  bulk.expression.matrix = bulked$expression_matrix
  bulk.metadata = bulked$metadata
  anno = bulked$annotation_table
  if (only.annotated){
    anno = anno[anno$adduct!='NA',]
    bulk.expression.matrix = bulk.expression.matrix[anno$m_z,]
    expression.matrix = as.data.frame(expression.matrix)[,anno$m_z]
  }
  return.list = c("expression.matrix",
                     "bulk.expression.matrix",
                     "metadata",
                     "bulk.metadata",
                     "anno")
  save(list=return.list,file=file.path(shiny.dir,'data.rda'))
  generateAppFile(shiny.dir,organism)
}


generateAppFile <- function(
    shiny.dir,
    organism
){
  lines.out <- c()
  
  packages.to.load <- c("MSIToolKit")
  code.load.packages <- paste0("library(", packages.to.load, ")")
  lines.out <- c(lines.out, code.load.packages, "")

  code.source.objects <- c(
    "rda.files <- list.files(pattern = '\\.rda$')",
    "for(fl in rda.files) load(fl)"
  )
  
  lines.out <- c(lines.out, code.source.objects, "")
  
  code.organism.set <- glue::glue("organism = '{organism}'")
  lines.out <- c(lines.out, code.organism.set, "")
  
  code.ui <- c(
    "ui <- function(request){",
    "navbarPage(",
    "'MSI toolkit',",
    "theme = shinythemes::shinytheme('flatly'),",
    "header = tags$head(tags$style('body {overflow-y: scroll;}')),",
    "footer = bookmarkButton(),",
    "tabPanel(title = 'Introduction',",
    "tabsetPanel(",
    "IntroAnnopanelUI(id='Anno', bulk.metadata = bulk.metadata),",
    "IntroSpatialVisPanelUI(id='spatialVis', bulk.metadata = bulk.metadata, full.metadata = metadata),",
    ")",
    "),",
    "tabPanel(title = 'Pseudobulk Analysis',",
    "tabsetPanel(",
    "BulkQCpanelUI(id='BulkQC', bulk.metadata = bulk.metadata),",
    "BulkDEpanelUI(id='BulkDE', bulk.metadata = bulk.metadata),",
    "BulkDESummaryPanelUI(id='BulkSummaryDE', bulk.metadata = bulk.metadata),",
    "BulkORAPanelUI(id='BulkORA', bulk.metadata = bulk.metadata)",
    ")",
    "),",
    "tabPanel(title = 'Region Analysis',",
    "tabsetPanel(",
    "RegionClusterPanelUI(id='RegionCluster', bulk.metadata = bulk.metadata, full.metadata = metadata),",
    "RegionDEpanelUI(id='RegionDE', full.metadata = metadata, bulk.metadata),",
    "BulkORAPanelUI(id='RegionORA', bulk.metadata = bulk.metadata),",
    "RegionNMFPanelUI(id='RegionNMF', bulk.metadata = bulk.metadata, full.metadata = metadata),",
    ")",
    ")",
    ")",
    "}"
  )

  lines.out <- c(lines.out, code.ui, "")
  
  code.server <- c(
    "server <- function(input, output, session) {",
    "IntroAnnopanelServer(id='Anno', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "IntroSpatialVisPanelServer(id='spatialVis', bulk.metadata = bulk.metadata, full.expression.matrix = expression.matrix, full.metadata = metadata, anno = anno)",
    "BulkQCpanelServer(id='BulkQC', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "bulkDEres <- BulkDEpanelServer(id='BulkDE', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "BulkDESummaryPanelServer(id='BulkSummaryDE', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = bulkDEres)",
    "BulkORAPanelServer(id='BulkORA', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = bulkDEres)",
    "clusters <- RegionClusterPanelServer(id='RegionCluster', full.expression.matrix = as.data.frame(t(expression.matrix)), full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno)",
    "regionDEres <- RegionDEpanelServer(id='RegionDE', full.expression.matrix = expression.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, region.clusters = clusters, anno = anno)",
    "BulkORAPanelServer(id='RegionORA', bulk.expression.matrix = bulk.expression.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = regionDEres)",
    "RegionNMFPanelServer(id='RegionNMF', full.expression.matrix = expression.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno)",
    "}"
  )
  lines.out <- c(lines.out, code.server, "")
  
  lines.out <- c(lines.out, "shinyApp(ui, server, enableBookmarking = 'url')")
  
  lines.out <- gsub("\\\\", "\\\\\\\\", lines.out)
  
  shiny.dir <- normalizePath(shiny.dir)
  write(lines.out, paste0(shiny.dir, "/app.R"))
  
}
