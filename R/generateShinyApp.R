#' Generate all files required for an autonomous MSI toolkit shiny app
#' @description This function creates an app.R file and all required objects
#' to run the app in .rda format in the target directory. Basic checks are
#' performed to avoid problems with input data format. The app directory
#' is standalone and can be used on another platform, as long as the MSI toolkit
#' package is installed there.
#' @param shiny.dir directory to store the shiny app
#' @param intensity.matrix the intensity matrix, a dataframe where rows correspond to
#' pixels (named by pixel ID) and columns correspond to peaks (named by m/z value).
#' @param metadata a data frame containing metadata for each pixel contained
#' in the intensity.matrix; must contain at minimum four columns:
#' the first column must contain the row names of the intensity matrix, the pixel IDs,
#' the second and third columns should be named x and y and contain (x,y) coordinates for each pixel,
#' and another column should correspond to the sample IDs which each pixel comes from
#' @param sample.id.column name of the column in metadata containing sample ID information where
#' each pixel comes from
#' @param sample.wide.columns name of any other columns which describe sample-wide information, e.g.
#' treatment, disease etc
#' @param only.annotated boolean argument describing whether only peaks which can be mapped to metabolites
#' should be included in the app (default: TRUE)
#' @param organism organism name to be used for metabolite mapping and pathway analysis.
#' Options are Human, Mouse or Rat (default: Mouse)
#' @param mode Whether the data comes from positive or negative mode (default: negative)
#' @param adducts Vector of adducts which should be used for metabolite mapping. Options are  M-H [1-],M-2H [2-],
#' M-3H [3-],M-H2O-H [1-],M-H+O [1-],M+K-2H [1-],M+Na-2H [1-],M+Cl [1-],M+Cl37 [1-],M+FA-H [1-],M+Hac-H [1-] ,
#' M+Br [1-],M+Br81 [1-],M+TFA-H [1-],M+ACN-H [1-],M+HCOO [1-],M+CH3COO [1-],2M-H [1-],2M+FA-H [1-],2M+Hac-H [1-],3M-H [1-],
#' M(C13)-H [1-],M(S34)-H [1-],M(Cl37)-H [1-] for negative mode and  M [1+],M+H [1+],M+2H [2+],M+3H [3+],M+Na [1+],M+2Na [2+]
#' M+3Na [3+],M+H+Na [2+],M+H+2Na [3+],M+2H+Na [3+],M+2Na-H [1+],M+NaCl [1+],M+K [1+],M+H+K [2+],M+ACN+H [1+],M+ACN+2H [2+],
#' M+ACN+Na [1+],M+2ACN+2H [2+],M+3ACN+2H [2+],M+2ACN+H [1+],M+H2O+H [1+],M-H2O+H [1+],M-H4O2+H [1+],M-HCOOH+H [1+],M+HCOONa [1+],
#' M-HCOONa+H [1+],M+HCOOK [1+],M-HCOOK+H [1+],M-CO+H [1+],M-CO2+H [1+],M-C3H4O2+H [1+],M+CH3OH+H [1+],M-NH3+H [1+],M+H+NH4 [2+],
#' M+NH4 [1+],M+IsoProp+H [1+],M+IsoProp+Na+H [1+],M+2K+H [1+],M+DMSO+H [1+],2M+H [1+],2M+NH4 [1+],2M+Na [1+],2M+3H2O+2H [2+],
#' 2M+K [1+],2M+ACN+H [1+],2M+ACN+Na [1+],M(C13)+H [1+],M(C13)+2H [2+],M(C13)+3H [3+],M(S34)+H [1+],M(Cl37)+H [1+]  for positive mode
#' @param ppm ppm tolerance for metabolite mapping (default: 5)
#' @export
generateShinyApp <- function(shiny.dir='MSIToolKitApp',
                             intensity.matrix,
                             metadata,
                             sample.id.column,
                             sample.wide.columns,
                             only.annotated=TRUE,
                             organism='Mouse',
                             adducts = c('M-H [1-]','M-H20-H [1-]','M+Cl [1-]'),
                             mode = 'Negative',
                             ppm = 5){
  # check inputs
  bulked <- create_bulk_exp(intensity.matrix,
                            metadata,
                            sample.id.column = sample.id.column,
                            sample.wide.columns = sample.wide.columns,
                            organism = organism,
                            adducts = adducts,
                            mode = mode,
                            ppm = ppm)
  # show the user the assigned identities and ask if they want to include everything or just the annotated ones
  bulk.intensity.matrix = bulked$intensity_matrix
  bulk.metadata = bulked$metadata
  anno = bulked$annotation_table
  if (only.annotated){
    anno = anno[anno$adduct!='NA',]
    bulk.intensity.matrix = bulk.intensity.matrix[anno$m_z,]
    intensity.matrix = as.data.frame(intensity.matrix)[,anno$m_z]
  }
  intensity.matrix = as.matrix(intensity.matrix)
  intensity.matrix.t = t(intensity.matrix)
  intensity.matrix.t.unique = unique(intensity.matrix.t)
  intensity.matrix = t(intensity.matrix.t.unique)
  bulk.intensity.matrix = bulk.intensity.matrix[colnames(intensity.matrix),]
  anno = anno[anno$m_z %in% colnames(intensity.matrix)]

  svm_identification <- run_svm(intensity.matrix,metadata,bulk.metadata)

  return.list = c("intensity.matrix",
                     "bulk.intensity.matrix",
                     "metadata",
                     "bulk.metadata",
                     "anno")
  save(list=return.list,file=file.path(shiny.dir,'data.rda'))
  save("svm_identification",file=file.path(shiny.dir,'svm_identification.rda'))
  generateAppFile(shiny.dir,organism)
}


generateAppFile <- function(
    shiny.dir,
    organism
){
  lines.out <- c()

  packages.to.load <- c("smew","shiny")
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
    "RegionDimRedPanelUI(id='RegionNMF', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix),",
    ")",
    "),",
    "tabPanel(title = 'Pixel-level Analysis',",
    "PixelSVMPanelUI(id='PixelSVM',bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    ")",
    ")",
    "}"
  )

  lines.out <- c(lines.out, code.ui, "")

  code.server <- c(
    "server <- function(input, output, session) {",
    "IntroAnnopanelServer(id='Anno', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "IntroSpatialVisPanelServer(id='spatialVis', bulk.metadata = bulk.metadata, full.intensity.matrix = intensity.matrix, full.metadata = metadata, anno = anno)",
    "BulkQCpanelServer(id='BulkQC', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "bulkDEres <- BulkDEpanelServer(id='BulkDE', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "BulkDESummaryPanelServer(id='BulkSummaryDE', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = bulkDEres)",
    "BulkORAPanelServer(id='BulkORA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = bulkDEres)",
    "clusters <- RegionClusterPanelServer(id='RegionCluster', full.intensity.matrix = as.data.frame(t(intensity.matrix)), full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno)",
    "regionDEres <- RegionDEpanelServer(id='RegionDE', full.intensity.matrix = intensity.matrix, bulk.intensity.matrix = bulk.intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, region.clusters = clusters, anno = anno)",
    "BulkORAPanelServer(id='RegionORA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno,DEresults = regionDEres, organism = organism)",
    "RegionDimRedPanelServer(id='RegionNMF', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno)",
    "PixelSVMPanelServer(id='PixelSVM',bulk.metadata = bulk.metadata,full.metadata = metadata, full.intensity.matrix = intensity.matrix,anno = anno, DEresults = bulkDEres, svm_identification = svm_identification,spatial.cross.cor = spatial.cross.cor)",
    "}"
  )
  lines.out <- c(lines.out, code.server, "")

  lines.out <- c(lines.out, "shinyApp(ui, server, enableBookmarking = 'url')")

  lines.out <- gsub("\\\\", "\\\\\\\\", lines.out)

  shiny.dir <- normalizePath(shiny.dir)
  write(lines.out, paste0(shiny.dir, "/app.R"))

}
