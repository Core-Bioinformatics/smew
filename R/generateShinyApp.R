#' Writes app.R Script
#'
#' This function writes a static app.R script using the provided data objects, matching the current example output of generateShinyApp.
#'
#' @param app_dir Directory to write the app.R file
#' @param run_enrichment Logical; whether to include enrichment panel
#' @param organism Character; organism name
#' @param run_autocorrelation Logical; whether to include autocorrelation panel
#' @param run_pixel_enrichment Logical; whether to include pixel enrichment panel
#' @param has_multi_modal Logical; whether multi-modal data is present
#' @param multi_modal_path Path to multi-modal data file (optional)
#' @return Path to the written app.R file
#' @keywords internal
preprocessing_write_app_file <- function(app_dir, run_enrichment = FALSE, organism = NULL, run_autocorrelation = FALSE, run_pixel_enrichment = FALSE, has_multi_modal = FALSE, multi_modal_path = NULL) {
  # Build UI and server code blocks based on parameters
  pseudobulk_panels <- c(
    "BulkPanel_QCTabUI(id = 'BulkQC', bulk.metadata = bulk.metadata)",
    "BulkPanel_DATabUI(id = 'BulkDA', bulk.metadata = bulk.metadata)",
    "BulkPanel_DASummaryTabUI(id = 'BulkSummaryDA', bulk.metadata = bulk.metadata)",
    if (run_enrichment) "BulkPanel_ORATabUI(id = 'BulkORA', bulk.metadata = bulk.metadata)",
    "BulkPanel_GRNTabUI(id = 'BulkGRN', bulk.metadata = bulk.metadata)",
    "BulkPanel_MultiCompareTabUI(id = 'BulkMultiCompare', bulk.metadata = bulk.metadata, bulk.intensity.matrix = bulk.intensity.matrix)",
    if (has_multi_modal) "BulkPanel_MultiGRNTabUI(id = 'BulkMultiGRN', bulk.metadata = bulk.metadata)"
  )
  region_panels <- c(
    "RegionPanel_DimRedTabUI(id = 'DimRed', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    "RegionPanel_ClusterTabUI(id = 'Cluster', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    "RegionPanel_VotingTabUI(id = 'Voting', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix, anno = anno)",
    "RegionPanel_HistologyTabUI(id = 'HistoTab', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    "RegionPanel_SpatialClustersTabUI(id = 'SpatialClusterPanel', bulk.metadata = bulk.metadata, full.metadata = metadata)",
    "RegionPanel_SpatialSmoothingTabUI(id = 'Smoothing')",
    "RegionPanel_ComparisonTabUI(id = 'Comparison')",
    "RegionPanel_DATabUI('ClusterDA', bulk.metadata = bulk.metadata, full.metadata = metadata)",
    "RegionPanel_RadialDistanceTabUI(id = 'Radial', bulk.metadata)",
    "RegionPanel_GRNTabUI(id = 'RegionGRN', bulk.metadata)"
  )
  pixel_panels <- c(
    if (run_autocorrelation | (run_pixel_enrichment & run_enrichment)) "shiny::tabPanel('Pixel-level Analysis',\nshiny::tabsetPanel(\n",
    if (run_autocorrelation) "PixelPanel_SVMTabUI(id = 'PixelSVM', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    if (run_autocorrelation & run_pixel_enrichment & run_enrichment) ",\n",
    if (run_pixel_enrichment & run_enrichment) "PixelPanel_EnrichmentTabUI(id = 'pixelEnrichment', bulk.metadata = bulk.metadata, full.metadata = metadata, pixel.enrichment = pixel_enrichment)",
    if (run_autocorrelation | (run_pixel_enrichment & run_enrichment)) "))\n"
  )

  pseudobulk_panels <- pseudobulk_panels[!sapply(pseudobulk_panels, is.null)]
  region_panels <- region_panels[!sapply(region_panels, is.null)]
  pixel_panels <- pixel_panels[!sapply(pixel_panels, is.null)]

  pseudobulk_server <- c(
    "BulkPanel_QCTabServer(id = 'BulkQC', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "de_res <- BulkPanel_DATabServer(id = 'BulkDA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "BulkPanel_DASummaryTabServer(id = 'BulkSummaryDA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno, de.results = de_res)",
    if (run_enrichment) "BulkPanel_ORATabServer(id = 'BulkORA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, de.results = de_res, anno = anno, organism = organism)",
    "BulkPanel_GRNTabServer('BulkGRN', bulk.intensity.matrix, bulk.metadata, anno)",
    "BulkPanel_MultiCompareTabServer('BulkMultiCompare', bulk.intensity.matrix, bulk.metadata, anno)",
    if (has_multi_modal) "BulkPanel_MultiGRNTabServer('BulkMultiGRN', bulk.intensity.matrix, bulk.metadata, anno, multi_modal)"
  )
  region_server <- c(
    "shared_data <- reactiveValues(updated.metadata = metadata)",
    "RegionPanel_DimRedTabServer(id = 'DimRed', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)",
    "RegionPanel_ClusterTabServer(id = 'Cluster', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)",
    "RegionPanel_VotingTabServer('Voting', bulk.metadata, metadata, intensity.matrix, anno, shared_data)",
    "RegionPanel_HistologyTabServer(id = 'HistoTab', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)",
    "RegionPanel_SpatialClustersTabServer(id = 'SpatialClusterPanel', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)",
    "RegionPanel_SpatialSmoothingTabServer('Smoothing', shared_data, metadata)",
    "RegionPanel_ComparisonTabServer('Comparison', shared_data, bulk.metadata)",
    "RegionPanel_DATabServer('ClusterDA', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, shared_data = shared_data, anno = anno)",
    "RegionPanel_RadialDistanceTabServer(id = 'Radial', intensity.matrix, metadata, bulk.metadata, shared_data, anno)",
    "RegionPanel_GRNTabServer(id = 'RegionGRN', intensity.matrix, metadata, anno, bulk.metadata, shared_data)"
  )
  pixel_server <- c(
    if (run_autocorrelation) "PixelPanel_SVMTabServer(id = 'PixelSVM', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix, anno = anno, de.results = de_res, svm.identification = svm_identification, spatial.cross.cor = spatial.cross.cor)",
    if (run_pixel_enrichment & run_enrichment) "PixelPanel_EnrichmentTabServer(id = 'pixelEnrichment', bulk.metadata = bulk.metadata, full.intensity.matrix = intensity.matrix, full.metadata = metadata, anno = anno, pixel.enrichment = pixel_enrichment)"
  )
  pseudobulk_server <- pseudobulk_server[!sapply(pseudobulk_server, is.null)]
  region_server <- region_server[!sapply(region_server, is.null)]
  pixel_server <- pixel_server[!sapply(pixel_server, is.null)]

  organism_literal <- if (!is.null(organism)) paste0('"', organism, '"') else 'NULL'
  app_r_code <- paste0(
"library(shiny)
library(smew)
bulked <- readRDS('processed_bulked.rds')
bulk.metadata = bulked$bulk_metadata
bulk.intensity.matrix = bulked$bulk_intensity
intensity.matrix <- readRDS('processed_data.rds')
metadata <- readRDS('processed_metadata.rds')
anno <- readRDS('processed_anno.rds')
",
if (has_multi_modal) "multi_modal <- readRDS('processed_multi_modal.rds')\n" else "",
if (run_autocorrelation) "spatial.cross.cor <- readRDS('spatial_cross_correlation.rds')\nsvm_identification <- readRDS('spatial_autocorrelation.rds')\n" else "",
if (run_pixel_enrichment & run_enrichment) "pixel_enrichment <- readRDS('pixel_level_enrichment.rds')\n" else "",
"
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
shiny::tabsetPanel(",
paste(pseudobulk_panels, collapse = ',\n'),
")
),
shiny::tabPanel('Region Analysis',
shiny::tabsetPanel(\n",
paste(region_panels, collapse = ",\n"),
")
),
",
paste(pixel_panels, collapse = "\n"),
")
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
",
paste(pseudobulk_server, collapse = "\n"), "

# --- Region-level Analysis Servers ---
",
paste(region_server, collapse = "\n"),
"

# --- Pixel-level Analysis Servers ---
",
paste(pixel_server, collapse = "\n"), "
}

# --- Organism parameter for reference ---
organism <- ", organism_literal, "

shiny::shinyApp(ui, server)
"
)
  app_r_path <- file.path(app_dir, "smew_app", "app.R")
  writeLines(app_r_code, app_r_path)
  message('Static app.R written to ', app_r_path)
  return(app_r_path)
}

#' Preprocesses data and creates files for SMEW app
#'
#' This is the main user-facing function to generate a ready-to-use SMEW Shiny app from your spatial metabolomics data. It performs all necessary preprocessing, saves processed data, and writes a static Shiny app (app.R) in the specified output directory. The resulting app enables interactive analysis and visualisation of spatial metabolomics data at bulk, region, and pixel levels.
#'
#' @param intensity_csv Path to the intensity matrix CSV file. The first column must be pixel IDs, and the remaining columns are m/z features (named as 'mz_<number>' or similar).
#' @param metadata_csv Path to the metadata CSV file. Must contain columns 'pixel_id' (matching intensity matrix), 'x', 'y', and 'Sample'.
#' @param output_dir Directory to save the processed data and the generated app. Will be created if it does not exist.
#' @param denoise Logical; whether to denoise the data (default FALSE).
#' @param anno Optional annotation data.frame or NULL. If provided, should map m/z features to metabolite names and KEGG IDs.
#' @param adducts Optional vector of adducts used for annotation. Options include M-H [1-], M-2H [2-], M-3H [3-], M-H2O-H [1-], M-H+O [1-], M+K-2H [1-], M+Na-2H [1-], M+Cl [1-], M+Cl37 [1-], M+FA-H [1-], M+Hac-H [1-], M+Br [1-], M+Br81 [1-], M+TFA-H [1-], M+ACN-H [1-], M+HCOO [1-], M+CH3COO [1-], 2M-H [1-], 2M+FA-H [1-], 2M+Hac-H [1-], 3M-H [1-], M(C13)-H [1-], M(S34)-H [1-], M(Cl37)-H [1-] for negative mode; and M [1+], M+H [1+], M+2H [2+], M+3H [3+], M+Na [1+], M+2Na [2+], M+3Na [3+], M+H+Na [2+], M+H+2Na [3+], M+2H+Na [3+], M+2Na-H [1+], M+NaCl [1+], M+K [1+], M+H+K [2+], M+ACN+H [1+], M+ACN+2H [2+], M+ACN+Na [1+], M+2ACN+2H [2+], M+3ACN+2H [2+], M+2ACN+H [1+], M+H2O+H [1+], M-H2O+H [1+], M-H4O2+H [1+], M-HCOOH+H [1+], M+HCOONa [1+], M-HCOONa+H [1+], M+HCOOK [1+], M-HCOOK+H [1+], M-CO+H [1+], M-CO2+H [1+], M-C3H4O2+H [1+], M+CH3OH+H [1+], M-NH3+H [1+], M+H+NH4 [2+], M+NH4 [1+], M+IsoProp+H [1+], M+IsoProp+Na+H [1+], M+2K+H [1+], M+DMSO+H [1+], 2M+H [1+], 2M+NH4 [1+], 2M+Na [1+], 2M+3H2O+2H [2+], 2M+K [1+], 2M+ACN+H [1+], 2M+ACN+Na [1+], M(C13)+H [1+], M(C13)+2H [2+], M(C13)+3H [3+], M(S34)+H [1+], M(Cl37)+H [1+] for positive mode.
#' @param ion_mode Ionisation mode ('Negative' or 'Positive'). Used for annotation.
#' @param organism Organism name (default 'Human'). Used for annotation and pathway enrichment.
#' @param ppm Numeric; mass accuracy in ppm (default 10). Used for annotation.
#' @param only_annotated Logical; if TRUE, only annotated peaks are used (default FALSE).
#' @param histology_images_dir Optional directory with histology images for overlay and visualisation.
#' @param run_autocorrelation Logical (default: FALSE); whether to run spatial autocorrelation and SVM analysis.
#' @param top_autocorrelated_peaks Integer; number of top autocorrelated peaks to use (default 10).
#' @param n_cores Number of cores for parallel processing (default 1, where processes are run sequentially).
#' @param run_pixel_enrichment Logical (default: FALSE); whether to run pixel-level pathway enrichment.
#' @param multi_modal_path Optional path to multi-modal data file (for multi-omics integration).
#' @param enrichment_controls Optional vector of control samples for enrichment analysis.
#' @param enrichment_comparisons Optional vector of comparison samples for enrichment analysis.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Checks and validates input files and columns.
#'   \item Orders metadata to match the intensity matrix.
#'   \item Maps m/z features to metabolite annotations (if possible).
#'   \item Optionally filters to only annotated peaks.
#'   \item Computes grid spacing and transforms coordinates for spatial analysis.
#'   \item Creates bulk-level and pixel-level data objects.
#'   \item Optionally denoises the data.
#'   \item Saves all processed data as .rds files in the output directory.
#'   \item Optionally overlays pixel coordinates on histology images and copies images to the app directory.
#'   \item Optionally runs spatial autocorrelation and cross-correlation analysis.
#'   \item Optionally runs pixel-level pathway enrichment analysis.
#'   \item Optionally integrates multi-modal data (e.g., transcriptomics, proteomics etc).
#'   \item Writes a static app.R file for a fully functional SMEW Shiny app which can be shared with collaborators or the community.
#' }
#'
#' The generated app supports interactive analysis at multiple levels, including quality control, differential analysis, pathway enrichment, network inference, clustering, spatial visualisation, and more. See the package documentation and vignette for a full description of the app features and input requirements.
#'
#' @return Invisibly returns the path to the generated app.R file. All processed data and app files are saved in the output directory.
#' @examples 
#' create_smew_app(
#'   intensity_csv = system.file("extdata", "bleo_sub_intensity.csv", package = "smew"),
#'   metadata_csv = system.file("extdata", "bleo_sub_meta.csv", package = "smew"),
#'   output_dir = tempdir(),
#'   organism = 'Mouse',
#' )
#' unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
#' @export
create_smew_app <- function(intensity_csv, metadata_csv, output_dir, denoise = FALSE, anno = NULL,
                                adducts = NULL, ion_mode = NULL, organism = 'Human',
                                ppm = 10, only_annotated = FALSE, histology_images_dir = NULL,
                                run_autocorrelation = FALSE, top_autocorrelated_peaks = 10, n_cores = 1,
                                run_pixel_enrichment = FALSE, multi_modal_path = NULL,
                                enrichment_controls = NULL, enrichment_comparisons = NULL) {

  # Create output directory if it doesn't exist
  ifelse(!dir.exists(file.path(output_dir)), dir.create(file.path(output_dir)), FALSE)

  # ---- Input Checks ----
  # 1. Check files exist and are readable
  message("Checking input files exist and are readable...")
  if (!file.exists(intensity_csv)) stop(paste("Intensity matrix file not found:", intensity_csv))
  if (!file.exists(metadata_csv)) stop(paste("Metadata file not found:", metadata_csv))

  # 2. Read files
  metadata <- tryCatch(utils::read.csv(metadata_csv, check.names = FALSE), error = function(e) stop("Failed to read metadata CSV: ", e$message))
  tryCatch({
    data.table::fread(intensity_csv, nrows = 5)   # small structural check
    TRUE
    }, error = function(e) stop("Failed to read intensity CSV: ", e$message))
  first_col_intensity <- data.table::fread(intensity_csv, select = 1)[[1]]

  # 3. Check required columns
  message("Checking required columns in input files...")
  if (!"Sample" %in% colnames(metadata)) stop("Metadata CSV must contain a 'Sample' column.")
  if (!"pixel_id" %in% colnames(metadata)) stop("Metadata CSV must contain a 'pixel_id' column.")
  if (!setequal(first_col_intensity, metadata$pixel_id)) stop("Intensity matrix rownames must match metadata pixel_id.")

  # 4. Check pixel_id overlap
  message("Checking pixel_id overlap between intensity matrix and metadata...")
  missing_pixels <- setdiff(first_col_intensity, metadata$pixel_id)
  if (length(missing_pixels) > 0) {
    stop(paste("The following pixel_id(s) in the intensity matrix are missing from metadata:", paste(utils::head(missing_pixels, 10), collapse=", "), if(length(missing_pixels)>10) "..." else ""))
  }

  # ----- Preprocessing Steps ----
  # 1. Ensure metadata is ordered to match intensity matrix
  message("Ordering metadata to match intensity matrix...")
  metadata <- metadata |>
    dplyr::filter(.data$pixel_id %in% first_col_intensity) |>
    dplyr::arrange(match(.data$pixel_id, first_col_intensity))

  # 2. Map m/z to feature names if not already provided
  message("Mapping m/z values to feature names...")
  intensity_header <- names(data.table::fread(intensity_csv, nrows = 0))[-1]

  # Here you would map m/z values to annotations and filter according to user preference
  if (!is.null(anno)) {
    message("Using user-supplied annotations")
    # Check that anno is a data.frame and has required columns
    required_anno_cols <- c("m_z", "adduct", "name", "kegg_id", "display_name")
    if (!is.data.frame(anno)) {
      stop("The supplied annotation (anno) must be a data.frame.")
    }
    missing_cols <- setdiff(required_anno_cols, colnames(anno))
    if (length(missing_cols) > 0) {
      stop(paste0("The supplied annotation table is missing required column(s): ", paste(missing_cols, collapse=", "))) 
    }
    # Optionally filter to only annotated ones if wanted
    # Check that KEGG IDs are valid (non-empty, non-NA) for enrichment
    if (only_annotated) {
      anno <- anno[!is.na(anno$name) & anno$name != "", ]
      message(paste0("Filtered to only the ", nrow(anno), " annotated peaks in user-supplied annotation table."))
    }
    if (!any(!is.na(anno$kegg_id) & anno$kegg_id != "")) {
      warning("No valid KEGG IDs found in user-supplied annotation table. Pathway enrichment will be disabled.")
      enrichment_possible <- FALSE
    } else {
      enrichment_possible <- TRUE
    }
  } else {
    message("No user-supplied annotations provided, using m/z values as feature names and mapping using internal annotation function")
    peak_list <- as.numeric(gsub('mz_', '', intensity_header))
    # If adducts, ion_mode, or organism is NULL, skip annotation and set anno to NA
    if (is.null(adducts) || is.null(ion_mode) || is.null(organism)) {
      message("adducts, ion_mode, or organism is NULL; skipping annotation. anno will have NA names.")
      peak_list <- intensity_header
      anno <- data.frame(exp_peak = peak_list, adduct = NA, name = NA, kegg_id = NA, display_name = peak_list)
      enrichment_possible <- FALSE
    } else if (any(is.na(peak_list))) {
      message("Intensity header is not in expected format 'mz_<number>', using as is.")
      peak_list <- intensity_header
      anno <- data.frame(exp_peak = peak_list, adduct = NA, name = NA, kegg_id = NA, display_name = peak_list)
      enrichment_possible <- FALSE
    } else {
      message("Intensity header is in expected format, proceeding with m/z mapping.")
      enrichment_possible <- TRUE
      mapped_peaks <- preprocessing_map_peaks_to_kegg(peak_list, kegg_db, adducts, ion_mode,
                                        neg_adduct_table, pos_adduct_table, ppm)
      anno <- preprocessing_combine_peak_annotations(peak_list, mapped_peaks, fields = c('adduct','compound_name','compound_id'), sep = ', ', organism = organism)
      message("Mapped peaks to KEGG: ", nrow(anno[!is.na(anno$name),]), " annotations generated.")
    }
  }

  # Filter to only annotated peaks if requested and if annotations are available
  if (only_annotated) {
    message(paste0("Filtering to only the ", nrow(anno[!is.na(anno$name), ]), " annotated peaks..."))
    anno <- anno[!is.na(anno$name), ]
  } else {
    message(paste0("Keeping all ", nrow(anno), " peaks for downstream analysis."))
  }

  # 3. Get GCD values for grid spacing
  message("Computing GCD values for grid spacing...")
  gcd = preprocessing_get_gcds(metadata)
  message("GCD values for grid spacing: ", paste(gcd, collapse = ", "))

  # 4. Transform coordinates to gcd 1 for app
  message("Transforming sample coordinates to gcd 1 for app...")
  metadata <- preprocessing_transform_sample_coordinates(metadata)

  # 5. Make bulk sample-level metadata and intensity matrix
  message("Creating bulk sample-level metadata and intensity matrix...")
  intensity <- data.table::fread(intensity_csv, select = (which(intensity_header %in% anno$m_z) + 1)) |> as.matrix()
  rownames(intensity) <- first_col_intensity
  bulked = suppressWarnings(preprocessing_create_bulk_matrix(metadata,intensity))

  # 5. Run denoising
  if (denoise) {
    message("Running denoising step...")
    processed_data <- preprocessing_denoise_all_samples(metadata, intensity)
  } else {
    message("Skipping denoising step...")
    processed_data <- intensity
  }
  # 6. Save processed data
    # ---- Copy package images to app folders ----
    pkg_www <- system.file("app/www", package = "smew")
    app_www <- file.path(output_dir, "smew_app", "www")
    if (dir.exists(pkg_www)) {
      if (!dir.exists(app_www)) dir.create(app_www, recursive = TRUE)
      file.copy(list.files(pkg_www, full.names = TRUE), app_www, overwrite = TRUE, recursive = FALSE)
      message("Copied package www images to ", app_www)
    }
    pkg_fig <- system.file("app/figures", package = "smew")
    app_fig <- file.path(output_dir, "smew_app", "figures")
    if (dir.exists(pkg_fig)) {
      if (!dir.exists(app_fig)) dir.create(app_fig, recursive = TRUE)
      file.copy(list.files(pkg_fig, full.names = TRUE), app_fig, overwrite = TRUE, recursive = FALSE)
      message("Copied package figures to ", app_fig)
    }
  message("Saving processed data for app...")
  ifelse(!dir.exists(file.path(output_dir, 'smew_app')), dir.create(file.path(output_dir, 'smew_app')), FALSE)
  saveRDS(processed_data, file = file.path(output_dir, 'smew_app', 'processed_data.rds'))
  saveRDS(metadata, file = file.path(output_dir,'smew_app','processed_metadata.rds'))
  saveRDS(bulked, file = file.path(output_dir, 'smew_app', 'processed_bulked.rds'))
  saveRDS(anno, file = file.path(output_dir, 'smew_app', 'processed_anno.rds'))

  # 7. Check histology images overlay with coordinates if available
    if (!is.null(histology_images_dir)) {
      message('Checking and visualizing pixel overlays on sample images...')
      samples <- unique(metadata$Sample)
      ifelse(!dir.exists(file.path(output_dir, 'histology_overlays')), dir.create(file.path(output_dir, 'histology_overlays')), FALSE)

      for (sample_name in samples) {
        image_folder <- file.path(histology_images_dir, sample_name)
        image_file <- file.path(image_folder, 'histology_aligned.jpg')
        coords <- metadata[metadata$Sample == sample_name, c('x', 'y', 'Sample')]
        if (file.exists(image_file)) {
          overlay_plot <- tryCatch({
            preprocessing_plot_pixel_overlay_on_image(coords, sample_name, image_file,alpha=0.5, point_size = 0.5)
          }, error = function(e) {
            message(paste('Overlay failed for sample', sample_name, ':', e$message))
            NULL
          })
          if (!is.null(overlay_plot)) {
            overlay_file <- file.path(output_dir, 'histology_overlays', paste0('overlay_', sample_name, '.png'))
            ggplot2::ggsave(overlay_file, plot = overlay_plot, width = 6, height = 6)
            message(paste('Overlay saved for sample', sample_name, 'at', overlay_file))
          }
        } else {
          message(paste('No image or coordinates found for sample', sample_name))
        }
      }

      # ---- Copy entire histology_images_dir into smew_app/histology_images ----

      if (!dir.exists(file.path(output_dir, 'smew_app', 'images'))) dir.create(file.path(output_dir, 'smew_app', 'images'))
      message("Copying histology images over to app directory...")
      matching_images = intersect(list.files(histology_images_dir), bulked$bulk_metadata$Sample)
      new.dirs = setdiff(matching_images, list.files(file.path(output_dir, 'smew_app', 'images')))
      for (new.path in new.dirs) {
        if (!dir.exists(file.path(output_dir, 'smew_app', 'images', basename(new.path)))) dir.create(file.path(output_dir, 'smew_app', 'images', basename(new.path)))
      }
      for (image in matching_images){
        file.copy(file.path(histology_images_dir, image), file.path(output_dir, 'smew_app', 'images'), recursive = TRUE)
      }
    }

    # --- Multi-modal check before pixel enrichment and autocorrelation ---
    has_multi_modal <- FALSE
    if (!is.null(multi_modal_path)) {
      if (!file.exists(multi_modal_path)) {
        stop(paste("multi_modal_path does not exist:", multi_modal_path))
      }
      multi_modal <- tryCatch(
        {
          mm <- utils::read.csv(multi_modal_path, row.names = 1, check.names = FALSE)
          mm
        },
        error = function(e) stop("Failed to read multi-omics data: ", e$message)
      )
      # Check column names match bulk.metadata
      bulk_cols <- bulked$bulk_metadata$Sample
      mm_cols <- colnames(multi_modal)
      if (!setequal(bulk_cols, mm_cols)) {
        stop("Column names in multi_modal_path do not match bulk metadata.")
      }
      # Reorder columns to match bulk.metadata
      multi_modal <- multi_modal[, bulk_cols, drop = FALSE]
      message("multi_modal_path loaded and columns reordered to match bulk metadata.")
      has_multi_modal <- TRUE
      saveRDS(multi_modal, file = file.path(output_dir, 'smew_app', 'processed_multi_modal.rds'))
    }

    # 8. Run SVM step is requested
    if (run_autocorrelation) {
      message("Running spatial autocorrelation step...")
      svm_results <- preprocessing_spatial_autocorrelation_pipeline(intensity, metadata)
      saveRDS(svm_results, file = file.path(output_dir, 'smew_app', 'spatial_autocorrelation.rds'))
      svm_peaks = base::Reduce(base::union,lapply(FUN = function(x)utils::head(x$peak,top_autocorrelated_peaks),X=svm_results))
      scc = preprocessing_spatial_cross_cor(svm_peaks, intensity, metadata, ncores = n_cores)
      saveRDS(scc, file = file.path(output_dir, 'smew_app', 'spatial_cross_correlation.rds'))
    }

    if (run_pixel_enrichment){
      if (!enrichment_possible) {
        message("Pixel-level enrichment analysis cannot be run because m/z values could not be mapped to annotations. Please check your input data and annotation parameters.")
      } else {
        message("Running pixel-leavel enrichment analysis step...")
        enrichment_results <- preprocessing_run_pixel_enrichment(
          intensity.matrix = intensity,
          bulk.intensity.matrix = bulked$bulk_intensity,
          metadata = metadata,
          anno = anno,
          control.samples = enrichment_controls, # Define control samples as needed
          comparison.samples = enrichment_comparisons, # Define comparison samples as needed
          ncores = n_cores,
          organism = organism
        )
        saveRDS(enrichment_results, file = file.path(output_dir, 'smew_app', 'pixel_level_enrichment.rds'))
      }
    }

    # ---- Write app.R after preprocessing ----
    preprocessing_write_app_file(
      app_dir = output_dir,
      run_enrichment = enrichment_possible,
      organism = organism,
      run_autocorrelation = run_autocorrelation,
      run_pixel_enrichment = run_pixel_enrichment,
      has_multi_modal = has_multi_modal,
      multi_modal_path = multi_modal_path
    )
}



