#' Write Example app.R Script
#'
#' This function writes a static app.R script using the provided data objects, matching the current example output of generateShinyApp.
#'
#' @param app_dir Directory to write the app.R file
#' @param bulk_intensity_matrix Path to or object for bulk intensity matrix
#' @param bulk_metadata Path to or object for bulk metadata
#' @param intensity_matrix Path to or object for intensity matrix
#' @param metadata Path to or object for metadata
#' @param anno Path to or object for annotation data
#' @return Path to the written app.R file
write_static_app_r <- function(app_dir, run_enrichment = FALSE, organism = NULL, run_autocorrelation = FALSE, run_pixel_enrichment = FALSE, has_multi_modal = FALSE, multi_modal_path = NULL) {
  # Build UI and server code blocks based on parameters
  pseudobulk_panels <- c(
    "BulkQCUI(id = 'BulkQC', bulk.metadata = bulk.metadata)",
    "BulkDEUI(id = 'BulkDE', bulk.metadata = bulk.metadata)",
    "BulkDESummaryUI(id = 'BulkSummaryDE', bulk.metadata = bulk.metadata)",
    if (run_enrichment) "BulkORAUI(id = 'BulkORA', bulk.metadata = bulk.metadata)",
    "BulkGRNUI(id = 'BulkGRN', bulk.metadata = bulk.metadata)",
    "BulkMultiCompareUI(id = 'BulkMultiCompare', bulk.metadata = bulk.metadata, bulk.intensity.matrix = bulk.intensity.matrix)",
    if (has_multi_modal) "BulkMultiGRNUI(id = 'BulkMultiGRN', bulk.metadata = bulk.metadata)"
  )
  region_panels <- c(
    "RegionDimRedUI(id = 'DimRed', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    "RegionClusterUI(id = 'Cluster', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    "RegionVotingUI(id = 'Voting', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix, anno = anno, show = TRUE)",
    "RegionHistologyUI(id = 'HistoTab', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    "RegionSpatialClusterPanelUI(id = 'SpatialClusterPanel', bulk.metadata = bulk.metadata, full.metadata = metadata, show = TRUE)",
    "RegionSmoothingUI(id = 'Smoothing')",
    "RegionComparisonUI(id = 'Comparison')",
    "RegionClusterDEUI('ClusterDE', bulk.metadata = bulk.metadata, full.metadata = metadata)",
    "RegionRadialDistanceUI(id = 'Radial', bulk.metadata)",
    "RegionGRNUI(id = 'RegionGRN', bulk.metadata)"
  )
  pixel_panels <- c(
    if (run_autocorrelation | (run_pixel_enrichment & run_enrichment)) "shiny::tabPanel('Pixel-level Analysis',\nshiny::tabsetPanel(\n",
    if (run_autocorrelation) "PixelSVMUI(id = 'PixelSVM', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix)",
    if (run_autocorrelation & run_pixel_enrichment & run_enrichment) ",\n",
    if (run_pixel_enrichment & run_enrichment) "PixelEnrichmentUI(id = 'pixelEnrichment', bulk.metadata = bulk.metadata, full.metadata = metadata, pixel_enrichment = pixel_enrichment)",
    if (run_autocorrelation | (run_pixel_enrichment & run_enrichment)) "))\n"
  )

  pseudobulk_panels <- pseudobulk_panels[!sapply(pseudobulk_panels, is.null)]
  region_panels <- region_panels[!sapply(region_panels, is.null)]
  pixel_panels <- pixel_panels[!sapply(pixel_panels, is.null)]

  pseudobulk_server <- c(
    "BulkQCServer(id = 'BulkQC', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "bulkDEres <- BulkDEServer(id = 'BulkDE', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno)",
    "BulkDESummaryServer(id = 'BulkSummaryDE', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, anno = anno, DEresults = bulkDEres)",
    if (run_enrichment) "BulkORAServer(id = 'BulkORA', bulk.intensity.matrix = bulk.intensity.matrix, bulk.metadata = bulk.metadata, DEresults = bulkDEres, anno = anno, organism = organism)",
    "BulkGRNServer('BulkGRN', bulk.intensity.matrix, bulk.metadata, anno)",
    "BulkMultiCompareServer('BulkMultiCompare', bulk.intensity.matrix, bulk.metadata, anno)",
    if (has_multi_modal) "BulkMultiGRNServer('BulkMultiGRN', bulk.intensity.matrix, bulk.metadata, anno, multi_modal)"
  )
  region_server <- c(
    "shared_data <- reactiveValues(updated.metadata = metadata)",
    "RegionDimRedServer(id = 'DimRed', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)",
    "RegionClusterServer(id = 'Cluster', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)",
    "RegionVotingServer('Voting', bulk.metadata, metadata, intensity.matrix, anno, shared_data)",
    "RegionHistologyServer(id = 'HistoTab', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)",
    "RegionSpatialClusterPanelServer(id = 'SpatialClusterPanel', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, anno = anno, shared_data = shared_data)",
    "RegionSmoothingServer('Smoothing', shared_data, metadata)",
    "RegionComparisonServer('Comparison', shared_data, bulk.metadata)",
    "RegionClusterDEServer('ClusterDE', full.intensity.matrix = intensity.matrix, full.metadata = metadata, bulk.metadata = bulk.metadata, shared_data = shared_data, anno = anno)",
    "RegionRadialDistanceServer(id = 'Radial', intensity.matrix, metadata, bulk.metadata, shared_data, anno)",
    "RegionGRNServer(id = 'RegionGRN', intensity.matrix, metadata, anno, bulk.metadata, shared_data)"
  )
  pixel_server <- c(
    if (run_autocorrelation) "PixelSVMServer(id = 'PixelSVM', bulk.metadata = bulk.metadata, full.metadata = metadata, full.intensity.matrix = intensity.matrix, anno = anno, DEresults = bulkDEres, svm_identification = svm_identification, spatial.cross.cor = spatial.cross.cor)",
    if (run_pixel_enrichment & run_enrichment) "PixelEnrichmentServer(id = 'pixelEnrichment', bulk.metadata = bulk.metadata, full.intensity.matrix = intensity.matrix, full.metadata = metadata, anno = anno, pixel_enrichment = pixel_enrichment)"
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
title = tags$img(src='logo without text.png',width='70px'),
window_title = 'SMEW: Spatial Metabolomics Enhanced Workflow',
theme = bslib::bs_theme(bootswatch = 'flatly'),
shiny::tabPanel('Introduction',
shiny::tabsetPanel(
IntroOverviewUI(id = 'Overview'),
IntroAnnoUI(id = 'Anno', bulk.metadata = bulk.metadata),
IntroSpatialVisUI(id = 'spatialVis', bulk.metadata = bulk.metadata, full.metadata = metadata)
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
IntroOverviewServer(id = 'Overview')
IntroAnnoServer(id = 'Anno',
bulk.intensity.matrix = bulk.intensity.matrix,
bulk.metadata = bulk.metadata, anno = anno)
IntroSpatialVisServer(id = 'spatialVis',
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

# ---- Preprocessing Function ----
##' Preprocess Data for App Generation
##'
##' This function runs all necessary preprocessing steps on the input data for the SMEW app.
##'
##' @param intensity_csv Path to the intensity matrix CSV file
##' @param metadata_csv Path to the metadata CSV file
##' @param output_dir Directory to save processed data
##' @param denoise Logical; whether to denoise the data (default FALSE)
##' @param anno Optional annotation data.frame or NULL
##' @param adducts Optional adducts table or NULL
##' @param ion_mode Ionization mode (character or NULL)
##' @param organism Organism name (default 'Human')
##' @param ppm Numeric; mass accuracy in ppm (default 10)
##' @param only_annotated Logical; if TRUE, only annotated peaks are used (default FALSE)
##' @param histology_images_dir Optional directory with histology images
##' @param run_autocorrelation Logical; whether to run spatial autocorrelation (default FALSE)
##' @param top_autocorrelated_peaks Integer; number of top autocorrelated peaks to use (default 10)
##' @param n_cores Number of cores for parallel processing (default 1)
##' @param run_pixel_enrichment Logical; whether to run pixel-level enrichment (default FALSE)
##' @return List or path to processed data
##' @export
preprocess_for_app <- function(intensity_csv, metadata_csv, output_dir, denoise = F, anno = NULL, 
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

  ## Here you would map m/z values to annotations and filter according to user preference
  if (!is.null(anno)) {
    message("Using user-supplied annotations")
    ## Make sure anno table is of the right format and filter to only annotated ones if wanted
    # intensity <- intensity |> select(Sample, matches(paste(anno, collapse = "|")))
    # Check that they are valid KEGG IDs and if not then set all enrichment to false and warn the user
    # STUFF TO DO HERE!!
  } else {
    message("No user-supplied annotations provided, using m/z values as feature names and mapping using internal annotation function")
    peak_list <- as.numeric(gsub('mz_', '', intensity_header))
    if (any(is.na(peak_list))) {
      message("Intensity header is not in expected format 'mz_<number>', using as is.")
      peak_list <- intensity_header
      anno <- data.frame(exp_peak = peak_list, adduct = NA, name = NA, kegg_id = NA, display_name = peak_list)
      enrichment_possible <- FALSE
  } else {
    message("Intensity header is in expected format, proceeding with m/z mapping.")
    enrichment_possible <- TRUE
    mapped_peaks <- map_peaks_to_kegg(peak_list, kegg_db, adducts, ion_mode,
                                      neg_adduct_table, pos_adduct_table, ppm)
    anno <- combine_peak_annotations(peak_list, mapped_peaks, fields = c('adduct','compound_name','compound_id'), sep = ', ', organism = organism)
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
  gcd = get_gcds(metadata)
  message("GCD values for grid spacing: ", paste(gcd, collapse = ", "))

  # 4. Transform coordinates to gcd 1 for app
  message("Transforming sample coordinates to gcd 1 for app...")
  metadata <- transform_sample_coordinates(metadata)

  # 5. Make bulk sample-level metadata and intensity matrix
  message("Creating bulk sample-level metadata and intensity matrix...")
  intensity <- data.table::fread(intensity_csv, select = (which(intensity_header %in% anno$m_z) + 1)) |> as.matrix()
  rownames(intensity) <- first_col_intensity
  bulked = create_bulk_matrix(metadata,intensity)

  # 5. Run denoising
  if (denoise) {
    message("Running denoising step...")
    processed_data <- denoise_all_samples(metadata, intensity)
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
        image_file <- file.path(image_folder, 'MSI_HE_aligned.jpg')
        coords <- metadata[metadata$Sample == sample_name, c('x', 'y', 'Sample')]
        if (file.exists(image_file)) {
          overlay_plot <- tryCatch({
            plot_pixel_overlay_on_image(coords, sample_name, image_file,alpha=0.5, point_size = 0.5)
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
      svm_results <- spatial_autocorrelation_pipeline(intensity, metadata)
      saveRDS(svm_results, file = file.path(output_dir, 'smew_app', 'spatial_autocorrelation.rds'))
      svm_peaks = base::Reduce(base::union,lapply(FUN = function(x)utils::head(x$peak,top_autocorrelated_peaks),X=svm_results))
      scc = spatial_cross_cor(svm_peaks, intensity, metadata, ncores = n_cores)
      saveRDS(scc, file = file.path(output_dir, 'smew_app', 'spatial_cross_correlation.rds'))
    }

    if (run_pixel_enrichment){
      if (!enrichment_possible) {
        message("Pixel-level enrichment analysis cannot be run because m/z values could not be mapped to annotations. Please check your input data and annotation parameters.")
      } else {
        message("Running pixel-leavel enrichment analysis step...")
        enrichment_results <- run_pixellevel_pipeline_parallel(
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
    write_static_app_r(
      app_dir = output_dir,
      run_enrichment = enrichment_possible,
      organism = organism,
      run_autocorrelation = run_autocorrelation,
      run_pixel_enrichment = run_pixel_enrichment,
      has_multi_modal = has_multi_modal,
      multi_modal_path = multi_modal_path
    )
}



