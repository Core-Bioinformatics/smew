utils::globalVariables(c("complete_compound_mass", "mw", "exp_peak", "allowed_tolerance"))

##' Transform sample coordinates to normalized grid
##'
##' For each sample, shifts x and y so the minimum is 0 and scales so the grid cell distance (gcd) is 1.
##' Adds x_tf and y_tf columns to the metadata.
##'
##' @param metadata Data frame with columns Sample, x, y, and spot_id.
##' @return Data frame with new columns x_tf and y_tf, normalized per sample.
##' @keywords internal
preprocessing_transform_sample_coordinates <- function(metadata) {
  metadata$x_tf <- NA
  metadata$y_tf <- NA
  gcd = preprocessing_get_gcds(metadata)
  for (sample in unique(metadata$Sample)) {
    idx <- which(metadata$Sample == sample)
    x = metadata$x[idx]
    y = metadata$y[idx]
    x0 <- min(x)
    y0 <- min(y)
    if (is.na(gcd) || gcd == 0) gcd <- 1
    metadata$x_tf[idx] <- round((x - x0) / gcd) + 1
    metadata$y_tf[idx] <- round((y - y0) / gcd) + 1
  }
  return(metadata)
}

##' Compute adduct m/z values for molecular weights
##'
##' Calculates theoretical m/z values for different ionization adducts based on molecular weights and adduct formulas.
##'
##' @param adduct_table Data frame with columns 'Ion_Name' (adduct names) and 'Ion_Mass' (adduct formulas as character expressions, using "PROTON" for H+).
##' @param mw Numeric vector of molecular weights to calculate adducts for.
##' @return Data frame with molecular weight in first column and calculated adduct m/z values in subsequent columns (one per adduct).
##' @keywords internal
preprocessing_compute_adduct_weights <- function(adduct_table, mw) {
  # Extract adduct information
  ion_name <- adduct_table$Ion_Name
  ion_mass <- adduct_table$Ion_Mass

  # Convert to list and evaluate formulas
  mass_list <- as.list(ion_mass)
  mass_user <- lapply(mass_list, function(x) {
    eval(parse(text = paste(gsub("PROTON", 1.00727646677, x))))
  })

  mass_df <- as.data.frame(mass_user)
  colnames(mass_df) <- ion_name
  mass_df <- cbind(mw, mass_df)

  return(mass_df)
}

##' Match experimental m/z values to theoretical adducts
##'
##' Matches experimental m/z values against theoretical masses with adduct modifications, returning compounds within specified ppm tolerance.
##'
##' @param adducts_table Data frame with molecular weights and adduct m/z values (output of preprocessing_compute_adduct_weights).
##' @param mz_list Numeric vector of experimental m/z values to match.
##' @param allowed_dppm Numeric; mass tolerance in parts per million (ppm).
##' @return Data frame with columns: exp_peak, allowed_tolerance, observed_difference, theoretical_mass, adduct, theoretical_mass_w_adduct, ppm_error.
##' @keywords internal
preprocessing_get_matched_comps <- function(adducts_table, mz_list, allowed_dppm) {
  # Calculate tolerance window for each peak
  tolerance <- mz_list * allowed_dppm * 1e-06

  # Get unique compounds
  adducts_dedup <- dplyr::distinct(adducts_table)

  # Convert to long format
  long_adducts_table <- adducts_dedup |>
    tidyr::pivot_longer(cols = -c(mw), names_to = 'adduct',
                        values_to = 'adduct_mass')

  # Create difference matrix (experimental vs theoretical)
  diff_df <- data.frame(outer(mz_list, long_adducts_table$adduct_mass, "-"))

  # Label columns with compound/adduct information
  colnames(diff_df) <- paste(long_adducts_table$mw, long_adducts_table$adduct,
                              long_adducts_table$adduct_mass, sep = '_')

  # Add peak and tolerance information
  diff_df$exp_peak <- mz_list
  diff_df$allowed_tolerance <- tolerance

  # Convert to long format
  diff_df_long <- tidyr::pivot_longer(diff_df,
                                      cols = -c(exp_peak, allowed_tolerance),
                                      names_to = 'mw_adduct_adductmass',
                                      values_to = 'observed_difference')

  colnames(diff_df_long) <- c("exp_peak", "allowed_tolerance",
                               "mw_adduct_adductmass", "observed_difference")

  # Calculate absolute mass difference
  diff_df_long$abs_difference <- abs(diff_df_long$observed_difference)

  # Filter to matches within tolerance
  diff_df_long$match <- diff_df_long$abs_difference < diff_df_long$allowed_tolerance
  diff_df_long <- diff_df_long[diff_df_long$match, ]

  # Parse compound information
  diff_df_long <- tidyr::separate_wider_delim(
    diff_df_long,
    cols = "mw_adduct_adductmass",
    names = c('theoretical_mass', 'adduct', 'theoretical_mass_w_adduct'),
    delim = '_'
  )

  # Calculate ppm error (parts per million error)
  # ppm = (observed_difference / exp_peak) * 1e6
  diff_df_long$ppm_error <- (diff_df_long$observed_difference / diff_df_long$exp_peak) * 1e6

  # Select and format output columns
  diff_df_long <- diff_df_long[, c('exp_peak', 'allowed_tolerance',
                                    'observed_difference', 'theoretical_mass',
                                    'adduct', 'theoretical_mass_w_adduct',
                                    'ppm_error')]

  # Convert to numeric and apply significance rounding
  diff_df_long <- diff_df_long |>
    dplyr::mutate(
      "exp_peak" = as.numeric(.data$exp_peak),
      "allowed_tolerance" = signif(as.numeric(.data$allowed_tolerance), 7),
      "observed_difference" = signif(as.numeric(.data$observed_difference), 7),
      "theoretical_mass" = as.numeric(.data$theoretical_mass),
      "theoretical_mass_w_adduct" = as.numeric(.data$theoretical_mass_w_adduct),
      "ppm_error" = signif(as.numeric(.data$ppm_error), 7)
    )

  return(diff_df_long)
}

##' Match experimental peaks to known compounds
##'
##' Matches experimental m/z peaks against a metabolite database using specified ionization mode and adduct formulas. Validates all input parameters.
##'
##' @param metabolite_db Data frame with known compound information. Expected columns are MetaboliteID, ExactMass, and MetaboliteName.
##' @param peak_list Numeric vector of experimental m/z values.
##' @param ppm Numeric; mass tolerance in parts per million (default: 5).
##' @param ion_mode Character; ionisation mode - either "Positive" or "Negative".
##' @param adducts Character vector of adduct names to use.
##' @param neg_adduct_formulas,pos_adduct_formulas Data frames with negative/positive ionisation adduct formulas.
##' @return Data frame of matched peaks with compound annotations and adduct information.
##' @keywords internal
preprocessing_get_matched_peaks <- function(metabolite_db = NULL, peak_list = NULL,
                               ppm = 5, ion_mode = NULL,
                               adducts = NULL, neg_adduct_formulas = NULL,
                               pos_adduct_formulas = NULL) {
  # ========================================================================
  # INPUT VALIDATION
  # ========================================================================
  if (is.null(metabolite_db)) {
    stop('Provide metabolite_db')
  }
  if (is.null(peak_list)) {
    stop('Provide peak_list (experimental m/z values)')
  }
  if (is.null(ppm)) {
    stop('Provide ppm (mass tolerance)')
  }
  if (is.null(ion_mode)) {
    stop('Provide ion_mode (ionisation mode)')
  }
  if (!ion_mode %in% c('Negative', 'Positive')) {
    stop('Ionisation mode must be either "Positive" or "Negative".')
  }

  required_cols <- c('MetaboliteID', 'ExactMass', 'MetaboliteName')
  missing_cols <- setdiff(required_cols, colnames(metabolite_db))
  if (length(missing_cols) > 0) {
    stop('metabolite_db is missing required columns: ', paste(missing_cols, collapse = ', '))
  }

  # Standardise to internal fields used by matching functions.
  metabolite_db$compound_id <- as.character(metabolite_db$MetaboliteID)
  metabolite_db$compound_name <- as.character(metabolite_db$MetaboliteName)
  metabolite_db$complete_compound_mass <- suppressWarnings(as.numeric(metabolite_db$ExactMass))

  exp_peak_list = unique(as.numeric(peak_list))
  # select adducts of interest
  if (ion_mode == 'Negative') {
    my_adduct_formulas = neg_adduct_formulas[neg_adduct_formulas$Ion_Name %in% adducts,]

  } else if (ion_mode == 'Positive') {
    my_adduct_formulas = pos_adduct_formulas[pos_adduct_formulas$Ion_Name %in% adducts,]
  }

  filt_metabolite_db = metabolite_db

  # drop rows without given mass
  filt_metabolite_db = filt_metabolite_db |> tidyr::drop_na(complete_compound_mass)

  # compute theoretical masses for all compounds+adducts in the kegg dataset
  adducts_computed = preprocessing_compute_adduct_weights(my_adduct_formulas,filt_metabolite_db$complete_compound_mass)

  # match experimental peaks to calculated adducts
  matched_comps = preprocessing_get_matched_comps(adducts_table = adducts_computed,
                                    mz_list = exp_peak_list,
                                    allowed_dppm = ppm)

  # merge with metabolite data for complete annotation (join on theoretical mass)
  matched_comps_annot = dplyr::left_join(matched_comps, filt_metabolite_db, by = c('theoretical_mass' = 'complete_compound_mass'), relationship = 'many-to-many')

  # convert theoretical_mass to numeric
  matched_comps_annot$theoretical_mass = as.numeric(matched_comps_annot$theoretical_mass)

  # add mode column
  matched_comps_annot$ion_mode = ion_mode

  # Keep core output columns for downstream annotation.
  col_order <- c(
    'exp_peak',
    'allowed_tolerance',
    'observed_difference',
    'theoretical_mass_w_adduct',
    'theoretical_mass',
    'adduct',
    'ppm_error',
    'ion_mode',
    'compound_id',
    'compound_name'
  )
  col_order <- col_order[col_order %in% colnames(matched_comps_annot)]
  matched_comps_annot <- matched_comps_annot[, col_order]

  return(matched_comps_annot)
}

##' Greatest common divisor for grid spacing
##'
##' Computes the greatest common divisor for a numeric vector, used for grid spacing.
##'
##' @param x Numeric vector.
##' @return Integer; greatest common divisor.
##' @keywords internal
preprocessing_calculate_gcd <- function(x) {
  m = min(x)

  while (any(x %% m > 0)){
    m = m - 1
  }

  return(m)
}

##' Compute grid spacing (GCD) for each sample
##'
##' Calculates the grid spacing (GCD) for each sample in the metadata, ensuring consistency across samples.
##'
##' @param metadata Data frame with 'Sample', 'spot_id', 'x', and 'y' columns.
##' @return Numeric; unique GCD value for the grid spacing.
##' @keywords internal
preprocessing_get_gcds <- function(metadata){
  gcd.values = list()
  for (sample in unique(metadata$Sample)){
    metadata.sub = metadata |> dplyr::filter(Sample==sample)
    pos <- metadata.sub[,c('x','y')]
    rownames(pos)=metadata.sub$spot_id
    rownames(metadata.sub)=metadata.sub$spot_id
    # calculating spacing between points
    min.x = min(pos$x)
    min.y = min(pos$y)
    pos$x = pos$x - min.x
    pos$y = pos$y - min.y
    gcd.x = preprocessing_calculate_gcd(round(unique(pos$x)[unique(pos$x)!=0]))
    gcd.y = preprocessing_calculate_gcd(round(unique(pos$y)[unique(pos$y)!=0]))
    if (gcd.x!=gcd.y){
      message('Mismatching gap!')
    } else {
      gcd.values[[sample]]<-gcd.x
    }
  }
  if (length(unique(gcd.values))!=1){
    stop('Mismatching gap across samples!')
  }
  return(unique(unlist(gcd.values)))
}

##' Denoise pixel intensity matrix for a sample
##'
##' Applies spatial and correlation-based denoising to pixel-wise intensity data for a given sample.
##'
##' @param metadata Data frame containing pixel metadata, including 'Sample', 'pixel_id', 'x', and 'y'.
##' @param intensity.matrix Data frame of intensity values (pixels as rows, features as columns).
##' @param sample Character or factor specifying the sample/group to denoise.
##' @param ncores Integer; number of cores for parallelisation (default 1).
##' @return Data frame of denoised intensities for the selected sample (pixels x features).
##' @keywords internal
preprocessing_denoise.intensity <- function(metadata, intensity.matrix, sample, ncores = 1) {
  # Subset metadata and intensity matrix for the sample
  current.metadata <- metadata |> dplyr::filter(Sample == sample)
  rownames(current.metadata) <- current.metadata$pixel_id

  current.intensity <- intensity.matrix[current.metadata$pixel_id, ]

  # Remove non-varying features
  current.intensity_pca <- current.intensity[, !apply(current.intensity, 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]

  # Create SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = t(current.intensity_pca)),
    colData = current.metadata
  )
  # Run PCA using scater
  sce <- scater::runPCA(sce, exprs_values = "counts", ncomponents = 50, scale = TRUE)

  pca_mat <- SingleCellExperiment::reducedDim(sce, "PCA")

  # Correlation matrix on PCA
  cor.matrix <- stats::cor(t(pca_mat[, 1:50]))

  # Spatial distance matrix
  distance.matrix <- as.matrix(stats::dist(current.metadata[, c('x_tf', 'y_tf')] ))

  # Weight matrix
  w <- cor.matrix * distance.matrix

  # Denoise by neighbourhood

  all.intensities <- pbapply::pbsapply(
    rownames(current.intensity),
    simplify = TRUE,
    USE.NAMES = TRUE,
    FUN = function(x) preprocessing_pick.neighbours.per.pixel(x, cor.matrix, distance.matrix, current.intensity, 1, w)
  )
  all.intensities <- data.frame(t(all.intensities))
  return(all.intensities)
}

##' Denoise a single pixel by neighbourhood
##'
##' Computes denoised intensity for a pixel by weighted averaging of its most correlated spatial neighbours.
##'
##' @param i Pixel identifier (row name).
##' @param cor.matrix Correlation matrix between pixels.
##' @param distance.matrix Spatial distance matrix between pixels.
##' @param current.intensity Intensity matrix for current sample.
##' @param gcd Numeric; grid spacing value for the sample.
##' @param w Weight matrix (correlation * distance).
##' @return Numeric vector of denoised intensities for pixel i.
##' @keywords internal
preprocessing_pick.neighbours.per.pixel <- function(i,cor.matrix,distance.matrix,current.intensity,gcd,w){
  # Get all neighbour indices within 3*gcd (excluding self)
  dists <- distance.matrix[i, ]
  valid_idx <- which(dists < 3 * gcd & names(dists) != i)
  if (length(valid_idx) == 0) {
    # No neighbours, return original intensity
    out <- current.intensity[i,]
    if (!is.numeric(out)) warning(paste("Output for", i, "is not numeric!"))
    if (is.null(names(out))) names(out) <- colnames(current.intensity)
    if (length(out) != ncol(current.intensity)) warning(paste("Output for", i, "has wrong length:", length(out)))
    if (any(is.na(out))) warning(paste("Output for", i, "contains NA values."))
    return(unlist(out))
  }
  # Get top 3 most correlated neighbours among valid_idx
  cors <- cor.matrix[i, valid_idx]
  if (length(cors) > 3) {
    top_idx <- valid_idx[order(cors, decreasing = TRUE)[1:3]]
  } else {
    top_idx <- valid_idx[order(cors, decreasing = TRUE)]
  }
  # Vectorized weighted sum
  weights <- w[i, top_idx]
  weight.sum <- sum(weights)
  if (weight.sum == 0) weight.sum <- 1
  new.intensity <- 0.5 * current.intensity[i, ] +
    0.5 * colSums(sweep(current.intensity[top_idx, , drop = FALSE], 1, weights, '*')) / weight.sum
 return(unlist(new.intensity))
}

##' Denoise all samples and combine
##'
##' Runs preprocessing_denoise.intensity for each sample and combines the results.
##'
##' @param metadata Data frame with pixel metadata (must include Sample column).
##' @param intensity.matrix Full intensity matrix (pixels x features).
##' @return Data frame of combined denoised intensities (pixels x features).
##' @keywords internal
preprocessing_denoise_all_samples <- function(metadata, intensity.matrix) {
  samples <- unique(metadata$Sample)
  denoised_list <- vector("list", length(samples))
  names(denoised_list) <- samples
  for (i in seq_along(samples)) {
    sample <- samples[i]
    message(paste("Denoising sample:", sample))
    denoised_list[[i]] <- preprocessing_denoise.intensity(metadata, intensity.matrix, sample)
  }
  combined <- do.call(rbind, denoised_list)
  return(combined)
}

##' Create bulk matrix and metadata
##'
##' Aggregates intensity matrix per sample and extracts bulk metadata columns with only one value per sample.
##'
##' @param metadata Data frame with pixel metadata (must include Sample column).
##' @param intensity.matrix Full intensity matrix (pixels x features).
##' @param agg_fun Aggregation function (e.g., mean, median).
##' @return List with bulk_intensity (features x samples) and bulk_metadata (samples x columns).
##' @keywords internal
preprocessing_create_bulk_matrix <- function(metadata, intensity.matrix, agg_fun = mean) {
  samples <- unique(metadata$Sample)
  # Identify columns with only one value per sample
  single_value_cols <- sapply(metadata, function(col) {
    tapply(col, metadata$Sample, function(x) length(unique(x)) == 1)
  })
  # Keep columns where all samples have only one value
  bulk_cols <- names(which(colSums(single_value_cols) == length(samples)))
  bulk_metadata <- unique(metadata[, unique(c("Sample", bulk_cols)), drop = FALSE])
  # Aggregate intensities per sample
  agg_list <- lapply(samples, function(s) {
    idx <- which(metadata$Sample == s)
    if (length(idx) == 1) {
      as.numeric(intensity.matrix[idx, ])
    } else {
      apply(intensity.matrix[idx, , drop = FALSE], 2, agg_fun, na.rm = TRUE)
    }
  })
  bulk_intensity <- do.call(cbind, agg_list)
  colnames(bulk_intensity) <- samples
  rownames(bulk_intensity) <- colnames(intensity.matrix)
  return(list(bulk_intensity = bulk_intensity, bulk_metadata = bulk_metadata))
}

##' Map peak masses to user-supplied metabolite IDs
##'
##' Maps a vector of m/z peak values to a user-supplied metabolite table
##' (MetaboliteID, ExactMass, MetaboliteName) and selected adducts using
##' \\code{\\link{preprocessing_get_matched_peaks}}.
##'
##' @param peak_list Numeric vector of experimental m/z values.
##' @param metabolite_table Data frame with columns MetaboliteID, ExactMass, MetaboliteName.
##' @param adducts Character vector of adduct names to use.
##' @param ion_mode Character; ionization mode - either "Positive" or "Negative".
##' @param neg_adduct_formulas,pos_adduct_formulas Data frames with negative/positive adduct formulas.
##' @param ppm Numeric; mass tolerance in ppm (default: 5).
##' @return Data frame mapping each peak to metabolite IDs/names and adduct information.
##' @keywords internal
preprocessing_map_peaks_to_ids <- function(peak_list,
                                           metabolite_table,
                                           adducts,
                                           ion_mode,
                                           neg_adduct_formulas,
                                           pos_adduct_formulas,
                                           ppm = 5) {

  matched_peaks <- preprocessing_get_matched_peaks(
    metabolite_db = metabolite_table,
    peak_list = peak_list,
    ppm = ppm,
    ion_mode = ion_mode,
    adducts = adducts,
    neg_adduct_formulas = neg_adduct_formulas,
    pos_adduct_formulas = pos_adduct_formulas
  )

  out <- matched_peaks[, c(
    'exp_peak',
    'adduct',
    'compound_id',
    'compound_name',
    'theoretical_mass',
    'theoretical_mass_w_adduct',
    'observed_difference',
    'allowed_tolerance',
    'ppm_error',
    'ion_mode'
  )]

  # Provide explicit generic-ID aliases for downstream use.
  out$MetaboliteID <- out$compound_id
  out$MetaboliteName <- out$compound_name

  out <- unique(out)
  return(out)
}

##' Combine annotations for peak list
##'
##' For a given peak list and annotation table, create a combined annotation table in the order of the input peaks. Multiple annotations are pasted together, and NAs are included if no match exists.
##'
##' @param peak_list Numeric vector of original peaks (order preserved).
##' @param anno_table Data frame from preprocessing_map_peaks_to_ids.
##' @param fields Character vector of annotation columns to combine (default: c('adduct','compound_name','compound_id')).
##' @param sep Separator for multiple annotations (default: ', ').
##' @return Data frame with one row per peak, columns for each annotation field, and display_name.
##' @keywords internal
preprocessing_combine_peak_annotations <- function(peak_list, anno_table, fields = c('adduct','compound_name','compound_id'), sep = ', ') {
  anno_table = unique(anno_table[,c('exp_peak',fields,'ppm_error')])
  res <- lapply(peak_list, function(pk) {
    # Subset matches for this peak
    matches <- anno_table[anno_table$exp_peak == pk, , drop = FALSE]
    # Order by ppm_error (smallest first) if present
    matches <- matches[order(abs(matches$ppm_error), na.last = TRUE), ]
    if (nrow(matches) == 0) {
      stats::setNames(as.character(rep(NA, length(fields))), fields)
    } else {
      sapply(fields, function(f) {
        vals <- unique(stats::na.omit(matches[[f]]))
        if (length(vals) == 0) NA else paste(vals, collapse = sep)
      }, USE.NAMES = TRUE)
    }
  })
  res_df <- as.data.frame(do.call(rbind, res), stringsAsFactors = FALSE)
  colnames(res_df) <- c('adduct','name','metabolite_id')
  res_df <- cbind(m_z = paste0('mz_', peak_list), res_df)
  res_df$display_name <- ifelse(is.na(res_df$name), res_df$m_z, paste0(res_df$m_z, '_', res_df$name))

  rownames(res_df) <- NULL
  return(res_df)
}

##' Overlay pixel coordinates on sample image
##'
##' Plots pixel coordinates as points overlaid on a sample image using ggplot2.
##'
##' @param coords Data frame with columns x, y, Sample.
##' @param sample_name Name of sample (matches folder).
##' @param image_file Path to image file.
##' @param point_size Numeric; size of overlay points.
##' @param alpha Numeric; transparency of overlay points.
##' @return ggplot object with overlay.
##' @keywords internal
preprocessing_plot_pixel_overlay_on_image <- function(coords, sample_name, image_file, point_size = 1, alpha = 0.7) {

  if (file.exists(image_file)){
    img <- jpeg::readJPEG(image_file)
    } else {
      stop(paste("Image file not found:", image_file))
    }

      static_coords = coords
      static_coords['x'] = (static_coords['x']-min(static_coords['x']))/2
      static_coords['y'] = (static_coords['y']-min(static_coords['y']))/2

  # Flip y axis for overlay
  p <- ggplot2::ggplot(static_coords, ggplot2::aes(x = .data$x, y = -.data$y)) +
    ggplot2::annotation_raster(img, 0, ncol(img), -nrow(img), 0) +
    ggplot2::geom_point(size = point_size, alpha = alpha, color = 'black') +
    ggplot2::theme_void() + ggplot2::coord_fixed()
return(p)
}

