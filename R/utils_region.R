##' Create pseudobulk region-level expression matrix
##'
##' Aggregates intensity values by user-defined regions and samples, producing a pseudobulk expression matrix and corresponding metadata for region-level differential analysis.
##'
##' @param intensity.matrix Data frame of intensity values (features x samples)
##' @param metadata Data frame of sample metadata (must include 'Sample' column)
##' @param bulk.metadata Data frame of bulk sample metadata
##' @param sample.id.metadata.column Name of the sample ID column in metadata
##' @param region.ids Vector of region IDs for each sample (length must match nrow(metadata))
##' @param minimum.pixels Minimum number of pixels per region per sample to be considered (default 3)
##' @return A list with:
##'   \item{intensity.matrix}{Matrix of aggregated intensities (features x pseudobulked samples)}
##'   \item{metadata}{Data frame of metadata for pseudobulked samples}
##' @details
##' For each unique region within each sample, the function computes the mean intensity across all pixels in that region, provided the region contains at least \code{minimum.pixels} pixels. The resulting matrix and metadata are suitable for downstream region-level differential analysis.
##' @keywords internal
region_utils_create_bulk_exp_regions <- function(intensity.matrix,
                                    metadata,
                                    bulk.metadata,
                                    sample.id.metadata.column,
                                    region.ids,
                                    minimum.pixels = 3
) {
  # Add extra checks
  intensity.matrix$sample <- metadata$Sample
  intensity.matrix$region <- factor(region.ids)

  intensity.matrix.mean <- dplyr::filter(intensity.matrix, !is.na(region)) |>
    dplyr::group_by(sample, region) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), mean), n = dplyr::n(), .groups = 'drop')

  intensity.matrix.mean = intensity.matrix.mean[intensity.matrix.mean$n > minimum.pixels, ]
  sample.names = paste0(intensity.matrix.mean$sample, '_', intensity.matrix.mean$region)
  # Remove 'n' column before transpose so it never becomes a row
  intensity.matrix.mean <- intensity.matrix.mean[, !(colnames(intensity.matrix.mean) %in% c('sample', 'region', 'n'))]
  intensity.matrix.mean <- dplyr::rename_with(
    as.data.frame(
      t(
        as.matrix(intensity.matrix.mean)
      )
    ),
    ~sample.names
  )

  for (region in unique(region.ids[!is.na(region.ids)])) {
    current.metadata = bulk.metadata
    current.metadata$AllSamples = 'AllSamples'
    current.metadata = data.frame(lapply(current.metadata, function(x) paste(x, region, sep = "_")))
    if (region == unique(region.ids[!is.na(region.ids)])[1]) {
      metadata.mean = current.metadata
    } else {
      metadata.mean = rbind(metadata.mean,current.metadata)
    }
  }
  metadata.mean = metadata.mean[metadata.mean[,1]%in%colnames(intensity.matrix.mean),]
  intensity.matrix.mean = intensity.matrix.mean[,metadata.mean[,1]]
  return(list('intensity.matrix'=as.matrix(intensity.matrix.mean),'metadata'=metadata.mean))
}

##' Colour blender
##'
##' Blend numeric feature/dimension values into RGB colour codes for visualisation.
##'
##' Given a matrix or data frame of 2 or 3 numeric features/dimensions, and a vector specifying which colour channel (red, green, blue) to assign to each, this function returns RGB colour codes for each row. Used for coloured plotting in region analysis.
##'
##' @param data Matrix or data frame of numeric values (n x 2 or n x 3).
##' @param channels.use Character vector of length 2 or 3 specifying which colour channel to assign to each column of data. Must be a permutation of 'red', 'green', 'blue'.
##' @return Character vector of RGB colour codes (length n).
##' @details
##' If data has 2 columns, the third channel is set to zero. Checks for valid and unique channel assignments. Used for coloured plotting in region analysis.
##' @keywords internal
region_utils_colour_blender <- function (
  data,
  channels.use = NULL
) {
  rgb.order <- stats::setNames(1:3, c("red", "green", "blue"))
  if (!length(channels.use) == ncol(data)) {
    stop(paste0("channels.use must be same length as number of features or dimensions"))
  } else if (!all(channels.use %in% names(rgb.order))) {
    stop("Invalid colour names in channels.use. Valid options are: 'red', 'green' and 'blue'")
  } else if (sum(duplicated(channels.use))){
    stop("Duplicate colour names are not allowed in channels.use")
  }
  col.order <- rgb.order[channels.use]

  if (ncol(data) == 2) {
    first_vec <- data[, 1]
    second_vec <- data[, 2]
    data <- matrix(data = 0, nrow = nrow(data), ncol = 3)
    data[, col.order[1]] <- first_vec; data[, col.order[2]] <- second_vec
  } else if (ncol(data) == 3) {
    data <- data[, col.order]
  }
  return(grDevices::rgb(data))
}


#' Smooth cluster labels across spatial neighborhoods
#'
#' For each pixel, replaces its cluster label with the most common label among its 8 spatial neighbors if a majority exists, otherwise retains the original label. Used to reduce noise in spatial cluster assignments.
#'
#' @param metadata.sub Data frame of sample metadata for a single sample (must include columns 'Sample', 'pixel_id', 'x_tf', 'y_tf', and the cluster column).
#' @param cluster.col Name of the column containing cluster labels to smooth.
#' @return Named character vector of smoothed cluster labels (names are pixel IDs).
#' @keywords internal
#' @details
#' Uses parallel processing to speed up smoothing across all pixels. Calls region_utils_select_majority_label for each pixel.
region_utils_smooth_cluster <- function(metadata.sub,cluster.col){
  sample = base::unique(metadata.sub$Sample)
  metadata.sub[,cluster.col] = base::as.character(metadata.sub[,cluster.col])
  base::rownames(metadata.sub)=metadata.sub$pixel_id
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  parallel::clusterExport(cl, c('metadata.sub','region_utils_select_majority_label'),envir=base::environment())
    smoothed.clusters = pbapply::pbsapply(base::rownames(metadata.sub),FUN=function(x)region_utils_select_majority_label(x,metadata.sub,cluster.col),simplify = TRUE, cl=cl)
  parallel::stopCluster(cl)
  base::names(smoothed.clusters)=metadata.sub$pixel_id
  return(smoothed.clusters)
}


#' Select majority label from spatial neighborhood
#'
#' Select the smoothed cluster label for a given pixel based on its spatial neighborhood.
#'
#' @param i Pixel ID (row name in metadata).
#' @param current.metadata Data frame of sample metadata (must include 'x_tf', 'y_tf', and the cluster column).
#' @param cluster.col Name of the column containing cluster labels.
#' @keywords internal
#' @return Smoothed cluster label for the pixel (character).
#' @details
#' Returns the most common cluster label among the 8 neighbors if it constitutes a majority (more than half), otherwise returns the original label.
region_utils_select_majority_label = function(i, current.metadata, cluster.col){
  current.col = current.metadata[i,'x_tf']
  current.row = current.metadata[i,'y_tf']
    selected.cluster = dplyr::select(dplyr::filter(current.metadata, .data$x_tf == current.col & .data$y_tf == current.row), dplyr::one_of(cluster.col))
  neighbour.clusters = dplyr::select(
    dplyr::filter(
      dplyr::filter(
        dplyr::filter(current.metadata, .data$x_tf %in% base::seq(current.col-1,current.col+1)),
        .data$y_tf %in% base::seq(current.row-1,current.row+1)),
      .data$x_tf!=current.col | .data$y_tf!=current.row),
    dplyr::one_of(cluster.col))
  if (base::nrow(neighbour.clusters)<5){
    return(selected.cluster[,cluster.col])
  } else {
    top.cluster = base::sort(base::table(neighbour.clusters[,cluster.col]),decreasing = TRUE)[1]
    if (top.cluster > (0.5 * base::nrow(neighbour.clusters))){
      return(base::names(top.cluster))
    } else {
      return(selected.cluster[,cluster.col])
    }
  }
}