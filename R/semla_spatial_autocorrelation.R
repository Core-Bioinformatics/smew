CorSpatialFeatures.default <- function (
    object,
    spatnet,
    across_all = FALSE,
    nCores = NULL,
    verbose = TRUE,
    ...
) {

  # Set global variables to NULL
  from <- to <- NULL

  # # Define multicore lapply function depending on OS
  # if (!.Platform$OS.type %in% c("windows", "unix")) {
  #   # Skip threading and use lapply instead
  #   if (verbose) cli_alert_danger("Threading not supported. Using single thread.")
    parLapplier <-  function(X, FUN, nCores) {
      res <- lapply(X, FUN)
      return(res)
    }
  #  }
  # else {
  #   parLapplier <- switch(.Platform$OS.type,
  #                         "windows" = .winLapply,
  #                         "unix" = .unixLapply)
  # }

  # get all spatial network spot barcode IDs
  all_spatial_network_spots <- Reduce(c, lapply(spatnet, function(x) x$from))

  # Check that objects match
  spots_in_spatnets <- !unique(all_spatial_network_spots) %in% rownames(object)
  if (any(spots_in_spatnets)) abort(glue("{sum(spots_in_spatnets)} spots in the spatial networks could not be found in the feature matrix.",
                                         "i" = "Make sure that the spatial networks share spot IDs with the feature matrix."))

  results <- lapply(seq_along(spatnet), function (i) {

    # pivot spatial network in long format to a wide format
    wide_spatial_network <- pivot_wider(spatnet[[i]] |> select(from, to) |> mutate(value = 1),
                                        names_from = "from", values_from = "value", values_fill = 0)

    # Convert wide spatial network to a matrix
    CN <- as.matrix(wide_spatial_network[, 2:ncol(wide_spatial_network)])
    rownames(CN) <- wide_spatial_network$to
    CN <- CN[colnames(CN), ]
    rm(wide_spatial_network)
    CN <- as(CN, "dgCMatrix")
    # Subset feature data to only include spots with neighbors
    x_subset <- object[colnames(CN), ]

    # Calculate lag matrixp
    lagMat <- (CN %*% x_subset) / Matrix::rowSums(CN)

    # return lagMat if the autocorrelation should be calculated across all samples
    if (across_all) {
      return(lagMat)
    }

    # Calculate spatial autocorrelation for each gene
    spatial_autocorrelation <- .colCors(x_subset, lagMat)

    # Summarize results
    results <- tibble(gene = names(spatial_autocorrelation), cor = spatial_autocorrelation) |>
      arrange(-cor)
  })

  # If across_all is set, calculate autocorrelations across all samples instead
  if (across_all) {
    lagMat <- do.call(rbind, results)

    # Calculate spatial autocorrelation for each gene
    spatial_autocorrelation <- .colCors(object[rownames(lagMat), ], lagMat)

    # Summarize results
    results <- tibble(gene = colnames(object), cor = spatial_autocorrelation) |>
      arrange(-cor)
  }

  return(results)
}

#' Calculate pairwise correlation across two matrices
#'
#' @param x,y Numeric matrices with identical dimensions
#'
#' @importFrom Matrix colMeans
#'
#' @return A numeric vector with correlation scores
#'
#' @noRd
.colCors = function(x, y) {
  x = sweep(x, 2, Matrix::colMeans(x))
  y = sweep(y, 2, Matrix::colMeans(y))
  cor = Matrix::colSums(x*y) / sqrt(Matrix::colSums(x*x)*Matrix::colSums(y*y))
  return(cor)
}

#' Find features with high spatial autocorrelation
#'
#' @param object An object (see details)
#' @param ... Arguments passed to other methods
#'
#' @family network-methods
#' @rdname cor-features
#'
#' @export
#'
#' @md
CorSpatialFeatures <- function(object, ...) {
  UseMethod(generic = 'CorSpatialFeatures', object = object)
}
