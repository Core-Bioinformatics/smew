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


#' @param features A character vector with features present in \code{Seurat} object. These
#' features need to be accessible with \code{\link{FetchData}}
#' @param assay_use Select assay to use for computation. If not specified, the default
#' assay will be used.
#' @param slot_use Select slot to use from assay object.
#'
#' @importFrom Seurat FetchData VariableFeatures GetAssayData
#' @importFrom rlang %||%
#'
#' @rdname cor-features
#' @family spatial-methods
#'
#' @author Ludvig Larsson
#'
#' @examples
#'
#' se_mbrain <- readRDS(system.file("extdata/mousebrain",
#'                                  "se_mbrain",
#'                                  package = "semla"))
#' se_mbrain <- se_mbrain |>
#'   ScaleData() |>
#'   RunPCA()
#'
#' # Compute spatial autocorrelation for variable features
#' spatgenes <- CorSpatialFeatures(se_mbrain,
#'                                 features = VariableFeatures(se_mbrain),
#'                                 nCores = 1)
#'
#' # Check genes with highest spatial autocorrelation
#' head(spatgenes[[1]])
#'
#' # Note that the top variable genes are blood related (hemoglobin genes)
#' # These genes have lower spatial autocorrelation since blood vessels
#' # typically only cover a few spots and more randomly dispersed throughput the tissue
#' head(VariableFeatures(se_mbrain))
#'
#' \donttest{
#' # The same principle can be used to estimate spatial autocorrelation for other features,
#' # for example dimensionality reduction vectors
#' spatpcs <- CorSpatialFeatures(se_mbrain,
#'                              features = paste0("PC_", 1:10),
#'                              nCores = 1)
#'
#' # Calculate spatial autocorrelation scores for principal components
#' head(spatpcs[[1]])
#'
#' # Compute spatial autocorrelation scores for multiple datasets
#' se_mcolon <- readRDS(system.file("extdata/mousecolon",
#'                                  "se_mcolon",
#'                                  package = "semla"))
#' se_merged <- MergeSTData(se_mbrain, se_mcolon) |>
#'   FindVariableFeatures()
#'
#' spatgenes <- CorSpatialFeatures(se_merged,
#'                                 features = VariableFeatures(se_merged),
#'                                 nCores = 1)
#'
#' # Check spatial autocorrelation scores mouse brain data
#' head(spatgenes[[1]])
#' # Check spatial autocorrelation scores mouse colon data
#' head(spatgenes[[2]])
#' }
#'
#' @export
#'
CorSpatialFeatures.Seurat <- function (
    object,
    features = NULL,
    assay_use = NULL,
    slot_use = "data",
    across_all = FALSE,
    nCores = NULL,
    verbose = TRUE,
    ...
) {

  # Validate Seurat object
  .check_seurat_object(object)

  # Get variable features if features=NULL
  features <- features %||% VariableFeatures(object)

  # Fetch features
  if (!is.null(assay_use)) {
    if (!requireNamespace("MatrixExtra", quietly = TRUE)) {
      abort("Package 'MatrixExtra' required.")
    }
    stopifnot(is.character(assay_use),
              length(assay_use) == 1)
    featureMat <- GetAssayData(object, assay = assay_use, slot = slot_use)
    featureMat <- MatrixExtra::t(featureMat[features, ])
  } else {
    featureMat <- FetchData(object, vars = features)
  }

  # Obtain spatial networks
  spatnet <- GetSpatialNetwork(object)

  # Compute spatial autocorrelation
  spatfeatures <- CorSpatialFeatures(featureMat, spatnet, across_all, nCores, verbose, ...)

  return(spatfeatures)
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
