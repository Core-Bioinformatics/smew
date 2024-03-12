# these functions are from the semla package

GetSpatialNetwork <- function (
    object,
    nNeighbors = 6,
    maxDist = NULL,
    minK = 0,
    ...
) {
  
  # Set global variables to NULL
  barcode <- x <- y <- sampleID <- from <- distance <- NULL
  
  # Check object object class
  if (!any(class(object) %in% c("data.frame", "matrix", "tbl")))
    abort(glue::glue("Invalid class '{class(object)}'."))
  if (ncol(object) != 4)
    abort(glue::glue("Invalid number of columns '{ncol(object)}'. Expected 4."))
  # match column names
  col_checks <- colnames(object) %in% c("barcode", "x", "y", "sampleID")
  if (any(!col_checks)) abort(glue::glue("Invalid column name(s) '{paste(colnames(object)[!col_checks], collapse = ', ')}'. ",
                                   "Expected 'barcode', 'x', 'y' and 'sampleID'"))
  if (!all(
    object |> summarize(
      check_barcode = is.character(barcode),
      check_x = is.numeric(x),
      check_y = is.numeric(y),
      check_sample = is.numeric(sampleID)
    ) |>
    unlist()
  )) {
    abort(glue::glue("Invalid column class(es)."))
  }
  
  # Set number of nearest neighbors if NULL
  nNeighbors <- nNeighbors %||% 6
  # Split coordinates by sample
  xys.list <- object |>
    group_by(sampleID) |>
    group_split() |>
    setNames(paste0(unique(object$sampleID)))
  
  # Compute network
  knn_long.list <- setNames(lapply(names(xys.list), function(sampleID) {
    
    xys_subset <- xys.list[[sampleID]]
    spotnames <- setNames(xys_subset$barcode, nm = c(1:nrow(xys_subset)) |> paste0())
    knn_spatial <- dbscan::kNN(x = xys_subset[, c("x", "y")] |> as.matrix(), k = nNeighbors)
    maxDist <- maxDist %||% (apply(knn_spatial$dist, 1, min) |> min())*1.2
    knn_long <- tibble::tibble(from = rep(1:nrow(knn_spatial$id), nNeighbors),
                               to = as.vector(knn_spatial$id),
                               distance = as.vector(knn_spatial$dist))
    
    knn_long$from <- spotnames[knn_long$from]
    knn_long$to <- spotnames[knn_long$to]
    knn_long <- knn_long |>
      add_count(from) |>
      filter(distance <= maxDist, n > minK) |>
      add_count(from, name = "nn") |>
      dplyr::select(-n)
    
    # Merge with coordinates
    knn_long <-
      left_join(
        x = knn_long,
        y = setNames(xys_subset[, 1:3], nm = c("barcode", "x_start", "y_start")),
        by = c("from" = "barcode")
      )
    knn_long <-
      left_join(
        x = knn_long,
        y = setNames(xys_subset[, 1:3], nm = c("barcode", "x_end", "y_end")),
        by = c("to" = "barcode")
      )
    
    return(knn_long)
  }), nm = names(xys.list))
  
  return(knn_long.list)
}

CorSpatialFeatures <- function (
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
  if (any(spots_in_spatnets)) abort(glue::glue("{sum(spots_in_spatnets)} spots in the spatial networks could not be found in the feature matrix.",
                                         "i" = "Make sure that the spatial networks share spot IDs with the feature matrix."))
  
  results <- lapply(seq_along(spatnet), function (i) {
    
    # pivot spatial network in long format to a wide format
    wide_spatial_network <- pivot_wider(spatnet[[i]] |> dplyr::select(from, to) |> mutate(value = 1),
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

.colCors = function(x, y) {
  x = sweep(x, 2, Matrix::colMeans(x))
  y = sweep(y, 2, Matrix::colMeans(y))
  cor = Matrix::colSums(x*y) / sqrt(Matrix::colSums(x*x)*Matrix::colSums(y*y))
  return(cor)
}

run_svm <- function(intensity.matrix,metadata,bulk.metadata){
  coords = metadata[,c('spot_id','x','y',colnames(bulk.metadata)[1])]
  rownames(intensity.matrix)=metadata$spot_id
  colnames(coords)=c('barcode','x','y','sampleID')
  coords$barcode = as.character(coords$barcode)
  coords = as.tibble(coords)
  coords$sampleID = as.numeric(as.factor(coords$sampleID))
  spatnet = GetSpatialNetwork(coords)
  spatgenes = CorSpatialFeatures(intensity.matrix,spatnet,nCores=1)
  names(spatgenes)=unique(metadata[,colnames(bulk.metadata)[1]])
  return(spatgenes)
}
