##' Grid-based spatial autocorrelation pipeline
##'
##' Computes spatial autocorrelation for features using grid-based neighbours (gcd 1). No external neighbour-finding packages required; based on the semla spatial lag concept.
##'

##' Find grid neighbours for each spot
##'
##' @param metadata Data frame with columns: pixel_id, x_tf, y_tf, Sample.
##' @return List of neighbour indices per spot (named by pixel_id).
##' @keywords spatial, neighbours
find_grid_neighbours <- function(metadata) {
  spots <- metadata$pixel_id
  coords <- metadata[, c('x_tf', 'y_tf')]
  neighbour_list <- lapply(seq_len(nrow(coords)), function(i) {
    xi <- coords$x_tf[i]
    yi <- coords$y_tf[i]
    neighbours <- which(
      abs(coords$x_tf - xi) <= 1 &
      abs(coords$y_tf - yi) <= 1 &
      !(coords$x_tf == xi & coords$y_tf == yi)
    )
    spots[neighbours]
  })
  names(neighbour_list) <- spots
  neighbour_list
}

##' Compute spatial lag for each spot
##'
##' @param intensity Matrix (spots x features), rownames = spot_id.
##' @param neighbour_list List of neighbour spot_ids per spot.
##' @return Matrix of spatial lags (spots x features).
##' @keywords spatial, lag
compute_spatial_lag <- function(intensity, neighbour_list) {
  lag_mat <- matrix(NA, nrow = nrow(intensity), ncol = ncol(intensity),
                   dimnames = dimnames(intensity))
  for (i in seq_len(nrow(intensity))) {
    neighbours <- neighbour_list[[rownames(intensity)[i]]]
    if (length(neighbours) > 0) {
      lag_mat[i, ] <- colMeans(intensity[neighbours, , drop = FALSE], na.rm = TRUE)
    } else {
      lag_mat[i, ] <- NA
    }
  }
  lag_mat
}

##' Compute spatial autocorrelation (correlation between feature and spatial lag)
##'
##' @param intensity Matrix (spots x features).
##' @param lag_mat Matrix (spots x features).
##' @return Named vector of autocorrelation per feature.
##' @keywords spatial, autocorrelation
compute_autocorrelation <- function(intensity, lag_mat) {
  n <- ncol(intensity)
  cors <- numeric(n)
  for (i in seq_len(n)) {
    cors[i] <- stats::cor(intensity[, i], lag_mat[, i], use = 'everything', method = 'pearson')
  }
  names(cors) <- colnames(intensity)
  cors
}

##' Main pipeline for spatial autocorrelation
##'
##' @param intensity Matrix (spots x features), rownames = pixel_id.
##' @param metadata Data frame with pixel_id, x_tf, y_tf, Sample.
##' @return Named list of autocorrelation vectors per Sample.
##' @keywords spatial, autocorrelation, pipeline
spatial_autocorrelation_pipeline <- function(intensity, metadata) {
  results <- pbapply::pblapply(unique(metadata$Sample), function(s) {
    idx <- which(metadata$Sample == s)
    sub_meta <- metadata[idx, ]
    sub_int <- intensity[sub_meta$pixel_id, , drop = FALSE]
    neighbour_list <- find_grid_neighbours(sub_meta)
    lag_mat <- compute_spatial_lag(sub_int, neighbour_list)
    cors <- compute_autocorrelation(sub_int, lag_mat)
    df <- data.frame(peak = names(cors), cor = cors, stringsAsFactors = FALSE)
    df <- df[order(-df$cor), ]
    rownames(df) <- NULL
    df
  })
  names(results) <- unique(metadata$Sample)
  results
}

##' Calculate grid-based spatial cross-correlation for a sample
##'
##' @param feature1 Name of the first feature (column) to correlate
##' @param feature2 Name of the second feature (column) to correlate
##' @param sample Name or ID of the sample to process
##' @param intensity.matrix Matrix of intensities (rows = spots/pixels, columns = features)
##' @param metadata Data frame with pixel_id, x_tf, y_tf, Sample, etc.
##' @param list.of.adjacency List of adjacency matrices (one per sample), each a sparse matrix
##' @return Numeric value: spatial cross-correlation for the two features in the given sample
cross.cor.persample <- function(feature1, feature2, sample, intensity.matrix, metadata, list.of.adjacency) {
  weight = list.of.adjacency[[sample]]
  idx = which(metadata$Sample == sample)
  x = intensity.matrix[idx, feature1]
  y = intensity.matrix[idx, feature2]
  names(x) <- rownames(intensity.matrix)[idx]
  names(y) <- rownames(intensity.matrix)[idx]
  # Demean x and y within sample
  dx = x - mean(x, na.rm=TRUE)
  dy = y - mean(y, na.rm=TRUE)
  # Efficiently sum only at nonzero entries of weight (sparse)
  # weight is dgCMatrix: use @i (row), @p (col ptr), @x (values)
  # For each column j, entries from weight@i[(weight@p[j]+1):weight@p[j+1]] are nonzero rows
  cv = 0
  for (j in seq_len(ncol(weight))) {
    col_start = weight@p[j] + 1
    col_end = weight@p[j+1]
    if (col_end >= col_start) {
      rows = weight@i[col_start:col_end] + 1
      wvals = weight@x[col_start:col_end]
      # dx[rows], dy[j], dy[rows], dx[j]
      cv = cv + sum(wvals * (dx[rows] * dy[j] + dy[rows] * dx[j]))
    }
  }
  return(cv)
}

apply.each.feature <- function(i, N, W, all.combos, intensity.matrix, metadata, list.of.adjacency){
  feature1 = all.combos[i,'Var1']
  feature2 = all.combos[i,'Var2']
  cv_sum = 0
  v_sum = 0
  for (sample in unique(metadata$Sample)){
    idx = which(metadata$Sample == sample)
    x = intensity.matrix[idx, feature1]
    y = intensity.matrix[idx, feature2]
    dx = x - mean(x, na.rm=TRUE)
    dy = y - mean(y, na.rm=TRUE)
    weight = list.of.adjacency[[sample]]
    cv = 0
    for (j in seq_len(ncol(weight))) {
      col_start = weight@p[j] + 1
      col_end = weight@p[j+1]
      if (col_end >= col_start) {
        rows = weight@i[col_start:col_end] + 1
        wvals = weight@x[col_start:col_end]
        cv = cv + sum(wvals * (dx[rows] * dy[j] + dy[rows] * dx[j]))
      }
    }
    v = sqrt(sum(dx^2, na.rm=TRUE) * sum(dy^2, na.rm=TRUE))
    cv_sum = cv_sum + cv
    v_sum = v_sum + v
  }
  SCC <- (N/W) * (cv_sum/v_sum) / 2.0
  return(SCC)
}

spatial_cross_cor <- function(top.peaks, intensity.matrix, metadata, ncores = 1) {
    list.of.adjacency <- list()
    W = 0
    message("Creating grid-based adjacency matrices:")
    for (sample in unique(metadata$Sample)) {
    pos <- metadata[metadata$Sample == sample, c('x_tf', 'y_tf')]
    spot_ids <- rownames(pos)
    n_spots <- nrow(pos)
    # Prepare i, j, x for sparseMatrix
    i_idx <- integer()
    j_idx <- integer()
    x_val <- numeric()
    for (i in seq_len(n_spots)) {
        xi <- pos$x_tf[i]
        yi <- pos$y_tf[i]
        neighbours <- which(
        abs(pos$x_tf - xi) <= 1 &
            abs(pos$y_tf - yi) <= 1 &
            !(pos$x_tf == xi & pos$y_tf == yi)
        )
        if (length(neighbours) > 0) {
        i_idx <- c(i_idx, rep(i, length(neighbours)))
        j_idx <- c(j_idx, neighbours)
        x_val <- c(x_val, rep(1, length(neighbours)))
        }
    }
    if (length(i_idx) > 0) {
        adj_mat <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val, dims = c(n_spots, n_spots), dimnames = list(spot_ids, spot_ids))
        rs <- Matrix::rowSums(adj_mat)
        rs[rs == 0] <- 1
        adj_mat <- adj_mat / rs
    } else {
        adj_mat <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n_spots, n_spots), dimnames = list(spot_ids, spot_ids))
    }
    list.of.adjacency[[sample]] <- adj_mat
    W = W + sum(adj_mat)
    message(paste0("Created grid-based adjacency matrix for ", sample))
    }
    N = nrow(metadata)

    # Get union of top spatially variable peaks across samples
    message(paste0("Calculating spatial cross-correlation for ",length(top.peaks)," peaks"))
    # Build all pairwise combinations
    all.combos <- base::expand.grid(top.peaks, top.peaks, stringsAsFactors = F)

    if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    message(paste0("Calculating spatial cross-correlation in parallel using ", ncores, " cores..."))
    parallel::clusterExport(cl, c('all.combos','intensity.matrix','metadata','N','W','apply.each.feature','list.of.adjacency','cross.cor.persample'),envir = environment())
    setup <- parallel::clusterEvalQ(cl, {library(Matrix)})
    all.combos$scc = pbapply::pbsapply(seq_len(nrow(all.combos)),FUN = function(x) apply.each.feature(x, N, W, all.combos, intensity.matrix, metadata, list.of.adjacency), simplify = TRUE, cl = cl)
    parallel::stopCluster(cl)
    } else {
    message("Calculating spatial cross-correlation sequentially...")
    all.combos$scc = pbapply::pbsapply(seq_len(nrow(all.combos)),FUN = function(x) apply.each.feature(x, N, W, all.combos, intensity.matrix, metadata, list.of.adjacency), simplify = TRUE)
    }

    # Reshape to wide matrix
    all.combos = as.data.frame(tidyr::pivot_wider(all.combos, names_from='Var2', values_from = 'scc'))
    rownames(all.combos) <- all.combos$Var1
    all.combos <- all.combos[,2:ncol(all.combos)]
    return(all.combos)
}