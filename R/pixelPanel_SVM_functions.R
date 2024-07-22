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
  nNeighbors <- nNeighbors %||% 8
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
    verbose = TRUE,
    ...
) {
  # Set global variables to NULL
  from <- to <- NULL

  # # Define multicore lapply function depending on OS
  # if (!.Platform$OS.type %in% c("windows", "unix")) {
  #   # Skip threading and use lapply instead
  #   if (verbose) cli_alert_danger("Threading not supported. Using single thread.")
  parLapplier <-  function(X, FUN) {
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
  my.fn = function (i) {
    # pivot spatial network in long format to a wide format
    step.by.step = spatnet[[i]] |> dplyr::select(from, to) |> mutate(value = 1)
    customers <- unique(step.by.step$from)
    products <- unique(step.by.step$to)
    step.by.step$row <- match(step.by.step$from, customers)
    step.by.step$col <- match(step.by.step$to, products)

    df_sparse <- sparseMatrix(
      i = step.by.step$row,
      j = step.by.step$col,
      x = step.by.step$value,
      dimnames = list(customers,
                      products)
    )
    CN = df_sparse
    # Convert wide spatial network to a matrix
    CN <- CN[colnames(CN), ]
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
    return(results)
  }
  # results = list()
  # for (i in seq_along(spatnet)){
  #   results[[i]]=my.fn(i)
  # }
  results <- lapply(seq_along(spatnet), my.fn)
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
  print('making network')
  spatnet = GetSpatialNetwork(coords)
  print('calculating spatial autocorrelation across network')
  spatgenes = CorSpatialFeatures(intensity.matrix,spatnet)
  names(spatgenes)=unique(metadata[,colnames(bulk.metadata)[1]])
  return(spatgenes)
}


GCD <- function(x) {
  m = min(x)

  while (any(x %% m > 0)){
    m = m - 1
  }

  return(m)
}


make.submatrix <- function(metadata.sub,intensity.sub,name){
  pos <- metadata.sub[,c('x','y')]
  rownames(pos)=metadata.sub$spot_id

  # calculating spacing between points
  min.x = min(pos$x)
  min.y = min(pos$y)
  pos$x = pos$x - min.x
  pos$y = pos$y - min.y
  gcd.x = GCD(round(unique(pos$x)[unique(pos$x)!=0]))
  gcd.y = GCD(round(unique(pos$y)[unique(pos$y)!=0]))
  pos$x = round((pos$x/gcd.x)+1)
  pos$y = round((pos$y/gcd.y)+1)

  # group each square of 9 points into 1
  not.started = T
  for (current.row in seq.int(ceiling(max(pos$y/3)))){
    for (current.col in seq.int(ceiling(max(pos$x/3)))){

      included.points = pos[(pos$x %in% seq(3*current.col,3*current.col+2))&
                              (pos$y %in% seq(3*current.row,3*current.row+2)),]
      included.intensity = intensity.sub[rownames(included.points),,drop=F]

      if (nrow(included.intensity)!=0){
        if (nrow(included.intensity)>1){
          # take the mean of all points within the square that have measurements
          squashed.intensity = colMeans(included.intensity)
        } else {
          # otherwise just take the one point
          squashed.intensity = included.intensity
        }
        if (not.started){
          overall.intensity = squashed.intensity
          # use the indices as the new x and y coordinates (these will not match the original points)
          overall.meta = data.frame('x'=current.col,'y'=current.row,'num_points'=nrow(included.intensity))
          not.started=F
        } else {

          overall.intensity = rbind(overall.intensity,squashed.intensity)
          overall.meta = rbind(overall.meta,c('x'=current.col,'y'=current.row,'num_points'=nrow(included.intensity)))
        }
      }
    }
  }
  rownames(overall.meta)=paste0(name,'_',1:nrow(overall.intensity))
  rownames(overall.intensity)=paste0(name,'_',1:nrow(overall.intensity))
  return(list('intensity'=overall.intensity,'meta'=overall.meta))
}

cross.cor.persample <- function(feature1,feature2,sample,combined.intensity,combined.metadata,list.of.adjacency){
  weight = list.of.adjacency[[sample]]
  x = combined.intensity[combined.metadata$sample==sample,feature1]
  y = combined.intensity[combined.metadata$sample==sample,feature2]
  cv1 <- x %o% y
  cv2 <- y %o% x
  cv1[is.na(cv1)] <- 0
  cv2[is.na(cv2)] <- 0
  return(sum(weight * ( cv1 + cv2 ), na.rm=TRUE))
}
apply.each.feature <- function(i,N,W,all.combos,combined.intensity,combined.metadata){
  feature1 = all.combos[i,'Var1']
  feature2 = all.combos[i,'Var2']
  v <- sqrt(sum(combined.intensity[,feature1]^2, na.rm=TRUE) * sum(combined.intensity[,feature2]^2, na.rm=TRUE))
  cv = sum(sapply(unique(combined.metadata$sample),function(x)cross.cor.persample(feature1,feature2,x,combined.intensity,combined.metadata,list.of.adjacency)))

  SCC <- (N/W) * (cv/v) / 2.0
  return(SCC)
}

apply.submatrix <- function(metadata,intensity.matrix){
  for (sample in unique(metadata$Group)){
    # do the submatrix bit
    metadata.sub = metadata[metadata$Group==sample,]
    intensity.sub = intensity.matrix[metadata$Group==sample,]
    rownames(intensity.sub)=metadata.sub$spot_id
    my.smoothed.output = make.submatrix(metadata.sub,intensity.sub,sample)
    meta = my.smoothed.output$meta
    meta$sample = sample
    message(paste0("Created reduced matrix for ",sample))
    if (sample==unique(metadata$Group)[1]){
      combined.intensity = my.smoothed.output$intensity
      combined.meta = meta
    } else {
      combined.intensity = rbind(combined.intensity,my.smoothed.output$intensity)
      combined.meta = rbind(combined.meta,meta)
    }
  }
  combined.intensity = scale(combined.intensity)
  return(list('intensity'=combined.intensity,'metadata'=combined.meta))
}

calculate.cross.cor <- function(metadata,intensity.matrix,num.genes,ncores,svm_identification){
  message("Creating reduced matrix:")
  submatrix = apply.submatrix(metadata,intensity.matrix)
  combined.meta = submatrix$metadata
  combined.intensity = submatrix$intensity

  # create weight adjacency matrices per sample
  list.of.adjacency <- list()
  W = 0
  message("Creating weight adjacency matrices:")
  for (sample in unique(combined.meta$sample)){
    #for (sample in 'bleo_d21_4a'){
    pos <- combined.meta[combined.meta$sample==sample,c('x','y')]
    weight <- MERINGUE::getSpatialNeighbors(pos,filterDist = sqrt(5))
    rs <- rowSums(weight)
    rs[rs == 0] <- 1
    weight <- weight/rs
    weight = weight[rownames(pos),rownames(pos)]
    list.of.adjacency[[sample]]<-weight
    W = W + sum(weight)
    message(paste0("Created weight adjacency matrix for ",sample))
  }
  N = nrow(combined.meta)

  svm.peaks = Reduce(union,lapply(FUN = function(x)head(x$gene,num.genes),X=svm_identification))
  message(paste0("Calculating spatial cross-correlation for ",length(svm.peaks)," peaks"))
  all.combos<-expand.grid(svm.peaks, svm.peaks,stringsAsFactors = F)
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, c('all.combos','combined.intensity','combined.meta','N','W','apply.each.feature','list.of.adjacency','cross.cor.persample'),envir = environment())
  all.combos$scc = pbapply::pbsapply(rownames(all.combos),FUN = function(x)apply.each.feature(x,N,W,all.combos,combined.intensity,combined.meta),simplify = T,cl = cl)
  parallel::stopCluster(cl)
  all.combos = as.data.frame(tidyr::pivot_wider(all.combos,names_from='Var2',values_from = 'scc'))
  rownames(all.combos)=all.combos$Var1
  all.combos = all.combos[,2:ncol(all.combos)]
  return(all.combos)
}

