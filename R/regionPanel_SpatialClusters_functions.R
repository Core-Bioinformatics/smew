select.label = function(i,current.metadata){
  current.col = current.metadata[i,'x']
  current.row = current.metadata[i,'y']
  selected.cluster = current.metadata |> dplyr::filter(x==current.col & y==current.row) |> dplyr::select(cluster)
  neighbour.clusters = current.metadata |> dplyr::filter(x %in% seq(current.col-1,current.col+1)) |>
    dplyr::filter(y %in% seq(current.row-1,current.row+1)) |>
    dplyr::filter(x!=current.col | y!=current.row) |>
    dplyr::select(cluster)
  if (nrow(neighbour.clusters)<5){
    return(selected.cluster$cluster)
  } else {
    top.cluster = sort(table(neighbour.clusters$cluster),decreasing = T)[1]
    if (top.cluster > 0.5*nrow(neighbour.clusters)){
      return(names(top.cluster))
    } else {
      return(selected.cluster$cluster)
    }
  }
}

smoothed.cluster <- function(metadata.sub,gcd.values){
  # this should only be 1 sample
  # if (!is.factor(metadata$cluster)){
  #   metadata$cluster = factor(metadata$cluster)
  # }
  sample = unique(metadata.sub$Group)
  metadata.sub$cluster = as.character(metadata.sub$cluster)
  rownames(metadata.sub)=metadata.sub$spot_id
#  smoothed.clusters = c()
#  for (sample in unique(metadata$Group)){
#    metadata.sub = metadata[metadata$Group==sample,]
    pos <- metadata.sub[,c('x','y','Group')]
    # calculating spacing between points
    min.x = min(pos$x)
    min.y = min(pos$y)
    pos$x = pos$x - min.x
    pos$y = pos$y - min.y
    gcd = gcd.values[[sample]]
    pos$x = round((pos$x/gcd)+1)
    pos$y = round((pos$y/gcd)+1)
    metadata.sub$x = pos$x
    metadata.sub$y = pos$y
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    parallel::clusterExport(cl, c('metadata.sub','select.label'),envir=environment())
    smoothed.clusters = pbapply::pbsapply(rownames(metadata.sub),FUN=function(x)select.label(x,metadata.sub),simplify = T,cl=cl)
    parallel::stopCluster(cl)
    names(smoothed.clusters)=rownames(pos)
 #   smoothed.clusters = c(smoothed.clusters,new.clusters)
 # }
#  return(factor(smoothed.clusters,levels=levels(metadata$cluster)))
  return(smoothed.clusters)
  # group each square of 9 points into 1
}

