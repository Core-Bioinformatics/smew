GCD <- function(x) {
  m = min(x)
  
  while (any(x %% m > 0)){
    m = m - 1
  }
  
  return(m)
}

select.label = function(i,current.metadata){
  current.col = current.metadata[i,'x']
  current.row = current.metadata[i,'y']
  selected.cluster = current.metadata |> dplyr::filter(x==current.col & y==current.row) |> dplyr::select(cluster)
  neighbour.clusters = current.metadata |> dplyr::filter(x %in% seq(current.col-1,current.col+1)) |> 
    dplyr::filter(y %in% seq(current.row-1,current.row+1)) |>
    dplyr::filter(x!=current.col | y!=current.row) |> 
    dplyr::select(cluster)
  if (length(neighbour.clusters)<5){
    return(as.character(selected.cluster))
  } else {
    top.cluster = sort(table(neighbour.clusters),decreasing = T)[1]
    if (top.cluster > 0.5*length(neighbour.clusters)){
      return(top.cluster)
    } else {
      return(selected.cluster)
    }
  }
}


# current.intensity.matrix <- t(intensity.matrix)
# current.metadata <- metadata
# current.metadata$Sample = current.metadata[,colnames(bulk.metadata)[1]]
# my_peak = anno[anno$m_z=='303.2324',]
# my_peak_expression <- t(current.intensity.matrix[my_peak$m_z,])
# quantiles = quantile(my_peak_expression, prob=c(25/100,(100-25)/100), type=1)
# current.metadata$cluster = factor(ifelse(my_peak_expression<=quantiles[1],'Low',
#                                          ifelse(my_peak_expression>=quantiles[2],'High','Medium')),levels=c('Low','Medium','High'))
# 
# current.intensity.matrix <- scale(x = current.intensity.matrix,center = T,scale = T)
# current.intensity.matrix <- current.intensity.matrix[complete.cases(current.intensity.matrix),]
# clusters <- kmeans(t(current.intensity.matrix),
#                    centers = 10,
#                    nstart = 1)
# current.metadata$cluster = clusters$cluster

# gcd.values = list()
# for (sample in unique(current.metadata$Group)){
#   metadata.sub = current.metadata |> filter(Group==sample)
#   pos <- metadata.sub[,c('x','y')]
#   rownames(pos)=metadata.sub$spot_id
#   rownames(metadata.sub)=metadata.sub$spot_id
#   # calculating spacing between points
#   min.x = min(pos$x)
#   min.y = min(pos$y)
#   pos$x = pos$x - min.x
#   pos$y = pos$y - min.y
#   gcd.x = GCD(round(unique(pos$x)[unique(pos$x)!=0]))
#   gcd.y = GCD(round(unique(pos$y)[unique(pos$y)!=0]))
#   if (gcd.x!=gcd.y){
#     message('Mismatching gap!')
#   } else {
#     gcd.values[[sample]]<-gcd.x
#   }
# }

smoothed.cluster <- function(metadata){
  rownames(metadata)=metadata$spot_id
  smoothed.clusters = c()
  for (sample in unique(metadata$Group)){
    print(sample)
    metadata.sub = metadata[metadata$Group==sample,]
    pos <- metadata.sub[,c('x','y','Group')]
    rownames(pos)=metadata.sub$spot_id
    rownames(metadata.sub)=metadata.sub$spot_id
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
    cl <- parallel::makeCluster(8)
    parallel::clusterExport(cl, c('metadata.sub','select.label'))
    new.clusters = pbapply::pbsapply(rownames(metadata.sub),FUN=function(x)select.label(x,metadata.sub),simplify = T,cl=cl)
    parallel::stopCluster(cl)
    names(new.clusters)=rownames(pos)
    smoothed.clusters = c(smoothed.clusters,new.clusters)
  }
  return(factor(smoothed.clusters,levels=levels(metadata$cluster)))
  # group each square of 9 points into 1
}

