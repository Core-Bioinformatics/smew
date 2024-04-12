denoise.intensity <- function(metadata,intensity.matrix,sample){
  current.metadata = metadata |> filter(Group==sample)
  current.intensity = intensity.matrix[metadata$Group==sample,]
  rownames(current.intensity)=current.metadata$spot_id
  rownames(current.metadata)=current.metadata$spot_id
  current.intensity = log2(current.intensity+1)
  current.intensity = current.intensity[,!apply(current.intensity, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
  gcd = gcd.values[[sample]]
  pca = prcomp(current.intensity,center = T,scale=T)
  cor.matrix = cor(t(pca$x[,1:50]))
  distance.matrix = as.matrix(dist(current.metadata[,c('x','y')]))
  w = cor.matrix * distance.matrix
  all.intensities = t(pbapply::pbsapply(rownames(current.intensity),simplify = T,USE.NAMES = T,FUN = function(x)pick.neighbours.per.pixel(x,cor.matrix,distance.matrix,current.intensity,gcd,w)))
}

pick.neighbours.per.pixel <- function(i,cor.matrix,distance.matrix,current.intensity,gcd,w){
  distance.neighbours = distance.matrix[i, ]
  distance.neighbours = names(distance.neighbours[distance.neighbours < 3*gcd])
  distance.neighbours = distance.neighbours[distance.neighbours != i]
  cor.neighbours = head(names(sort(cor.matrix[i,distance.neighbours],decreasing = T)),3)
  new.intensity = 0.5*current.intensity[i,]
  weight.sum = sum(w[i,cor.neighbours])
  for (neighbour in cor.neighbours){
    new.intensity = new.intensity + 0.5*(w[i,neighbour]*current.intensity[neighbour,]/weight.sum)
  }
  return(new.intensity)
}
