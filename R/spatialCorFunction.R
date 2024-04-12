# rda.files <- list.files(pattern = '\\.rda$')
# for(fl in rda.files) load(fl)

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

calculate.cross.cor <- function(metadata,intensity.matrix,num.genes,ncores){
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
    weight <- getSpatialNeighbors(pos,filterDist = sqrt(5))
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
  parallel::clusterExport(cl, c('all.combos','combined.intensity','combined.meta','N','W','apply.each.feature','list.of.adjacency','cross.cor.persample'))
  all.combos$scc = pbapply::pbsapply(rownames(all.combos),FUN = function(x)apply.each.feature(x,N,W,all.combos,combined.intensity,combined.meta),simplify = T,cl = cl)
  parallel::stopCluster(cl)
  all.combos = as.data.frame(tidyr::pivot_wider(all.combos,names_from='Var2',values_from = 'scc'))
  rownames(all.combos)=all.combos$Var1
  all.combos = all.combos[,2:ncol(all.combos)]
  return(all.combos)
}

# res = calculate.cross.cor(metadata,intensity.matrix,10,8)
# spatial.cross.cor = res
# save('spatial.cross.cor',file = 'spatial_cross_cor.rda')
# #ideas for tab in app
# # ggcorrplot-style heatmap with cluster annotations into modules
# # show the average cluster expression in each sample
# # show bubble plot with top peaks per sample showing overlap in spatial autocorrelation in each
# clust = hclust(as.dist(1-all.combos[,1:nrow(all.combos)]))
# corrplot::corrplot(as.matrix(all.combos[clust$order,clust$order]))
# # result <- microbenchmark::microbenchmark(
# #   pbapply::pbsapply(rownames(all.combos),FUN = function(x)apply.each.feature(x,N,W,all.combos,combined.intensity,combined.meta),simplify = T,cl = cl),
# #   pbapply::pbsapply(rownames(all.combos),FUN = function(x)apply.each.feature.loop(x,N,W,all.combos,combined.intensity,combined.meta),simplify = T,cl = cl),
# #   times = 10
# # )
# # Unit: seconds
# # expr
# # pbapply::pbsapply(rownames(all.combos), FUN = function(x) apply.each.feature(x,      N, W, all.combos, combined.intensity, combined.meta), simplify = T,      cl = cl)
# # pbapply::pbsapply(rownames(all.combos), FUN = function(x) apply.each.feature.loop(x,      N, W, all.combos, combined.intensity, combined.meta), simplify = T,      cl = cl)
# # min       lq     mean  median       uq      max neval
# # 3.939423 3.956582 4.035985 4.07579 4.083981 4.122176    10
# # 4.273183 4.279005 4.293323 4.28712 4.303191 4.327932    10
# 
# # clusters = dynamicTreeCut::cutreeDynamic(dendro=clust, distM=1-all.combos[,1:nrow(all.combos)], method='hybrid', minClusterSize=0, deepSplit=0)
# # names(clusters)=clust$labels
# # sort(clusters)
# # plot.cors = all.combos[clust$order,clust$order]
# # names(clusters)=paste0('mz_',names(clusters))
# # p<-ggcorrplot::ggcorrplot(plot.cors)
# # p$data$facets_x <- clusters[as.character(p$data$Var2)]
# # p$data$facets_y <- clusters[as.character(p$data$Var1)]
# # p$coordinates <- coord_cartesian()
# # p + facet_grid(facets_x~facets_y, scales = "free", switch = "both",space='free')
# # 
# # combined.meta$combined_cluster_1 = rowMeans(pmin(pmax(combined.intensity[,names(clusters[clusters==1])],(-3)),3))
# # ggplot(combined.meta,aes(x=x,y=y,fill=combined_cluster_5,color=combined_cluster_5))+geom_tile()+facet_wrap(~combined.meta$sample)+scale_fill_gradient2()+scale_color_gradient2()
# 
# num.genes = 10
# svm.peaks = Reduce(union,lapply(FUN = function(x)head(x$gene,num.genes),X=svm_identification))
# 
# 
# svm_identification_sub = lapply(svm_identification,function(x)x[x$gene%in%svm.peaks,])
# svm_identification_sub = dplyr::bind_rows(svm_identification_sub,.id = 'sample')
# svm_identification_sub_square = as.data.frame(tidyr::pivot_wider(svm_identification_sub,id_cols = 'sample',names_from = 'gene',values_from = 'cor'))
# rownames(svm_identification_sub_square)=svm_identification_sub_square$sample
# svm_identification_sub_square = as.data.frame(t(svm_identification_sub_square[,2:ncol(svm_identification_sub_square)]))
# peaks = rownames(svm_identification_sub_square)
# peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
# mat <- svm_identification_sub_square
# mat[] <- peaks
# heatmaply::heatmaply_cor(svm_identification_sub_square,limits = c(min(svm_identification_sub_square),1),node_type='scatter',point_size_mat=svm_identification_sub_square,custom_hovertext = mat)
# 
# peaks = rownames(spatial.cross.cor)
# peaks = stringr::str_wrap(anno[match(peaks,anno$m_z),]$name,30)
# peaks.combo = expand.grid(peaks,peaks,stringsAsFactors = F)
# peaks.combo = paste0('row:',peaks.combo$Var1,',\ncol:',peaks.combo$Var2)
# mat <- spatial.cross.cor
# mat[] <- peaks.combo
# a = heatmaply::heatmaply_cor(spatial.cross.cor,custom_hovertext = mat,k_col=10)
# d <- dist(spatial.cross.cor, method = "euclidean")
# clusters = cutree(hclust(d), k = 10)
# dend <- as.dendrogram(hclust(d, method = "complete"))
# dend <- dendextend::seriate_dendrogram(dend, d)
# rownames(spatial.cross.cor[order.dendrogram(dend),])
# clusters = clusters[rev(rownames(spatial.cross.cor[order.dendrogram(dend),]))]
# cluster.names = names(clusters)
# clusters = as.numeric(factor(clusters,levels=unique(clusters)))
# names(clusters)=cluster.names
# intensity.scale = scale(intensity.matrix)
# for (cluster in unique(clusters)){
#   cluster.peaks = names(clusters[clusters==cluster])
#   metadata$combined_cluster = rowMeans(intensity.scale[,cluster.peaks])
# }
# heatmaply::heatmaply_cor(spatial.cross.cor,custom_hovertext = mat,k_col=10,k_row=10)
# ggplot(metadata,aes(x=x,y=y,color=combined_cluster))+geom_tile()+scale_color_gradient2()+facet_wrap(~metadata$Group,scale='free')
