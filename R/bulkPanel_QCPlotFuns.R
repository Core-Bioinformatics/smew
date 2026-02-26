#' QC plotting helpers
#' @description Utilities for plotting PCA/PLS-DA and other QC plots used in the Bulk QC panel.
#' @keywords internal
plot_pca <- function(
    pca.res,
    metadata,
    annotation.id,
    show.confidence.ellipses = TRUE,
    label.force = 1
){
  annotation.name <- colnames(metadata)[annotation.id]

  # expr.PCA.list <- intensity.matrix |>
  #   as.data.frame() |>
  #   t()
  # expr.PCA.res <- expr.PCA.list[, apply(expr.PCA.list, 2, function(x) max(x) != min(x))] |>
  #   stats::prcomp(center = TRUE, scale = TRUE)
  expr.PCA <- dplyr::mutate(
                  as.data.frame(pca.res$x),
                  "name" = factor(metadata[, 1], levels = metadata[, 1]),
                  "condition" = if(!is.factor(metadata[,annotation.id])){
                    factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
                    } else {metadata[,annotation.id]}

  )
  pca.plot <- ggplot2::ggplot(expr.PCA, ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$condition, label = .data$name)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("PC1 (proportion of variance = ", summary(pca.res)$importance[2, 1] * 100, "%)"),
         y = paste0("PC2 (proportion of variance = ", summary(pca.res)$importance[2, 2] * 100, "%)"),
         colour = annotation.name)
  if(show.confidence.ellipses){
    pca.plot <- pca.plot +
      ggplot2::stat_ellipse(geom='polygon',alpha=0.3,ggplot2::aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)

  }
  pca.plot <- pca.plot + ggplot2::geom_point()
  list('plot'=pca.plot)
}

pca_contrib <- function(pca.res,
                        comp=1,
                        anno){
  contrib = as.data.frame(pca.res$rotation)
  contrib$metab = rownames(contrib)
  colnames(contrib)[comp]='PCAComp'
  contrib = contrib[order(-abs(contrib$PCAComp)),]
  contrib = data.frame(contrib)
  contrib$display_metab = stringr::str_wrap(anno$display_name[match(contrib$metab,anno$m_z)],30)
  contrib$metab = factor(contrib$metab,levels=rev(contrib$metab))
  contrib.plot = ggplot2::ggplot(utils::head(contrib,30),ggplot2::aes(x=.data$PCAComp,y=.data$metab,fill=.data$PCAComp>0,label=.data$display_metab))+
    ggplot2::geom_bar(stat='identity') +
    ggplot2::theme_minimal()+
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab('Contribution')+
    ggplot2::ylab('')
  return(list('plot'=contrib.plot,'table'=utils::head(contrib,30)))
}

perform_plsda <- function(intensity.matrix,
                          metadata,
                          separator.id
){
  return(mixOmics::plsda(t(intensity.matrix),as.vector(metadata[,separator.id]) , ncomp = 2))
}

plot_plsda <- function(
    intensity.matrix,
    metadata,
    separator.id,
    annotation.id,
    show.confidence.ellipses = TRUE,
    label.force = 1
){
  annotation.name <- colnames(metadata)[annotation.id]
  my.plsda <- perform_plsda(intensity.matrix,metadata,separator.id)
  coords = my.plsda$variates$X
  expr.plsda <- dplyr::mutate(
    as.data.frame(coords),
    "name" = factor(metadata[, 1], levels = metadata[, 1]),
    "condition" = if(!is.factor(metadata[,annotation.id])){
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
    } else {
      metadata[,annotation.id]
    }
  )
  plsda.plot <- ggplot2::ggplot(expr.plsda, ggplot2::aes(x = .data$comp1, y = .data$comp2, colour = .data$condition, label = .data$name)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("PLS-DA Comp1 (proportion of variance = ", round(my.plsda$prop_expl_var$X[1] * 100,digits = 1), "%)"),
         y = paste0("PLS-DA Comp2 (proportion of variance = ", round(my.plsda$prop_expl_var$X[2] * 100,digits=1), "%)"),
         colour = annotation.name)
  if(show.confidence.ellipses){
    plsda.plot <- plsda.plot +
      ggplot2::stat_ellipse(geom='polygon',alpha=0.3,ggplot2::aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
  }
  plsda.plot <- plsda.plot + ggplot2::geom_point()
  plsda.plot
}


plsda_contrib <- function(intensity.matrix,
                          metadata,
                          separator.id,
                          comp=1,
                          anno){
  my.plsda <- perform_plsda(intensity.matrix,metadata,separator.id)
  contrib = data.frame(my.plsda$loadings$X)
  contrib$metab = rownames(contrib)
  colnames(contrib)[comp]='PLSDAComp'
  contrib = contrib[order(-abs(contrib$PLSDAComp)),]
  contrib = data.frame(contrib)
  contrib$display_metab = stringr::str_wrap(anno$display_name[match(contrib$metab,anno$m_z)],30)
  contrib$metab = factor(contrib$metab,levels=rev(contrib$metab))
  contrib.plot = ggplot2::ggplot(utils::head(contrib,30),ggplot2::aes(x=.data$PLSDAComp,y=.data$metab,fill=.data$PLSDAComp>0,label=.data$display_metab))+
    ggplot2::geom_bar(stat='identity') +
    ggplot2::theme_minimal()+
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab('Contribution')+
    ggplot2::ylab('')
  return(list('plot'=contrib.plot,'table'=utils::head(contrib,30)))
}

peaks_barplot <- function(sub.intensity.matrix,
                          log.transformation = TRUE,
                          condition.vector){
  if (log.transformation){
    log.intensity.matrix <- data.frame(log2(as.matrix(sub.intensity.matrix) + 1))
  } else {
    log.intensity.matrix <- sub.intensity.matrix
  }
  log.intensity.matrix$peak <- stringr::word(rownames(log.intensity.matrix),sep='_',1,2)
  melted.intensity.matrix <- tidyr::pivot_longer(log.intensity.matrix,
                                                  cols = colnames(log.intensity.matrix)[1:(ncol(log.intensity.matrix)-1)])
  melted.intensity.matrix$condition <- rep(condition.vector,nrow(sub.intensity.matrix))
  p <- ggplot2::ggplot(melted.intensity.matrix, ggplot2::aes(x = .data$name,
                                            y = .data$value,
                                            fill = .data$condition)) +
    ggplot2::geom_bar(stat='identity',position='dodge') +
    ggplot2::theme_minimal() +
    ggplot2::ylab(ifelse(log.transformation,'log2 intensity','intensity')) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),legend.position = "bottom") +
    ggplot2::facet_wrap(~.data$peak,scales='free')
  p
}

peaks_boxplot <- function(sub.intensity.matrix,
                          metadata,
                          log.transformation = TRUE,
                          metadata.column = 2){
  if (log.transformation){
    log.intensity.matrix <- data.frame(log2(as.matrix(sub.intensity.matrix) + 1))
  } else {
    log.intensity.matrix <- sub.intensity.matrix
  }
  log.intensity.matrix$peak <- stringr::word(rownames(log.intensity.matrix),sep='_',1,2)
  melted.intensity.matrix <- tidyr::pivot_longer(log.intensity.matrix,
                                                  cols = colnames(log.intensity.matrix)[1:(ncol(log.intensity.matrix)-1)])
  melted.intensity.matrix$metadata = rep(metadata[,metadata.column],nrow(sub.intensity.matrix))
  melted.intensity.matrix$sample = rep(metadata[,1],nrow(sub.intensity.matrix))
  p <- ggplot2::ggplot(melted.intensity.matrix, ggplot2::aes(fill = .data$metadata,
                                            x = .data$metadata,
                                            y = .data$value,label=.data$sample)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width=0.1,color='black',fill='black') +
    ggplot2::theme_minimal() +
    ggplot2::ylab(ifelse(log.transformation,'log2 intensity','intensity')) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
    ggplot2::facet_wrap(~.data$peak,scales = 'free')
  p
  return.list = list('plot'=p,'table'=melted.intensity.matrix)
}
