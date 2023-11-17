plot_pca <- function(
    expression.matrix,
    metadata,
    annotation.id,
    n.abundant = NULL,
    show.labels = FALSE,
    show.ellipses = TRUE,
    show.confidence.ellipses = FALSE,
    label.force = 1
){
  annotation.name <- colnames(metadata)[annotation.id]
  n.abundant <- min(n.abundant, nrow(expression.matrix))
  
  expr.PCA.list <- expression.matrix |>
    as.data.frame() |>
    dplyr::filter(seq_len(nrow(expression.matrix)) %in%
                    utils::tail(order(rowSums(expression.matrix)), n.abundant)) |>
    t()
  expr.PCA.list <- expr.PCA.list[, apply(expr.PCA.list, 2, function(x) max(x) != min(x))] %>%
    stats::prcomp(center = TRUE, scale = TRUE)
  
  expr.PCA <- dplyr::mutate(
                  as.data.frame(expr.PCA.list$x),
                  name = factor(metadata[, 1], levels = metadata[, 1]),
                  condition = if(!is.factor(metadata[,annotation.id])){
                    factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
                    }else{metadata[,annotation.id]}
                  
  )
  if(min(table(metadata[, annotation.id])) <= 2){
    expr.PCA.2 <- expr.PCA
    expr.PCA.2$PC1 <- expr.PCA.2$PC1 * 1.001
    expr.PCA.2$PC2 <- expr.PCA.2$PC2 * 1.001
    expr.PCA.full <- rbind(expr.PCA, expr.PCA.2)
  }
  else {expr.PCA.full <- expr.PCA}
  pca.plot <- ggplot2::ggplot(expr.PCA.full, ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$condition)) +
    ggplot2::theme_minimal() +
    ggplot2::geom_point() +
    ggplot2::labs(x = paste0("PC1 (proportion of variance = ", summary(expr.PCA.list)$importance[2, 1] * 100, "%)"),
         y = paste0("PC2 (proportion of variance = ", summary(expr.PCA.list)$importance[2, 2] * 100, "%)"),
         colour = annotation.name)
  if(show.confidence.ellipses){
    pca.plot <- pca.plot +
      stat_ellipse(geom='polygon',alpha=0.3,aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
    
  } else if(show.ellipses){
    pca.plot <- pca.plot +
      ggforce::geom_mark_ellipse(aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
  }
  if(show.labels){
    pca.plot <- pca.plot +
      ggrepel::geom_label_repel(
        data = expr.PCA,
        mapping = aes(x = .data$PC1, y = .data$PC2, colour = .data$condition, label = .data$name),
        max.overlaps = nrow(expr.PCA),
        force = label.force,
        point.size = NA
      )
  }
  
  pca.plot
}

perform_plsda <- function(expression.matrix,
                          metadata,
                          separator.id
){
  return(mixOmics::plsda(t(expression.matrix),as.vector(metadata[,separator.id]) , ncomp = 2))
}

plot_plsda <- function(
    expression.matrix,
    metadata,
    separator.id,
    annotation.id,
    n.abundant = NULL,
    show.labels = FALSE,
    show.ellipses = TRUE,
    show.confidence.ellipses = FALSE,
    label.force = 1
){
  annotation.name <- colnames(metadata)[annotation.id]
  my.plsda <- perform_plsda(expression.matrix,metadata,separator.id)
  coords = my.plsda$variates$X
  expr.plsda <- dplyr::mutate(
    as.data.frame(coords),
    name = factor(metadata[, 1], levels = metadata[, 1]),
    condition = if(!is.factor(metadata[,annotation.id])){
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))}else{metadata[,annotation.id]}
    
  )
  plsda.plot <- ggplot2::ggplot(expr.plsda, ggplot2::aes(x = .data$comp1, y = .data$comp2, colour = .data$condition)) +
    ggplot2::theme_minimal() +
    ggplot2::geom_point() +
    ggplot2::labs(x = paste0("PLS-DA Comp1 (proportion of variance = ", round(my.plsda$prop_expl_var$X[1] * 100,digits = 1), "%)"),
         y = paste0("PLS-DA Comp2 (proportion of variance = ", round(my.plsda$prop_expl_var$X[2] * 100,digits=1), "%)"),
         colour = annotation.name)
  if(show.confidence.ellipses){
    plsda.plot <- plsda.plot +
      ggplot2::stat_ellipse(geom='polygon',alpha=0.3,aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
    
  } else if(show.ellipses){
    plsda.plot <- plsda.plot +
      ggforce::geom_mark_ellipse(aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
  }
  if(show.labels){
    plsda.plot <- plsda.plot +
      ggrepel::geom_label_repel(
        data = expr.plsda,
        mapping = aes(x = .data$comp1, y = .data$comp2, colour = .data$condition, label = .data$name),
        max.overlaps = nrow(expr.plsda),
        force = label.force,
        point.size = NA
      )
  }
  plsda.plot
}


plsda_contrib <- function(expression.matrix,
                          metadata,
                          separator.id,
                          comp=1,
                          anno){
  my.plsda <- perform_plsda(expression.matrix,metadata,separator.id)
  contrib = data.frame(my.plsda$loadings$X)
  contrib$metab = rownames(contrib)
  colnames(contrib)[comp]='PLSDAComp'
  contrib = contrib[order(-abs(contrib$PLSDAComp)),]
  contrib = data.frame(contrib)
  contrib$display_metab = anno$display_name[match(contrib$metab,anno$m_z)]
  contrib$metab = factor(contrib$metab,levels=rev(contrib$metab))
  contrib.plot = ggplot2::ggplot(head(contrib,30),ggplot2::aes(x=PLSDAComp,y=metab,fill=PLSDAComp>0))+
    ggplot2::geom_bar(stat='identity') +
    ggplot2::theme_minimal()+
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab('Contribution')+
    ggplot2::ylab('')
  return(contrib.plot)
}

peaks_barplot <- function(sub.expression.matrix,
                          log.transformation = TRUE,
                          condition.vector){
  if (log.transformation){
    log.expression.matrix <- data.frame(log2(as.matrix(sub.expression.matrix) + 1))
  } else {
    log.expression.matrix <- sub.expression.matrix
  }
  log.expression.matrix$peak <- stringr::word(rownames(log.expression.matrix),sep='_',1,1)
  melted.expression.matrix <- tidyr::pivot_longer(log.expression.matrix,
                                                  cols = colnames(log.expression.matrix)[1:(ncol(log.expression.matrix)-1)])
  melted.expression.matrix$condition <- rep(condition.vector,nrow(sub.expression.matrix))
  p <- ggplot2::ggplot(melted.expression.matrix, ggplot2::aes(x = .data$name,
                                            y = .data$value,
                                            fill = .data$condition)) +
    ggplot2::geom_bar(stat='identity',position='dodge') +
    ggplot2::theme_minimal() +
    ggplot2::ylab(ifelse(log.transformation,'log2 intensity','intensity')) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),legend.position = "bottom") +
    ggplot2::facet_wrap(~melted.expression.matrix$peak,scales='free')
  p
}

peaks_boxplot <- function(sub.expression.matrix,
                          metadata,
                          log.transformation = TRUE,
                          metadata.column = 2){
  if (log.transformation){
    log.expression.matrix <- data.frame(log2(as.matrix(sub.expression.matrix) + 1))
  } else {
    log.expression.matrix <- sub.expression.matrix
  }
  log.expression.matrix$peak <- stringr::word(rownames(log.expression.matrix),sep='_',1,1)
  melted.expression.matrix <- tidyr::pivot_longer(log.expression.matrix,
                                                  cols = colnames(log.expression.matrix)[1:(ncol(log.expression.matrix)-1)])
  melted.expression.matrix$metadata = rep(metadata[,metadata.column],nrow(sub.expression.matrix))
  melted.expression.matrix$sample = rep(metadata[,1],nrow(sub.expression.matrix))
  p <- ggplot2::ggplot(melted.expression.matrix, ggplot2::aes(fill = metadata,
                                            x = metadata,
                                            y = value)) +
    ggplot2::geom_boxplot() + 
    ggplot2::geom_jitter(width=0.1) +
    ggplot2::theme_minimal() +
    ggplot2::ylab(ifelse(log.transformation,'log2 intensity','intensity')) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) + 
    ggplot2::facet_wrap(~melted.expression.matrix$peak,scales = 'free')
  p
  return.list = list('plot'=p,'table'=melted.expression.matrix)
}
