create_bulk_exp_regions <- function(expression.matrix,
                            metadata,
                            bulk.metadata,
                            sample.id.metadata.column,
                            region.ids,
                            minimum.pixels = 3
) {
  # add extra checks
  print(head(metadata))
  expression.matrix$sample <- metadata[,sample.id.metadata.column]
  # passed here
  expression.matrix$region <- region.ids
  expression.matrix.mean <- expression.matrix |>
    dplyr::group_by(sample,region) |>
    dplyr::summarise(across(everything(), mean),n=n())
  expression.matrix.mean = expression.matrix.mean[expression.matrix.mean$n>minimum.pixels,]
  sample.names = paste0(expression.matrix.mean$sample,'_',expression.matrix.mean$region)
  expression.matrix.mean <- expression.matrix.mean[,!(colnames(expression.matrix.mean)%in%c('sample','region'))]
  expression.matrix.mean <- expression.matrix.mean |>
    dplyr::select(-dplyr::one_of(c('sample','region'))) |>
    as.matrix() |>
    t() |>
    as.data.frame() |>
    dplyr::rename_with(~sample.names)
  for (region in unique(region.ids)){
    current.metadata = bulk.metadata
    current.metadata = data.frame(lapply(current.metadata, function(x) paste(x,region, sep="_")))
    if (region == unique(region.ids)[1]){
      metadata.mean = current.metadata
    } else {
      metadata.mean = rbind(metadata.mean,current.metadata)
    }
  }
  metadata.mean = metadata.mean[metadata.mean[,1]%in%colnames(expression.matrix.mean),]
  return(list('expression.matrix'=as.matrix(expression.matrix.mean),'metadata'=metadata.mean))
}
