create_bulk_exp_regions <- function(expression.matrix,
                            metadata,
                            bulk.metadata,
                            sample.id.metadata.column,
                            region.ids,
                            minimum.pixels = 3
) {
  # add extra checks
  print(sample.id.metadata.column)
  print("I'm working on it")
  print(head(metadata[,sample.id.metadata.column]))
  print(head(region.ids))

  print(str(expression.matrix))
  expression.matrix$sample <- metadata[,sample.id.metadata.column]

  expression.matrix$region <- region.ids
  print(head(expression.matrix$sample))
  print(head(expression.matrix$region))
  expression.matrix.mean <- expression.matrix |>
    dplyr::group_by(sample,region) |>
    dplyr::summarise(across(everything(), mean),n=n())
  expression.matrix.mean = expression.matrix.mean[expression.matrix.mean$n>minimum.pixels,]
  print(expression.matrix.mean[1:5,1:5])
  print(summary(expression.matrix.mean$n))
  print(expression.matrix.mean[1:5,(ncol(expression.matrix.mean)-5):ncol(expression.matrix.mean)])
  print(summary(expression.matrix.mean$n))
  sample.names = paste0(expression.matrix.mean$sample,'_',expression.matrix.mean$region)
  print(sample.names)
  expression.matrix.mean <- expression.matrix.mean |>
    dplyr::select(-dplyr::one_of(c('sample','region'))) |>
    as.matrix() |>
    t() |>
    as.data.frame() |>
    dplyr::rename_with(~sample.names)
  print(expression.matrix.mean[1:5,1:5])
  for (region in unique(region.ids)){
    print(region)
    current.metadata = bulk.metadata
    current.metadata = data.frame(lapply(current.metadata, function(x) paste(x,region, sep="_")))
    print(current.metadata)
    if (region == unique(region.ids)[1]){
      metadata.mean = current.metadata
    } else {
      metadata.mean = rbind(metadata.mean,current.metadata)
    }
  }
  print(metadata.mean)
  metadata.mean = metadata.mean[metadata.mean[,1]%in%colnames(expression.matrix.mean),]
  return(list('expression.matrix'=expression.matrix.mean,'metadata'=metadata.mean))
}
