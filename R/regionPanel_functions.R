create_bulk_exp_regions <- function(intensity.matrix,
                            metadata,
                            bulk.metadata,
                            sample.id.metadata.column,
                            region.ids,
                            minimum.pixels = 3
) {
  # add extra checks
  intensity.matrix$sample <- metadata[,sample.id.metadata.column]
  # passed here
  intensity.matrix$region <- region.ids
  intensity.matrix.mean <- intensity.matrix |>
    dplyr::group_by(sample,region) |>
    dplyr::summarise(across(everything(), mean),n=n())
  intensity.matrix.mean = intensity.matrix.mean[intensity.matrix.mean$n>minimum.pixels,]
  sample.names = paste0(intensity.matrix.mean$sample,'_',intensity.matrix.mean$region)
  intensity.matrix.mean <- intensity.matrix.mean[,!(colnames(intensity.matrix.mean)%in%c('sample','region'))]
  intensity.matrix.mean <- intensity.matrix.mean |>
    as.matrix() |>
    t() |>
    as.data.frame() |>
    dplyr::rename_with(~sample.names)
  for (region in unique(region.ids)){
    current.metadata = bulk.metadata
    current.metadata$AllSamples = 'AllSamples'
    current.metadata = data.frame(lapply(current.metadata, function(x) paste(x,region, sep="_")))
    if (region == unique(region.ids)[1]){
      metadata.mean = current.metadata
    } else {
      metadata.mean = rbind(metadata.mean,current.metadata)
    }
  }
  metadata.mean = metadata.mean[metadata.mean[,1]%in%colnames(intensity.matrix.mean),]
  intensity.matrix.mean = intensity.matrix.mean[,metadata.mean[,1]]
  return(list('intensity.matrix'=as.matrix(intensity.matrix.mean),'metadata'=metadata.mean))
}
