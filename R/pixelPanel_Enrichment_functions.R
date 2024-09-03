make_control_samples <- function(control.samples,bulk.intensity.matrix,anno){
  anno = anno[anno$m_z%in%rownames(bulk.intensity.matrix),]
  anno = anno[anno$m_z!='NA',]
  return(bulk.intensity.matrix[anno$m_z,control.samples])
}


get_differential_peaks <- function(x,y,intensity.sub,metadata.sub,control.matrix,gap.x,gap.y,min.neighbours=3,pvalue.cutoff = 0.05,lfc.cutoff = 0.5,test='t-test',anno){

  # find the neighbours and check there are enough
  neighbours = metadata.sub[metadata.sub$x>(x-2*gap.x) &
                              metadata.sub$x<(x+2*gap.x) &
                              metadata.sub$y>(y-2*gap.y) &
                              metadata.sub$y<(y+2*gap.y),]
  if (nrow(neighbours)>3){
    intensity.sub = intensity.sub[metadata.sub$x>(x-2*gap.x) &
                                    metadata.sub$x<(x+2*gap.x) &
                                    metadata.sub$y>(y-2*gap.y) &
                                    metadata.sub$y<(y+2*gap.y),]
    intensity.sub = as.data.frame(intensity.sub)
    rownames(intensity.sub)=paste0('comparison_',1:nrow(intensity.sub))
    intensity.sub = t(intensity.sub)
    comparison.matrix = cbind(control.matrix,intensity.sub)
    de.table <- DEanalysis(intensity.matrix = comparison.matrix,
                           condition = c(rep('control',ncol(control.matrix)),rep('comparison',ncol(intensity.sub))),
                           var1 = 'control',
                           var2 = 'comparison',
                           test = test,
                           anno = anno)
    de.table = de.table[de.table$pvalAdj<pvalue.cutoff &
                          abs(de.table$lfc)>lfc.cutoff,]

    return(de.table)
  }
}
get_enriched_pathways_perpixel = function(row,intensity.sub,metadata.sub,control.table,
                                          min.neighbours=3,
                                          pvalue.cutoff = 0.05,
                                          lfc.cutoff = 0.5,
                                          test='t-test',
                                          background = NULL,
                                          min_path_size = 3,
                                          ora_pvalue_cutoff = 1.01,
                                          min_pathway_hits = 2,gap,
                                          anno,
                                          organism){

  de.table = get_differential_peaks(metadata.sub[row,'x'],metadata.sub[row,'y'],
                                    intensity.sub,
                                    metadata.sub,
                                    control.table,
                                    gap.x=gap,
                                    gap.y=gap,
                                    anno = anno)
  current.row = c(x=metadata.sub[row,'x'],y=metadata.sub[row,'y'],'significant'=paste(de.table$m_z,collapse = ','))
  if (!is.null(de.table)){
    ora <- execute_ora(de_peaks = de.table,
                       path_dict = NULL,
                       background = background,
                       min_path_size = min_path_size,
                       ora_pvalue_cutoff = ora_pvalue_cutoff,
                       min_pathway_hits = min_pathway_hits,
                       organism = organism,
                       anno = anno)
    if (!is.null(ora)){
      if (nrow(ora)!=0){
        ora$x = metadata.sub[row,'x']
        ora$y = metadata.sub[row,'y']
        return(ora)
      }
      else{
        return(NULL)
      }}else {
        return(NULL)
      }
  }
}

run_pixellevel_pipeline_parallel <- function(intensity.matrix,
                                             bulk.intensity.matrix,
                                             metadata,
                                             anno,
                                             control.samples,
                                             min.neighbours=3,
                                             pvalue.cutoff = 0.05,
                                             lfc.cutoff = 0.5,
                                             test='t-test',
                                             background = NULL,
                                             min_path_size = 3,
                                             ora_pvalue_cutoff = 1.01,
                                             min_pathway_hits = 2,
                                             ncores = 1,
                                             gcd.values,
                                             organism){
  control.table = make_control_samples(control.samples,bulk.intensity.matrix,anno)
  listofoutputs = list()
  message('Calculating pixel-level enrichment scores:')
  for (sample in unique(metadata$Group)){
    message(paste0('Calculating pixel-level enrichment scores for ',sample))
    # gap.x = find_pixel_size(metadata,'x')
    # gap.y = find_pixel_size(metadata,'y')
    gap.x = gap.y = gcd.values[[sample]]
    intensity.sub = intensity.matrix[metadata$Group==sample,]
    metadata.sub = metadata[metadata$Group==sample,]
    #  pb = txtProgressBar(min = 0, max = length(rownames(metadata)), initial = 0,style=3)
    cl <- parallel::makeCluster(ncores)
    clusterExport(cl, c("get_differential_peaks", "execute_ora", "metadata.sub","control.table","gap.x","DEanalysis","get_enriched_pathways_perpixel","anno","organism","kegg_db","get_ORA","intensity.sub"),envir = environment())
    outlist = pblapply(rownames(metadata.sub),FUN = function(x)get_enriched_pathways_perpixel(x,intensity.sub = intensity.sub,metadata.sub = metadata.sub,control.table = control.table,
                                                                                              min.neighbours=3,
                                                                                              pvalue.cutoff = 0.05,
                                                                                              lfc.cutoff = 0.5,
                                                                                              test='t-test',
                                                                                              background = NULL,
                                                                                              min_path_size = 3,
                                                                                              ora_pvalue_cutoff = 1.01,
                                                                                              min_pathway_hits = 2,
                                                                                              gap=gap.x,
                                                                                              anno=anno,
                                                                                              organism),cl=cl)
    names(outlist)=metadata.sub$spot_id
    outlist = outlist[-which(sapply(outlist, is.null))]
    outdf = dplyr::bind_rows(outlist, .id = "spot_id")
    outdf = merge(outdf,metadata.sub)
    stopCluster(cl)
    listofoutputs[[sample]]=outdf
  }
  return(listofoutputs)
}





