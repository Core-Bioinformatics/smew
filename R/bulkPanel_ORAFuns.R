get_ORA = function (pathways, metabolites, universe, minSize = 1, maxSize = length(universe) - 1, direction = c("up", "down")) {

  if (!is.list(pathways)) {
    stop("pathways should be a list with each element containing metabolites from the universe")
  }

  if (any(duplicated(universe))) {
    warning("There were duplicate metabolites in universe, they were collapsed")
    universe = unique(universe)
  }

  empty_ora_result = data.frame(
    pathway = character(),
    total = numeric(),
    expected = numeric(),
    hits = numeric(),
    Raw.p = numeric(),
    Holm.p = numeric(),
    FDR = numeric(),
    metabolites = character(),
    direction = character()
  )

  minSize = max(minSize, 1)
  pathwaysFiltered = lapply(pathways, function(p) {
    unique(p[p%in%universe])
  })

  pathwaysSizes = sapply(pathwaysFiltered, length)

  toKeep = which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)

  if (length(toKeep) == 0) {
    return(empty_ora_result)
  }

  pathwaysFiltered = pathwaysFiltered[toKeep]
  pathwaysSizes = pathwaysSizes[toKeep]
  if (!all(metabolites %in% universe)) {
    warning("Not all of the input metabolites belong to the universe, such metabolites were removed")
  }

  metabolitesFiltered = unique(na.omit(metabolites[metabolites%in%universe]))
  if (length(metabolitesFiltered) == 0) {
    warning("No metabolites from the input list belong to the universe")
    return(empty_ora_result)
  }

  overlaps = lapply(pathwaysFiltered, intersect, metabolitesFiltered)
  overlap_metabolites = lapply(overlaps, function(x) paste0(x, collapse = '; '))
  overlapsT = data.frame(
    q = sapply(overlaps, length),
    m = sapply(pathwaysFiltered, length), # set.num
    n = length(universe) - sapply(pathwaysFiltered, length),
    k = length(metabolitesFiltered)) # q.size
  pathways_pvals = with(overlapsT, phyper(q - 1, m, n, k, lower.tail = FALSE))
  expected = with(overlapsT, k * (m/length(universe)))
  res = data.frame(
    pathway = names(pathwaysFiltered),
    total = overlapsT$m,
    expected = expected,
    hits = overlapsT$q,
    Raw.p = pathways_pvals,
    Holm.p = p.adjust(pathways_pvals, method = "holm"),
    FDR = p.adjust(pathways_pvals, method = "BH") ,
    metabolites = unlist(overlap_metabolites,use.names = F),
    direction = direction
  )
  res = res[order(res$Raw.p), ]
  return(res)
}


execute_ora = function(de_peaks, path_dict, background, min_path_size, ora_pvalue_cutoff, min_pathway_hits, peak_direction,organism, anno) {

  # if (background == 'Whole metabolome') {
  #   background_peaks = unique(unlist(path_dict))
  # } else if (background == 'Detected metabolome') {
  #   background_peaks = unique(de_peaks$theoretical_mass)
  # }
  # get up and down regulated peaks
  # if (peak_direction == 'Up'){
  #   my_peaks = unique(strsplit(paste(de_peaks[de_peaks$lfc>0,]$kegg_id,collapse=', '),split = ', ')[[1]])
  # } else if (peak_direction == 'Down') {
  #   my_peaks = unique(strsplit(paste(de_peaks[de_peaks$lfc<0,]$kegg_id,collapse=', '),split = ', ')[[1]])
  # } else {
  #   my_peaks = unique(strsplit(paste(de_peaks$kegg_id,collapse=', '),split = ', ')[[1]])
  # }

  # get up and down regulated peaks
  up_peaks = unique(strsplit(paste(de_peaks[de_peaks$lfc>0,]$kegg_id,collapse=', '),split = ', ')[[1]])
  down_peaks = unique(strsplit(paste(de_peaks[de_peaks$lfc<0,]$kegg_id,collapse=', '),split = ', ')[[1]])

  # if no peaks are up AND down regulated, return NULL and display message
  if (length(c(up_peaks,down_peaks)) == 0) {
    return(NULL)
  }
  if (organism == 'Human'){
    kegg_db_filt = kegg_db[kegg_db$human_pathway == 'True',]
  } else if (organism == 'Mouse'){
    kegg_db_filt = kegg_db[kegg_db$mouse_pathway == 'True',]
  } else if (organism == 'Rat'){
    kegg_db_filt = kegg_db[kegg_db$rat_pathway == 'True',]
  } else {
    print('Organism not supported')
    return(NULL)
  }
  # prepare a pathway:peak dictionary
  path_list = unique(kegg_db_filt$pathway_name)
  names(path_list) = path_list

  pathway2peaks = lapply(path_list, function (x) {
    a = unique(kegg_db_filt[kegg_db_filt$pathway_name == x,]$compound_id)
    a[!is.na(a)]
  })

  # keep only pathways with < 500 and > 1 entries
  pathway2peaks = pathway2peaks[which(lapply(pathway2peaks, length) < 500 & lapply(pathway2peaks, length) > 1)]
  ora_up = get_ORA(
    pathways = pathway2peaks,
    metabolites = up_peaks,
    universe = unique(strsplit(paste(anno$kegg_id,collapse=', '),split = ', ')[[1]]),
    minSize = min_path_size,
    maxSize = 500,
    direction = 'up'
  ) |> as.data.frame()

  ora_down = get_ORA(
    pathways = pathway2peaks,
    metabolites = down_peaks,
    universe = unique(strsplit(paste(anno$kegg_id,collapse=', '),split = ', ')[[1]]),
    minSize = min_path_size,
    maxSize = 500,
    direction = 'down'
  ) |> as.data.frame()

  ora_combined = rbind(ora_up,ora_down)
  # calculate fc, log2fc
  # ora_combined$FC = ora_combined$hits/ora_combined$expected
  # ora_combined$log2FC = log2(ora_combined$FC)
  #
  # # change the sign of the fc and log2fc for downregulated pathways
  # ora_combined$FC[ora_combined$direction == 'down'] = -ora_combined$FC[ora_combined$direction == 'down']
  # ora_combined$log2FC[ora_combined$direction == 'down'] = -ora_combined$log2FC[ora_combined$direction == 'down']
# could move this further up to make more efficient!
  if (min_pathway_hits) {
    ora_combined = ora_combined[ora_combined$hits >= min_pathway_hits,]
  }
  return(ora_combined)
}

ora_volcano_plot <- function(
    ORA_results,
    pval.threshold = 0.05,
    selectedPathways
){
  print(selectedPathways)
  df = ORA_results |>
    dplyr::mutate(log10pval = log10(.data$FDR),
                  lfc = log2(.data$hits/.data$expected)) |>
    dplyr::filter(!is.na(.data$log10pval))
  df$lfc = ifelse(df$direction=='up',df$lfc,-df$lfc)
  df$significance <- ifelse(df$FDR<pval.threshold,'Significant','Non-significant')
  df.label = df[df$pathway %in% selectedPathways,]
  lfc <- NULL; log10pval <- NULL; significance <- NULL
  vp <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = lfc, y = -log10pval,color=significance)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::xlab("log2(FC)") +
    ggplot2::ylab("-log10(pval)") +
    ggplot2::scale_color_manual(values=c("Non-significant"="#999999", "Significant"="#FF0000"))+
    ggrepel::geom_text_repel(data = df.label, mapping = aes(x = lfc, y = -log10pval,label = pathway))


  return(vp)
}

