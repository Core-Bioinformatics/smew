compute_adduct_weights = function (adduct_table, mw) {
  # adduct table should contain only adducts for which the mw should be calculated
  ion.name <- adduct_table$Ion_Name
  ion.mass <- adduct_table$Ion_Mass
  mass.list <- as.list(ion.mass)
  mass.user <- lapply(mass.list, function(x) eval(parse(text = paste(gsub("PROTON", 1.00727646677, x)))))
  mass.df = as.data.frame(mass.user)
  colnames(mass.df) = ion.name
  mass.df = cbind(mw, mass.df)

  return(mass.df)
}

# function to extract compounds and corresponding associated adducts
get_matched_comps = function(adducts_table, mz_list, allowed_dppm) {

  # calculate the tolerance for all peaks
  tolerance = mz_list*allowed_dppm*1e-06

  # get unique entries from KEGG peaks
  adducts_dedup = dplyr::distinct(adducts_table)

  long_adducts_table = adducts_dedup |> pivot_longer(cols = -c(mw), names_to = 'adduct', values_to = 'adduct_mass')

  # get a subtraction matrix where rows are experimental peaks (mz_list) and columns are theoretical masses + adduct
  # values are the difference between the two
  diff_df = data.frame(outer(mz_list, long_adducts_table$adduct_mass, "-"))

  # rename columns to be the theoretical mass + adduct + adduct mass
  colnames(diff_df) = paste(long_adducts_table$mw, long_adducts_table$adduct, long_adducts_table$adduct_mass, sep = '_')

  # add the experimental peak and tolerance columns
  diff_df$exp_peak = mz_list
  diff_df$allowed_tolerance = tolerance

  # convert to long format
  diff_df_long = pivot_longer(diff_df, cols = -c(exp_peak, allowed_tolerance), names_to = 'mw_adduct_adductmass', values_to = 'observed_difference')
  # Melt the data table - this is a faster alternative to pivot_longer
#  diff_df_long = reshape2::melt(diff_df, id.vars = c("exp_peak", "allowed_tolerance"))
  # # Rename the variable columns
  colnames(diff_df_long) = c("exp_peak", "allowed_tolerance", "mw_adduct_adductmass", "observed_difference")

  # get absolute difference
  diff_df_long$abs_difference = abs(diff_df_long$observed_difference)

  # check if the difference is within the allowed tolerance
  diff_df_long$match = diff_df_long$abs_difference < diff_df_long$allowed_tolerance

  # keep only the matches
  diff_df_long = diff_df_long[diff_df_long$match,]

  # separate the theoretical mass, adduct and adduct mass
  diff_df_long = tidyr::separate_wider_delim(diff_df_long, cols = mw_adduct_adductmass, names = c('theoretical_mass', 'adduct', 'theoretical_mass_w_adduct'), delim = '_')

  # keep only the columns we need
  diff_df_long = diff_df_long[,c('exp_peak', 'allowed_tolerance', 'observed_difference', 'theoretical_mass', 'adduct', 'theoretical_mass_w_adduct')]

  # convert to numeric
  diff_df_long = diff_df_long |> dplyr::mutate(
    exp_peak = as.numeric(exp_peak),
    allowed_tolerance = signif(as.numeric(allowed_tolerance), 7),
    observed_difference = signif(as.numeric(observed_difference), 7),
    theoretical_mass = as.numeric(theoretical_mass),
    theoretical_mass_w_adduct = as.numeric(theoretical_mass_w_adduct)
  )

  return(diff_df_long)

}



get_matched_peaks = function(
    kegg_db = NULL, peak_list = NULL, #path_size = 2,
    ppm = 5, mode = NULL,
    adducts = NULL, neg_adduct_formulas = NULL, pos_adduct_formulas = NULL){

  if (is.null(kegg_db)) {
    message('Provide kegg db')
    stop()
  }
  if (is.null(peak_list)) {
    message('Provide the list of peaks')
    stop()
  }

  if (is.null(ppm)) {
    message('Provide ppm')
    stop()
  }
#  print(paste0('Using ppm = ', ppm))

  if (is.null(mode)) {
    message('Provide mode')
    stop()
  }
  if (!mode %in% c('Negative', 'Positive')) {
    message('Mode has to be either "Positive" or "Negative".')
    stop()
  }
#  print(paste0('Using mode = ', mode))
  exp_peak_list = unique(as.numeric(peak_list))
  # select adducts of interest
  if (mode == 'Negative') {
    my_adduct_formulas = neg_adduct_formulas[neg_adduct_formulas$Ion_Name %in% adducts,]

  } else if (mode == 'Positive') {
    my_adduct_formulas = pos_adduct_formulas[pos_adduct_formulas$Ion_Name %in% adducts,]
  }

  # keep only data for pathways with at least path_size compounds under it
  filt_kegg_db = kegg_db #%>% filter(pathway_id %in% filter(plyr::count(kegg_db$pathway_id), freq >= path_size)$x)

  # drop rows without given mass
  nbefore = n_distinct(filt_kegg_db$compound_id)
  filt_kegg_db = filt_kegg_db |> tidyr::drop_na(complete_compound_mass)
  nafter = n_distinct(filt_kegg_db$compound_id)

#  print(paste('Dropped', nbefore - nafter, 'compounds without mass.'))

#  print('Computing adduct masses...')
  # compute theoretical masses for all compounds+adducts in the kegg dataset
  adducts_computed = compute_adduct_weights(my_adduct_formulas,filt_kegg_db$complete_compound_mass)

#  print('Computing adduct masses...done')

#  print('Matching peaks to computed adducts...')

  # match experimental peaks to calculated adducts
  matched_comps = get_matched_comps(adducts_table = adducts_computed,
                                    mz_list = exp_peak_list,
                                    allowed_dppm = ppm)
#  print('Matching peaks to computed adducts...done')

  # merge with kegg data for complete annotation (join on theoretical mass)
  matched_comps_annot = left_join(matched_comps, filt_kegg_db, by = c('theoretical_mass' = 'complete_compound_mass'))

  # convert theoretical_mass to numeric
  matched_comps_annot$theoretical_mass = as.numeric(matched_comps_annot$theoretical_mass)

#  print(paste("Unique peaks mapped:", n_distinct(matched_comps_annot$exp_peak)))
  # print(paste(
  #   "Unique compounds mapped:",
  #   n_distinct(matched_comps_annot$compound_id)
  # ))
  # print(paste(
  #   "Unique compounds in LIPIDMAPS:",
  #   dplyr::n_distinct(
  #     matched_comps_annot |> dplyr::filter(is_lm == 'True') |> dplyr::pull(compound_id)
  #   )
  # ))

  # NOTE: some compounds have exactly the same molecular formula, hence they match the same peak.
  # This can inflate some pathway representation.

  # add mode column
  matched_comps_annot$mode = mode

  # reorder columns
  matched_comps_annot = matched_comps_annot[, c(
    'exp_peak',
    'allowed_tolerance',
    'observed_difference',
    'theoretical_mass_w_adduct',
    'theoretical_mass',
    'adduct',
    'mode',
    'compound_id',
    'compound_name',
    'compound_synonyms',
    'compound_formula',
    'is_hmdb',
    'is_lm',
    'lm_id',
    'lm_name',
    'lm_abbrev',
    'lm_formula',
    'pathway_id',
    'pathway_name',
    'human_pathway',
    'rat_pathway',
    'mouse_pathway'
  )]

  return(matched_comps_annot)
}


#' Create bulk dataset from pixel-wise intensity matrix
#'
#' @param intensity.matrix A data.frame containing m/z values on the columns and pixels/spots on the rows.
#' @param metadata A data.frame with the first column matching the rownames from intensity.matrix, the next 2
#' columns containing x and y coordinates and other columns describing attributes of the pixel, e.g. sample
#' of origin, timepoint, treatment, pathology annotations.
#' @param sample.id.column Name or index of column in metadata table which specifies the sample of origin
#' @param sample.wide.columns Name or indices of columns in metadata table other than sample.id.column which
#' contain sample-wide information, such as timepoint or treatment
#' @returns A numeric vector.
#' @examples
#' add(1, 1)
#' add(10, 1)
create_bulk_exp <- function(intensity.matrix,
                            metadata,
                            sample.id.column = 4,
                            sample.wide.columns = c(),
                            organism = 'Mouse',
                            adducts = c('M-H [1-]','M-H20-H [1-]','M+Cl [1-]'),
                            mode = 'Negative',
                            ppm = 5
                            ) {
  # add extra checks
#  intensity.matrix = t(unique(t(intensity.matrix)))
  intensity.matrix$sample <- metadata[,sample.id.column]
  intensity.matrix.mean <- intensity.matrix |>
          dplyr::group_by(sample) |>
          dplyr::summarise(across(everything(), mean))
  sample.names = intensity.matrix.mean$sample
  intensity.matrix.mean <- intensity.matrix.mean |>
          dplyr::select(-sample) |>
          as.matrix() |>
          t() |>
          as.data.frame() |>
          dplyr::rename_with(~sample.names)
  metadata.mean = unique(metadata[,c(sample.id.column,sample.wide.columns)])
  if (nrow(metadata.mean)!=ncol(intensity.matrix.mean)){
    stop('Please check all your sample.wide.columns are indeed sample-wide.')
  }

  #finished current stuff
  matched_peaks = get_matched_peaks(kegg_db = kegg_db,
                                    peak_list = gsub('X','',rownames(intensity.matrix.mean)),
                                    ppm = ppm,
                                    mode = mode,
                                    adducts = adducts,
                                    neg_adduct_formulas = neg_adduct_table,
                                    pos_adduct_formulas = pos_adduct_table)
  if (organism == 'Human'){
    matched_peaks = matched_peaks[matched_peaks$human_pathway == 'True',]
  } else if (organism == 'Mouse'){
    matched_peaks = matched_peaks[matched_peaks$mouse_pathway == 'True',]
  } else if (organism == 'Rat'){
    matched_peaks = matched_peaks[matched_peaks$rat_pathway == 'True',]
  } else {
    message('Organism supplied is not supported, no filtering applied')
    return(NULL)
  }
  annotation_table = matched_peaks[,c(1,6,8,9)]
  colnames(annotation_table)=c('m_z','adduct','kegg_id','name')
#  annotation_table$m_z = paste0('X',annotation_table$m_z)
  full_table = data.frame('m_z'=rownames(intensity.matrix.mean))
  annotation_table = merge(full_table,annotation_table,all.x=T)
  annotation_table = unique(annotation_table)
  annotation_table = annotation_table |>
    dplyr::group_by(m_z) |>
    dplyr::summarise(adduct = paste(adduct,collapse = ', '),
              kegg_id = paste(kegg_id,collapse = ', '),
              name = paste(name,collapse = ', '))
  annotation_table$display_name = ifelse(is.na(annotation_table$name),annotation_table$m_z,
                                               paste0(annotation_table$m_z,'_',annotation_table$name))
  return(list('intensity_matrix'=intensity.matrix.mean,
              'metadata'=metadata.mean,
              'matched_peaks'=matched_peaks,
              'annotation_table'=annotation_table))
}
