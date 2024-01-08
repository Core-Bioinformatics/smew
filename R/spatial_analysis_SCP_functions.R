#install.packages("tidyverse")
library(tidyverse)

#install.packages("rio")
library(rio)

#install.packages("data.table")
library(data.table)

#install.packages("janitor")
library(janitor)

#install.packages("fields")
library(fields)

#install.packages("paletteer")
library(paletteer)

#install.packages("pbapply")
library(pbapply)

# functions - rewrite to be MAR independent and to incoporate aleks' code instead
MAR_run <- function(trial_MAR){ 
  
  x <- data.frame(new, trial_MAR)
  colnames(x)[1:2] <- c("m.z", "pvalue")
  
  export(x, "trial_MAR.csv")
  
  #rm(list = ls(all.names = TRUE))
  mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
  mSet <- SetPeakFormat(mSet, "mp")
  mSet <- UpdateInstrumentParameters(mSet, 5, "negative");
  mSet <- Read.PeakListData(mSet, "trial_MAR.csv");
  mSet <- SanityCheckMummichogData(mSet)
  add.vec <- c("M-H [1-]","M-H2O-H [1-]","M-H+O [1-]","M+Cl [1-]")
  mSet<-Setup.AdductData(mSet, add.vec);
  mSet<-PerformAdductMapping(mSet, "negative")
  mSet<-SetPeakEnrichMethod(mSet, "mum", "v1")
  mSet<-SetMummichogPval(mSet, 0.05)
  mSet<-PerformPSEA(mSet, "mmu_kegg", "current", 3 , 100)
  
  MAR_output <- import("mummichog_pathway_enrichment.csv")
  MAR_output <- MAR_output[1:40,]
  
  final_table <- data.frame(matrix(nrow = 40, ncol = 2))
  print(final_table)
  final_table[,1] <- MAR_output[,"V1"]
  print(final_table)
  # run this line if the required output is pval enrichment
  final_table[,2] <- MAR_output[1:40,"FET"]
  print(final_table)
  
  return(final_table)
}
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}
MAR_run_silent <- function(trial_MAR){
  hush(MAR_run(trial_MAR))
}


MAR_run_FC <- function(trial_MAR){ 
  
  x <- data.frame(new, trial_MAR)
  colnames(x)[1:2] <- c("m.z", "pvalue")
  
  export(x, "trial_MAR.csv")
  
  #rm(list = ls(all.names = TRUE))
  mSet <- InitDataObjects("mass_all", "mummichog", FALSE)
  mSet <- SetPeakFormat(mSet, "mp")
  mSet <- UpdateInstrumentParameters(mSet, 5, "negative");
  mSet <- Read.PeakListData(mSet, "trial_MAR.csv");
  mSet <- SanityCheckMummichogData(mSet)
  add.vec <- c("M-H [1-]","M-H2O-H [1-]","M-H+O [1-]","M+Cl [1-]")
  mSet<-Setup.AdductData(mSet, add.vec);
  mSet<-PerformAdductMapping(mSet, "negative")
  mSet<-SetPeakEnrichMethod(mSet, "mum", "v1")
  mSet<-SetMummichogPval(mSet, 0.05)
  mSet<-PerformPSEA(mSet, "mmu_kegg", "current", 3 , 100)
  
  MAR_output <- import("mummichog_pathway_enrichment.csv")
  MAR_output <- MAR_output[1:40,]
  
  final_table <- data.frame(matrix(nrow = 40, ncol = 2))
  print(final_table)
  final_table[,1] <- MAR_output[,"V1"]
  print(final_table)
  final_table[,2] <- MAR_output$Hits.sig/  MAR_output$Expected 
  print(final_table)
  
  return(final_table)
}
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}
MAR_run_FC_silent <- function(trial_MAR){
  hush(MAR_run_FC(trial_MAR))
}

read_input_table_from_API = function(path, file) {
  # Check if 'path' and 'file' are provided
  if (missing(path) || missing(file)) {
    stop("Both 'path' and 'file' must be specified.")
  }
  
  # Construct the full file path
  full_path <- file.path(path, file)
  
  # Check if the file exists
  if (!file.exists(full_path)) {
    stop("File not found: ", full_path)
  }
  
  # Read the table from the API using fread
  sp_control_table <- data.table::fread(full_path, nThread = 5)
  
  # Check if the table is empty
  if (nrow(sp_control_table) == 0) {
    warning("The table is empty.")
    return(NULL)
  }
  
  # Perform column deduplication
  sp_control_table <- sp_control_table[, .SD, .SDcols = unique(names(sp_control_table))]
  
  sp_control_table = data.frame(sp_control_table)
  
  return(sp_control_table)
  # EW: further potential checks to make sure the structure of df is correct? 
}
path = "/Users/da2395/Desktop/work/az"
file = "controls.csv"
read_input_table_from_API(path,file)


prepare_table_for_analysis =  function(sp_control_table){
  if (is.null(sp_control_table)) {
    stop("Input table is NULL.")
  }
  
  # Remove metadata columns
  sp_control_table <- sp_control_table[, 2:ncol(sp_control_table)]
  
  # Get rid of "x" in front of m/z value
  sp_control_table_filtered <- sp_control_table %>%
    select(contains("X", ignore.case = FALSE)) %>%
    setNames(str_sub(names(.), 2))
  
  # Get control ID
  sp_control_table_labels <- sp_control_table %>%
    select(ID)
  
  # Join tables
  sp_control_table_new <- dplyr::bind_cols(sp_control_table_labels, sp_control_table_filtered)
  
  # Pivot long
  sp_control_table_long <- sp_control_table_new %>%
    tidyr::pivot_longer(!ID, names_to = "mass", values_to = "intensity")
  
  return(sp_control_table_long)
}
prepare_table_for_analysis(file)


process_sp_control_table <- function(sp_control_table_long) {
  # Check if sp_control_table_long is a data frame
  if (!is.data.frame(sp_control_table_long)) {
    stop("Input must be a data frame.")
  }

  # Check if 'mass' and 'ID' columns exist
  if (!all(c("mass", "ID") %in% names(sp_control_table_long))) {
    stop("Input data frame must contain 'mass' and 'ID' columns.")
  }

  # Check if 'intensity' column exists
  if (!"intensity" %in% names(sp_control_table_long)) {
    stop("Input data frame must contain 'intensity' column.")
  }

  # Transform intensity column
  sp_control_table_long <- sp_control_table_long %>%
    mutate(log_intensity = log1p(intensity)) %>%
    select(-intensity) %>%
    gc()

  # Calculating averages per control
  averages_per_control <- sp_control_table_long %>%
    group_by(mass, ID) %>%
    dplyr::summarise(mean_log_intensity = mean(log_intensity), n = n()) %>%
    ungroup()

  # Round mass values
  averages_per_control$mass <- as.numeric(round(averages_per_control$mass, digits = 4))

  return(averages_per_control)
}

filter_data_by_day <- function(averages_per_control, days = c(7, 21)) {
  # Validate input
  stopifnot(is.data.frame(averages_per_control))

  # Convert specified days to character vector
  day_strings <- paste0("C_", days, "_", 1:3)

  # Filter data
  filtered_data <- averages_per_control %>%
    filter(ID %in% day_strings) %>%
    arrange(mass) %>%
    select(-n)

  return(filtered_data)
}
filter_data_by_day(averages_per_control, days = 7)


# picking out only peaks that match metabolites, using Aleks's metabolic pathway analysis package
source(paste0(path, '/metabolic-pathway-analysis-master@f5f1dffc620/functions.R'))
kegg_pathway_path = paste0(path,'metabolic-pathway-analysis-master@f5f1dffc620/Data/')
path_dict = get_path_dict(kegg_path = paste0(kegg_pathway_path, "complete_KEGG_LM_map_db.csv"), organism = 'Mouse')
neg_adduct_formulas_df = read.csv(paste0(kegg_pathway_path, "neg_adduct_table.csv"))
pos_adduct_formulas_df = read.csv(paste0(kegg_pathway_path, "pos_adduct_table.csv"))
matched_peaks = get_matched_peaks(kegg_path = paste0(kegg_pathway_path,"complete_KEGG_LM_map_db.csv"), 
                                  peak_list = unique(control_table_clean$mass), 
                                  ppm = 5, 
                                  mode = 'Negative', 
                                  adducts = c('M-H [1-]','M-H20-H [1-]','M+Cl [1-]'), 
                                  neg_adduct_formulas = neg_adduct_formulas_df, 
                                  pos_adduct_formulas = pos_adduct_formulas_df)

matched_peaks = matched_peaks[matched_peaks$mouse_pathway=='True',]
write.csv(matched_peaks,paste0(path,'matched_peaks_controls.csv'))

metabolites <- import(paste0(path,"matched_peaks_controls.csv"))
metabolites$m_z <- round(metabolites$exp_peak, digits=4)
metabolites <- unique(metabolites$m_z) 

my_table <- control_table_clean %>% 
  filter(mass %in% metabolites) %>% 
  spread(mass, "mean(log_intensity)") %>% 
  select(-ID)

print("Control table finished pre-processing")

# making individual sample tables
pixel_table <- data.table::fread(paste0(path, "B_21_6a.csv"),nThread = 5,)
pixel_table = pixel_table[, .SD, .SDcols = unique(names(pixel_table))]
pixel_table = data.frame(pixel_table)
data_sub <- pixel_table %>% select(c("x", "y", contains("X", ignore.case = F)))
names(data_sub)[3:length(data_sub)] <- str_sub(names(data_sub)[3:length(data_sub)],2)

# picking out only peaks that match metabolites
colnames(data_sub)[3:ncol(data_sub)]=as.character(round(as.numeric(colnames(data_sub)[3:ncol(data_sub)]),digits=4))
trial <- data_sub[,colnames(data_sub) %in% as.character(metabolites)]
trial <- log1p(trial)
trial <- trial[,!is.na(colnames(trial))]
trial <- trial %>%
  add_column(x = data_sub$x, .before = 1) %>%
  add_column(y = data_sub$y, .before = 2)
print("Sample table finished pre-processing")

# for loop for the statistical testing using the sliding window 
gc(); output_df <- data.frame("x" = 0, "y" = 0); for(i in colnames(trial[3:ncol(trial)])){
  print(which(i==colnames(trial)))
  # sub-setting tissue table to only get the required information
  subset <- trial %>% select(x, y, i) %>% data.table()
  
  # converts table to a matrix, to conserve the spatial information
  matrix <- as.matrix(dcast.data.table(data=subset, x ~ y, value.var=i, fill=NA)[,, with=FALSE]) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = 'x')
  
  # df to store the output information
  averages <- data.frame()
  
  # sliding window
  for(n in seq(1:nrow(matrix))){
    for(o in seq(1:ncol(matrix))){
      
      averages[noquote(paste0(n+(o-1)*nrow(matrix))), 'x'] <- as.numeric(rownames(matrix[n,]))
      averages[noquote(paste0(n+(o-1)*nrow(matrix))), 'y'] <- as.numeric(colnames(matrix[o]))
      
      matrix_sub <- matrix[ifelse(n-1 < 1, 1 , n-1):ifelse(n+1 > nrow(matrix), nrow(matrix), n+1),
                           ifelse(o-1 < 1, 1, o-1):ifelse(o+1 > ncol(matrix), ncol(matrix), o+1)]
      
      # if a pixel was background it was filled with NA in the matrix making step
      x <- na.omit(c(t(matrix_sub)))
      y <- as.numeric(my_table[,colnames(my_table) == i])
      
      # checking pixel has at least 3 non NA neighbours
      # tryCatch therefore required because of those NA values
      if (length(x)> 4){
        test <- tryCatch({wilcox.test(x, y)},
                         error = function(cond){return(NA)},
                         silent = T)
        
        averages[noquote(paste0(n+(o-1)*nrow(matrix))), paste(i)] <- tryCatch({test$p.value}, 
                                                                              error = function(cond){return(NA)},
                                                                              silent = T)
        averages[noquote(paste0(n+(o-1)*nrow(matrix))), 'median pixel intensity'] <- tryCatch({median(x)}, 
                                                                                              error = function(cond){return(NA)},
                                                                                              silent = T)
        averages[noquote(paste0(n+(o-1)*nrow(matrix))), 'median control intensity'] <- median(y)
        averages[noquote(paste0(n+(o-1)*nrow(matrix))), paste('log2fc',i)] <- log2(median(x)/median(y))
        
      }
    }
  }
  
  output_df <- full_join(output_df,  averages[,c('x','y',i)])
  write.csv(output_df, paste0(path, "output_df_C_7_3.csv"))
}
print("Finished running the slidding window on all pixels")

output_df <- output_df[2:nrow(output_df),2:ncol(output_df)]
row.has.na <- apply(output_df, 1, function(x){!all(is.na(x[3:length(x)]))})
output_df = output_df[row.has.na,]
output_df <- output_df %>% mutate(across(c(x,y), as.character))%>% unite("id", x:y, sep = " x ")

trial_MAR <- data.frame(t(output_df[,2:length(output_df)])) 
colnames(trial_MAR)=output_df$id
new = rownames(trial_MAR)

library(MetaboAnalystR)

output_final <- pblapply(trial_MAR, MAR_run_silent)
output_file_store <- output_final
x <- ldply(output_final, rbind)
x <- pivot_wider(x, names_from = .id, values_from = X2)

export(x, paste0(path, "pixel_MAR_output_C_7_3_pval.csv"))

# this iterates through each sample and gets the pvalues and highlights only the ones with p-value<0.1 (change to 0.05 if you want that cutoff)
for (my.sample in c('bleo_d21_1','bleo_d21_2','bleo_d21_3','bleo_d21_4a','bleo_d21_5b','bleo_d21_6a','bleo_d7_1','bleo_d7_4b',
                    'control_d21_1','control_d21_2','control_d21_3','control_d7_1'))
{
  print(my.sample)
  output_final <- import(paste0("pixel_MAR_output_",my.sample,"_pval.csv"))
  
  x <- output_final
  x[is.na(x)] <- 0
  
  image_prep <- data.frame(t(x)) %>% row_to_names(1) %>%
    rownames_to_column("coords") 
  
  image_prep[,1] <- sub(" x ", " ", image_prep[,1])
  image_prep <- image_prep %>% separate(coords, c("x", "y"),sep = ' ')
  
  image_prep <- sapply(image_prep, as.numeric)
  
  image_prep <- image_prep %>% data.table()
  image_prep[image_prep<=1&image_prep>0.1]<-1
  
  #colours for pvalue
  paletteer_c("ggthemes::Red-Black Diverging", 100)
  my_colours <- c("black", "#AE123AFF", "#B1163BFF", "#B4193BFF", "#B71C3CFF", "#BB1F3CFF", "#BE223DFF", "#C1243DFF", "#C4273EFF", "#C7293FFF", "#CB2C3FFF", "#CE2F40FF", "#D13140FF", "#D43341FF", "#D83641FF", "#DB3842FF", "#DE3B42FF", "#E13D43FF", "#E44144FF", "#E54546FF", "#E64A49FF", "#E84E4BFF", "#E9534EFF", "#EB5750FF", "#EC5B52FF", "#ED5F55FF", "#EE6357FF","#F0675AFF", "#F16B5CFF", "#F26F5FFF", "#F37261FF", "#F57663FF", "#F67A66FF", "#F77D68FF", "#F8816BFF", "#F78771FF", "#F78C78FF", "#F6927EFF", "#F59885FF", "#F49D8BFF", "#F2A292FF", "#F1A898FF", "#EFAD9FFF", "#EDB2A6FF", "#EBB8ACFF", "#E9BDB3FF", "#E7C2BAFF", "#E4C7C1FF", "#E1CCC8FF", "#DED1CFFF", "#DBD6D6FF", "#D7D7D7FF", "#D4D4D4FF", "#D0D1D1FF", "#CDCECEFF", "#C9CBCBFF", "#C6C8C8FF", "#C2C5C5FF", "#BFC2C2FF", "#BBBFBFFF", "#B8BCBDFF", "#B4B9BAFF", "#B1B6B7FF", "#AEB3B4FF", "#AAB0B1FF", "#A7ADAEFF", "#A3AAABFF", "#A0A7A8FF", "#9DA4A6FF", "#9AA2A3FF", "#979FA1FF", "#949C9FFF", "#919A9DFF", "#8E979AFF", "#8B9498FF", "#889296FF", "#868F94FF", "#838D91FF", "#808A8FFF", "#7D878DFF", "#7A858BFF", "#778289FF", "#748086FF", "#717D84FF", "#6F7B82FF", "#6C7880FF", "#6A757DFF", "#68737BFF", "#657079FF", "#636E76FF", "#606B74FF", "#5E6972FF", "#5C6670FF", "#59636DFF", "#57616BFF", "#555E69FF", "#525C67FF", "#505965FF", "#4E5762FF", "#4B5460FF", "#49525EFF") 
  
  # for loop, creates images of the tissues for each pathway
  pdf(paste0(path,'Pathways_',my.sample,'_pval0.1.pdf'))
  for(i in sort(colnames(image_prep)[3:ncol(image_prep)])){
    matrix <- as.matrix(dcast.data.table(data=image_prep, x ~ y, value.var=i, fill=0)[,, with=FALSE]) %>% as.data.frame() %>% column_to_rownames(var = 'x')
    matrix_num <- sapply(matrix, as.numeric)
    
    img <- imagePlot(matrix_num, 
                     col = my_colours
                     #breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3) #for FC
                     #breaks = c(0, 0.0000001, 0.05, 0.25, 0.5, 0.75, 1)
    )
    mtext(line=3, side=1, i, outer=F)
    
  }
  dev.off()
}

# making a big table with all the p-values for all the samples
my_colours <- c("#4B5460FF", "#AE123AFF", "#B1163BFF", "#B4193BFF", "#B71C3CFF", "#BB1F3CFF", "#BE223DFF", "#C1243DFF", "#C4273EFF", "#C7293FFF", "#CB2C3FFF", "#CE2F40FF", "#D13140FF", "#D43341FF", "#D83641FF", "#DB3842FF", "#DE3B42FF", "#E13D43FF", "#E44144FF", "#E54546FF", "#E64A49FF", "#E84E4BFF", "#E9534EFF", "#EB5750FF", "#EC5B52FF", "#ED5F55FF", "#EE6357FF","#F0675AFF", "#F16B5CFF", "#F26F5FFF", "#F37261FF", "#F57663FF", "#F67A66FF", "#F77D68FF", "#F8816BFF", "#F78771FF", "#F78C78FF", "#F6927EFF", "#F59885FF", "#F49D8BFF", "#F2A292FF", "#F1A898FF", "#EFAD9FFF", "#EDB2A6FF", "#EBB8ACFF", "#E9BDB3FF", "#E7C2BAFF", "#E4C7C1FF", "#E1CCC8FF", "#DED1CFFF", "#DBD6D6FF", "#D7D7D7FF", "#D4D4D4FF", "#D0D1D1FF", "#CDCECEFF", "#C9CBCBFF", "#C6C8C8FF", "#C2C5C5FF", "#BFC2C2FF", "#BBBFBFFF", "#B8BCBDFF", "#B4B9BAFF", "#B1B6B7FF", "#AEB3B4FF", "#AAB0B1FF", "#A7ADAEFF", "#A3AAABFF", "#A0A7A8FF", "#9DA4A6FF", "#9AA2A3FF", "#979FA1FF", "#949C9FFF", "#919A9DFF", "#8E979AFF", "#8B9498FF", "#889296FF", "#868F94FF", "#838D91FF", "#808A8FFF", "#7D878DFF", "#7A858BFF", "#778289FF", "#748086FF", "#717D84FF", "#6F7B82FF", "#6C7880FF", "#6A757DFF", "#68737BFF", "#657079FF", "#636E76FF", "#606B74FF", "#5E6972FF", "#5C6670FF", "#59636DFF", "#57616BFF", "#555E69FF", "#525C67FF", "#505965FF", "#4E5762FF", "#4B5460FF", "#ECECEC") 

for (my.sample in c('bleo_d21_1','bleo_d21_2','bleo_d21_3','bleo_d21_4a','bleo_d21_5b','bleo_d21_6a','bleo_d7_1','bleo_d7_4b',
                    'control_d21_1','control_d21_2','control_d21_3','control_d7_1')){
  print(my.sample)
  output_final <- import(paste0(path,"pixel_MAR_output_",my.sample,"_pval.csv"))
  
  x <- output_final
  x[is.na(x)] <- 0
  
  image_prep <- data.frame(t(x)) %>% row_to_names(1) %>%
    rownames_to_column("coords") 
  
  image_prep[,1] <- sub(" x ", " ", image_prep[,1])
  image_prep <- image_prep %>% separate(coords, c("x", "y"),sep = ' ')
  
  image_prep <- sapply(image_prep, as.numeric)
  #image_prep[image_prep<=1&image_prep>0.1]<-1
  image_prep_ongoing = as.data.frame(image_prep)
  image_prep_ongoing <- image_prep_ongoing[,order(colnames(image_prep_ongoing))]
  image_prep_ongoing$sample = my.sample
  print(str(image_prep_ongoing))
  if (my.sample=='bleo_d21_1'){
    overall_image_prep = image_prep_ongoing
  } else {
    overall_image_prep = overall_image_prep[,intersect(colnames(overall_image_prep),colnames(image_prep_ongoing))]
    image_prep_ongoing = image_prep_ongoing[,intersect(colnames(overall_image_prep),colnames(image_prep_ongoing))]
    overall_image_prep = rbind(overall_image_prep,image_prep_ongoing)
  }
}
overall_image_prep$sample = factor(overall_image_prep$sample,levels=c('control_d21_1','control_d21_2','control_d21_3','control_d7_1',
                                                                      'bleo_d21_1','bleo_d21_2','bleo_d21_3','bleo_d7_1','bleo_d21_4a',
                                                                      'bleo_d21_5b','bleo_d21_6a','bleo_d7_4b'
))

# for each pathway plot the overall p-val distribution and the significant ones under 2 thresholds
pdf(paste0(path,'data/new/AllPathways.pdf'),width=10,height=7)
for (pathway in colnames(overall_image_prep)[1:61]){
  print(pathway)
  overall_image_prep$pval = overall_image_prep[,pathway]
  print(ggplot(overall_image_prep,aes(x=x,y=y,fill=pval,color=pval))+
          geom_tile()+
          facet_wrap(~overall_image_prep$sample,ncol=4,scales = 'free')+scale_fill_gradientn(colors=my_colours,values = c(seq(0,1,0.01)))+
          scale_color_gradientn(colors=my_colours,values = c(seq(0,1,0.01)))+theme_classic()+ggtitle(pathway))
  
  overall_image_prep_sig = overall_image_prep
  overall_image_prep_sig[overall_image_prep_sig<=1&(overall_image_prep_sig>0.05)]<-1
  overall_image_prep_sig[overall_image_prep_sig<1&(overall_image_prep_sig==0)]<-1
  print(ggplot(overall_image_prep_sig,aes(x=x,y=y,fill=pval,color=pval))+
          geom_tile()+
          facet_wrap(~overall_image_prep_sig$sample,ncol=4,scales = 'free')+scale_fill_gradientn(colors=my_colours,values = c(seq(0,1,0.01)),na.value = '#ECECEC')+
          scale_color_gradientn(colors=my_colours,values = c(seq(0,1,0.01)),na.value = '#ECECEC')+theme_classic()+ggtitle(pathway))
  
  overall_image_prep_sig = overall_image_prep
  overall_image_prep_sig[overall_image_prep_sig<=1&(overall_image_prep_sig>0.1)]<-1
  overall_image_prep_sig[overall_image_prep_sig<1&(overall_image_prep_sig==0)]<-1
  
  print(ggplot(overall_image_prep_sig,aes(x=x,y=y,fill=pval,color=pval))+
          geom_tile()+
          facet_wrap(~overall_image_prep_sig$sample,ncol=4,scales = 'free')+scale_fill_gradientn(colors=my_colours,values = c(seq(0,1,0.01)),na.value = '#ECECEC')+
          scale_color_gradientn(colors=my_colours,values = c(seq(0,1,0.01)),na.value = '#ECECEC')+theme_classic()+ggtitle(pathway))
}
dev.off()
write.csv(overall_image_prep,paste0(path,'data/new/AllPathways_pvals.csv'))

# compare number significant and non-significant in each sample
for (pathway in colnames(overall_image_prep)[1:61]){
  a = as.data.frame(table((overall_image_prep[,pathway]<0.05&overall_image_prep[,pathway]!=0),overall_image_prep$sample))
  if (length(unique(a[,1]))!=1){
    a = tidyr::pivot_wider(a,id_cols = 'Var2',names_from = 'Var1',values_from = 'Freq')
    colnames(a)=c('sample','nonsignificant','significant')
    a$proprtion_significant = round(100*a$significant/(a$nonsignificant+a$significant),digits = 2)
    write.csv(a,paste0(path,'proportion_significant/',pathway,'_pval0.05.csv'),row.names=F)
  }
  a = as.data.frame(table(overall_image_prep[,pathway]<0.1&overall_image_prep[,pathway]!=0,overall_image_prep$sample))
  if (length(unique(a[,1]))!=1){
    a = tidyr::pivot_wider(a,id_cols = 'Var2',names_from = 'Var1',values_from = 'Freq')
    colnames(a)=c('sample','nonsignificant','significant')
    a$proprtion_significant = round(100*a$significant/(a$nonsignificant+a$significant),digits = 2)
    write.csv(a,paste0(path,'proportion_significant/',pathway,'_pval0.1.csv'),row.names=F)
  }
}

# make plots showing significance level of difference for a few pathways of interest
overall_image_prep = read.csv(paste0(path,'data/new/AllPathways_pvals.csv',row.names = 1))
colnames(overall_image_prep)=gsub('\\.',' ',colnames(overall_image_prep))
pdf(paste0(path,'PathwayBoxplots.pdf'),width=4,height=4)
for (pathway in c('Biosynthesis of unsaturated fatty acids','Citrate cycle  TCA cycle ','Fructose and mannose metabolism','Glyoxylate and dicarboxylate metabolism','Inositol phosphate metabolism','Phosphatidylinositol signaling system','Purine metabolism','Arachidonic acid metabolism','Ascorbate and aldarate metabolism')){
  a = as.data.frame(table((overall_image_prep[,pathway]<0.05&overall_image_prep[,pathway]!=0),overall_image_prep$sample))
  a = tidyr::pivot_wider(a,id_cols = 'Var2',names_from = 'Var1',values_from = 'Freq')
  colnames(a)=c('sample','nonsignificant','significant')
  a$proprtion_significant = round(100*a$significant/(a$nonsignificant+a$significant),digits = 2)
  a$type = stringr::word(a$sample,sep='_',1,1)
  print(ggplot(a,aes(x=type,y=proprtion_significant,fill=type))+geom_boxplot()+ylab('% significant')+theme_classic()+stat_compare_means()+ggtitle(pathway))
}
dev.off()
