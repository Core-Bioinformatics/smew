DEanalysis <-function(expression.matrix,condition,var1,var2,test='t-test',anno){

  # calculate FC
  expression.matrix = expression.matrix[(matrixStats::rowMins(as.matrix(expression.matrix))!=matrixStats::rowMaxs(as.matrix(expression.matrix))) &
                                          ((matrixStats::rowMins(as.matrix(expression.matrix[,condition==var1]))!=matrixStats::rowMaxs(as.matrix(expression.matrix[,condition==var1]))) |
                                          (matrixStats::rowMins(as.matrix(expression.matrix[,condition==var2]))!=matrixStats::rowMaxs(as.matrix(expression.matrix[,condition==var2])))),]
  fc_results=as.data.frame(expression.matrix)
  fc_results$log2_intensity = log2(rowMeans(expression.matrix))
  fc_results$group1_mean = as.numeric(rowMeans(expression.matrix[, condition==var1]))
  fc_results$group2_mean = as.numeric(rowMeans(expression.matrix[, condition==var2]))
  fc_results$fc = fc_results$group1_mean / fc_results$group2_mean
  fc_results$log2fc = log2(fc_results$fc)
  # drop NA
  na_index = apply(is.na(fc_results), 1, any)
  fc_results = fc_results[!na_index,]
  if (test=='t-test'){
    fc_results$pvalue = sapply(1:nrow(fc_results), function(i) t.test(as.numeric(fc_results[i, condition==var1]), as.numeric(fc_results[i, condition==var2]), paired = FALSE)$p.value)
  } else {
    fc_results$pvalue = sapply(1:nrow(fc_results), function(i) wilcox.test(as.numeric(fc_results[i, condition==var1]), as.numeric(fc_results[i, condition==var2]), paired = FALSE)$p.value)
  }
  fc_results$pvalue_adj = p.adjust(fc_results$pvalue, method = "BH")
  names = rownames(fc_results)
  # convert all columns to numeric (do not use tidyverse)
  fc_results = as.data.frame(apply(fc_results, 2, as.numeric))
  rownames(fc_results)=names
  # check if there are any NA values
  if (any(is.na(fc_results$pvalue))){
    print("NA values found in pvalue column")
  }
  if (any(is.na(fc_results$pvalue_adj))){
    print("NA values found in pvalue_adj column")
  }
  # drop NA values
  print(paste("Number of rows before dropping NA values:", nrow(fc_results)))
  fc_results = fc_results[apply(fc_results, 1, function(x) !any(is.na(x))),]
  print(paste("Number of rows after dropping NA values:", nrow(fc_results)))

  # drop duplicates
   print(paste("Number of rows before dropping duplicates:", nrow(fc_results)))
   fc_results = dplyr::distinct(fc_results)
   print(paste("Number of rows after dropping duplicates:", nrow(fc_results)))
  # rearrange columns
  fc_results = fc_results[,c('group1_mean', 'group2_mean', 'pvalue', 'pvalue_adj', 'log2fc','log2_intensity')] # , colnames_1, colnames_2
  print(paste("#### Number of rows after rearranging columns:", nrow(fc_results)))
  mytable = tibble::tibble('m_z'=rownames(fc_results),
                           'pval'=fc_results$pvalue,
                           'pvalAdj'=fc_results$pvalue_adj,
                           'lfc'=fc_results$log2fc,
                           'log2_intensity'=fc_results$log2_intensity)
  mytable = merge(mytable,anno[,colnames(anno)!='name'],all.x=T)
  return(mytable)
}
