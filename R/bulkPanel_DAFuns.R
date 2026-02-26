#' Differential Analysis Helper Functions
#'
#' Internal functions supporting differential intensity analysis computation,
#' including statistical testing (t-test and Wilcox rank-sum), fold-change
#' calculation, and result formatting.
#'
#' @keywords internal
#' @name DAFuns
NULL

#' Perform Differential Analysis on Intensity Matrix
#'
#' @description Computes log2 fold-change and statistical test p-values comparing
#' two groups in an intensity matrix. Filters to varying features and handles
#' NA values and duplicate removal.
#'
#' @param intensity_matrix Numeric matrix of intensity values (peaks × samples).
#' @param condition Factor or character vector of group assignments (length = ncol(intensity_matrix)).
#' @param var1 Character; first group to compare.
#' @param var2 Character; second group to compare.
#' @param test Character; statistical test: "t-test" or "wilcox" (default: "t-test").
#' @param anno Data frame with peak annotation including 'm_z' and other columns.
#'
#' @return Data frame with columns:
#'   - `m_z`: Peak m/z value (from rownames)
#'   - `pval`: Unadjusted p-value
#'   - `pvalAdj`: BH-adjusted p-value
#'   - `lfc`: log2 fold-change (group1/group2)
#'   - `log2_intensity`: log2-transformed average intensity across all samples
#'   - Additional columns from `anno`
#'
#' @details
#' Only peaks with variance in the overall matrix OR variance within at least
#' one group are included. NA values are removed prior to analysis.
#'
#' @keywords internal
#' @export
de_analysis <- function(intensity_matrix, condition, var1, var2, test = 't-test', anno) {
  # ========================================================================
  # CALCULATE FOLD CHANGES
  # ========================================================================
  group1_names <- colnames(intensity_matrix[, condition == var1])
  group2_names <- colnames(intensity_matrix[, condition == var2])

  # Filter to peaks with variance in at least one group
  intensity_matrix <- intensity_matrix[
    (matrixStats::rowMins(as.matrix(intensity_matrix)) !=
       matrixStats::rowMaxs(as.matrix(intensity_matrix))) &
      ((matrixStats::rowMins(as.matrix(intensity_matrix[, condition == var1])) !=
          matrixStats::rowMaxs(as.matrix(intensity_matrix[, condition == var1]))) |
         (matrixStats::rowMins(as.matrix(intensity_matrix[, condition == var2])) !=
            matrixStats::rowMaxs(as.matrix(intensity_matrix[, condition == var2])))),
    ]

  fc_results <- as.data.frame(intensity_matrix)
  fc_results$log2_intensity <- log2(rowMeans(intensity_matrix) + 1)
  fc_results$group1_mean <- as.numeric(rowMeans(intensity_matrix[, condition == var1]))
  fc_results$group2_mean <- as.numeric(rowMeans(intensity_matrix[, condition == var2]))
  fc_results$fc <- fc_results$group1_mean / fc_results$group2_mean
  fc_results$lfc <- log2(fc_results$fc)

  # ========================================================================
  # Remove NA values
  # ========================================================================
  na_index <- apply(is.na(fc_results), 1, any)
  fc_results <- fc_results[!na_index, ]

  # ========================================================================
  # Statistical testing
  # ========================================================================
  if (test == 't-test') {
    fc_results$pvalue <- sapply(1:nrow(fc_results), function(i) {
      stats::t.test(as.numeric(fc_results[i, group1_names]),
             as.numeric(fc_results[i, group2_names]),
             paired = FALSE)$p.value
    })
  } else if (test == 'Wilcox rank sum') {
    fc_results$pvalue <- sapply(1:nrow(fc_results), function(i) {
      stats::wilcox.test(as.numeric(fc_results[i, group1_names]),
                  as.numeric(fc_results[i, group2_names]),
                  paired = FALSE)$p.value
    })
  }

  fc_results$pvalue_adj <- stats::p.adjust(fc_results$pvalue, method = "BH")
  names <- rownames(fc_results)

  # ========================================================================
  # Convert to numeric and clean up
  # ========================================================================
  fc_results <- as.data.frame(apply(fc_results, 2, as.numeric))
  rownames(fc_results) <- names

  # Filter to rows without NA values
  fc_results <- fc_results[apply(fc_results, 1, function(x) !any(is.na(x))), ]

  # Remove duplicate rows
  fc_results <- dplyr::distinct(fc_results)

  # ========================================================================
  # Format and merge with annotation
  # ========================================================================
  # Reorder columns for output
  fc_results <- fc_results[, c('group1_mean', 'group2_mean', 'pvalue',
                                'pvalue_adj', 'lfc', 'log2_intensity')]

  mytable <- tibble::tibble(
    'm_z' = rownames(fc_results),
    'pval' = fc_results$pvalue,
    'pvalAdj' = fc_results$pvalue_adj,
    'lfc' = fc_results$lfc,
    'log2_intensity' = fc_results$log2_intensity
  )

  mytable <- merge(mytable, anno[, colnames(anno) != 'name'], all.x = TRUE)
  return(mytable)
}
