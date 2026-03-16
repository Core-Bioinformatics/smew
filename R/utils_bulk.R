utils::globalVariables(c(".data"))

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
bulk_utils_DA <- function(intensity_matrix, condition, var1, var2, test = 't-test', anno) {
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

#' Create a volcano plot visualising differential analysis (DA) results
#' @description This function creates a volcano plot to visualise the results
#' of a DA analysis.
#' @param peaks.de.results the table of DA peaks, usually generated by
#' \code{\link{bulk_utils_DA}}
#' @param pval.threshold,lfc.threshold the p-value and/or log2(fold-change)
#' thresholds to determine whether a peak is DA
#' @param alpha the transparency of points; ignored for DA peaks if
#' add.intensity.colour.gradient is TRUE; default is 0.1
#' @param xlims a single value to create (symmetric) x-axis limits; by default
#' inferred from the data
#' @param log10pval.cap whether to cap the log10(p-value at -10); any p-values
#' lower that 10^(-10) are set to the cap for plotting
#' @param add.colours whether to colour peaks based on their log2(fold-change)
#' and -log10(p-value); default is TRUE
#' @param add.intensity.colour.gradient whether to add a colour gradient
#' for DA peaks to present their log2(intensity); default is TRUE
#' @param add.guide.lines whether to add vertical and horizontal guide lines
#' to the plot to highlight the thresholds; default is TRUE
#' @param add.labels.auto whether to automatically label peaks with the
#' highest |log2(fold-change)| and intensity; default is TRUE
#' @param add.labels.custom whether to add labels to user-specified peaks;
#' the parameter peaks.to.label must also be specified; default is FALSE
#' @param ... parameters passed on to \code{\link{bulk_utils_volcano_enhance}}
#' @keywords internal
#' @return The volcano plot as a ggplot object.
bulk_utils_volcano_plot <- function(
    peaks.de.results,
    pval.threshold = 0.05,
    lfc.threshold = 1,
    alpha = 0.1,
    xlims = NULL,
    log10pval.cap = TRUE,
    add.colours = TRUE,
    add.intensity.colour.gradient = TRUE,
    add.guide.lines = TRUE,
    add.labels.auto = TRUE,
    add.labels.custom = FALSE,
    ...
){
  df = peaks.de.results |>
    dplyr::mutate(peak = .data$m_z, log10pval = log10(.data$pvalAdj), name = stringr::str_wrap(.data$display_name)) |>
    dplyr::filter(!is.na(.data$log10pval))

  if(all(df$log10pval >= -10)) log10pval.cap <- FALSE
  if(log10pval.cap) df$log10pval[df$log10pval < -10] <- -10

  lfc <- NULL; log10pval <- NULL
  vp <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = .data$lfc, y = -.data$log10pval,label=.data$name)) +
    ggplot2::theme_minimal() +
    ggplot2::xlab("log2(FC)") +
    ggplot2::ylab("-log10(pval)")

  if(is.null(xlims)){
    max.abs.lfc = max(abs(df[df$log10pval > -Inf,]$lfc))
    vp <- vp + ggplot2::xlim(-max.abs.lfc, max.abs.lfc)
  }else{
    vp <- vp + ggplot2::xlim(-abs(xlims), abs(xlims))
  }

  if(log10pval.cap){
    vp <- vp + ggplot2::scale_y_continuous(labels=c("0.0", "2.5", "5.0", "7.5", ">10"))
  }

  if(any(add.colours,
         add.intensity.colour.gradient,
         add.guide.lines,
         add.labels.auto,
         add.labels.custom)){
    vp <- bulk_utils_volcano_enhance(
      vp = vp,
      df = df,
      pval.threshold = pval.threshold,
      lfc.threshold = lfc.threshold,
      alpha = alpha,
      add.colours = add.colours,
      add.intensity.colour.gradient = add.intensity.colour.gradient,
      add.guide.lines = add.guide.lines,
      add.labels.auto = add.labels.auto,
      add.labels.custom = add.labels.custom,
      ...
    )
  }

  return(vp)

}

#' @description \code{\link{bulk_utils_volcano_enhance}} is called indirectly by
#' \code{\link{bulk_utils_volcano_plot}} to add extra features.
#' @param vp volcano plot as a ggplot object (usually passed by \code{\link{bulk_utils_volcano_plot}})
#' @param df data frame of DE results for all peaks (usually passed by
#' \code{\link{bulk_utils_volcano_plot}})
#' @param point.colours a vector of 4 colours to colour peaks with both pval
#' and lfc under thresholds, just pval under threshold, just lfc under threshold,
#' both pval and lfc over threshold (DE peaks) respectively; only used if
#' add.colours is TRUE
#' @param raster whether to rasterize non-DE peaks with ggraster to reduce
#' memory usage; particularly useful when saving plots to files
#' @param colour.gradient.scale a vector of two colours to create a colour
#' gradient for colouring the DA peaks based on intensity; a named list with
#' components left and right can be supplied to use two different colour scales;
#' only used if add.intensity.colour.gradient is TRUE
#' @param colour.gradient.breaks,colour.gradient.limits parameters to customise
#' the legend of the colour gradient scale; especially useful if creating
#' multiple plots or a plot with two scales;
#' only used if add.intensity.colour.gradient is TRUE
#' @param guide.line.colours a vector with two colours to be used to colour
#' the guide lines; the first colour is used for the p-value and log2(fold-change)
#' thresholds and the second for double those values
#' @param annotation annotation data frame containing a match between the peak
#' field of df (usually ENSEMBL IDs) and the peak names that should be shown
#' in the plot labels; not necessary if df already contains peak names
#' @param n.labels.auto a integer vector of length 3 denoting the number of
#' peaks that should be automatically labelled; the first entry corresponds to
#' DE peaks with the lowest p-value, the second to those with highest absolute
#' log2(fold-change) and the third to those with highest intensity; a single
#' integer can also be specified, to be used for all 3 entries; default is 5
#' @param peaks.to.label a vector of peak names to be labelled in the
#' plot; if names are present those are shown as the labels (but the values are
#' the ones matched - this is to allow custom peak names to be presented)
#' @param seed the random seed to be used for reproducibility; only used for
#' ggrepel::geom_label_repel if labels are present
#' @param label.force passed to the force argument of ggrepel::geom_label_repel;
#' higher values make labels overlap less (at the cost of them being further
#' away from the points they are labelling)
#' @return The enhanced volcano plot as a ggplot object.
#' @keywords internal
#' @rdname bulk_utils_volcano_plot
bulk_utils_volcano_enhance <- function(
    vp,
    df,
    pval.threshold,
    lfc.threshold,
    alpha,
    add.colours,
    point.colours = c("#bfbfbf", "orange", "red", "blue"),
    raster = FALSE,
    add.intensity.colour.gradient,
    colour.gradient.scale = list(left  = c("#99e6ff", "#000066"),
                                 right = c("#99e6ff", "#000066")),
    colour.gradient.breaks = ggplot2::waiver(),
    colour.gradient.limits = NULL,
    add.guide.lines,
    guide.line.colours = c("green", "blue"),
    add.labels.auto,
    add.labels.custom,
    annotation = NULL,
    n.labels.auto = c(5, 5, 5),
    peaks.to.label = NULL,
    seed = 0,
    label.force = 1
){

  logp.threshold = log10(pval.threshold)

  if(add.colours){
    colours = vector(length=nrow(df))
    colours[] = point.colours[1]
    colours[abs(df$lfc) > lfc.threshold] = point.colours[2]
    colours[df$log10pval < logp.threshold] = point.colours[3]
    colours[abs(df$lfc) > lfc.threshold & df$log10pval < logp.threshold] = point.colours[4]
    df$colours <- colours

    if(raster){
      vp <- vp + ggrastr::rasterise(ggplot2::geom_point(alpha = alpha, colour = colours))
    }else{
      vp <- vp + ggplot2::geom_point(alpha = alpha, colour = colours, fill = colours)
    }
  }

  if(add.intensity.colour.gradient){
    df.colour.gradient <- df |>
      dplyr::filter(abs(.data$lfc) > lfc.threshold & .data$log10pval < logp.threshold) |>
      dplyr::arrange(.data$log2_intensity)
    if(identical(colour.gradient.scale$left, colour.gradient.scale$right)){
      vp <- vp +
        ggplot2::geom_point(data = df.colour.gradient,
                   mapping = ggplot2::aes(x = .data$lfc, y = -.data$log10pval, colour = .data$log2_intensity)) +
        ggplot2::scale_color_gradient(low = colour.gradient.scale$left[1],
                             high = colour.gradient.scale$left[2],
                             breaks = colour.gradient.breaks,
                             limits = colour.gradient.limits) +
        ggplot2::labs(colour = "log2(intensity)")
    }else{
      vp <- vp +
        ggplot2::geom_point(data = dplyr::filter(df.colour.gradient, .data$lfc < 0),
                   mapping = ggplot2::aes(x = .data$lfc, y = -.data$log10pval, colour = .data$log2_intensity)) +
        ggplot2::scale_color_gradient(low = colour.gradient.scale$left[1],
                             high = colour.gradient.scale$left[2],
                             breaks = colour.gradient.breaks,
                             limits = colour.gradient.limits) +
        ggplot2::labs(colour = "log2(intensity)") +
        ggnewscale::new_scale_colour() +
        ggplot2::geom_point(data = dplyr::filter(df.colour.gradient, .data$lfc > 0),
                   mapping = ggplot2::aes(x = .data$lfc, y = -.data$log10pval, colour = .data$log2_intensity)) +
        ggplot2::scale_colour_gradient(low = colour.gradient.scale$right[1],
                              high = colour.gradient.scale$right[2],
                              breaks = colour.gradient.breaks,
                              limits = colour.gradient.limits) +
        ggplot2::labs(colour = "log2(intensity)")
    }
  }

  if(add.guide.lines){
    vp <- vp +
      ggplot2::geom_vline(xintercept = lfc.threshold,  colour = guide.line.colours[1]) +
      ggplot2::geom_vline(xintercept = -lfc.threshold,  colour = guide.line.colours[1]) +
      ggplot2::geom_vline(xintercept =  2 * lfc.threshold,  colour = guide.line.colours[2]) +
      ggplot2::geom_vline(xintercept = -2 * lfc.threshold,  colour = guide.line.colours[2]) +
      ggplot2::geom_hline(yintercept = -logp.threshold, colour = guide.line.colours[1]) +
      ggplot2::geom_hline(yintercept = -2 * logp.threshold, colour = guide.line.colours[2])
  }

  if(add.labels.auto | add.labels.custom){
    if(!is.null(annotation)){
      df <- df |>
        dplyr::mutate(
          symbol = .data$peak,
          name = ifelse(is.na(.data$symbol), .data$peak, .data$symbol)
        ) |>
        dplyr::select(-.data$symbol)
    }else{
      df <- df |> dplyr::mutate(name = .data$peak)
    }

    df.label <- tibble::tibble()
    if(add.labels.custom){
      peaks.to.rename <- peaks.to.label[names(peaks.to.label) != ""]
      peaks.to.label <- df$name[(match(peaks.to.label, c(df$name, df$peak)) - 1) %% nrow(df) + 1]
      peaks.to.label <- unique(peaks.to.label[!is.na(peaks.to.label)])
      peaks.to.rename <- peaks.to.rename[peaks.to.rename %in% peaks.to.label]
      df.label <- dplyr::filter(df, .data$name %in% peaks.to.label)
      df.label$name[match(peaks.to.rename, df.label$name)] <- names(peaks.to.rename)
      if(nrow(df.label) == 0){
        message(paste0("add.labels.custom was TRUE but no peaks specified; ",
                       "did you forget to supply peaks.to.label or annotation?"))
      }
    }

    if(add.labels.auto){
      if(length(n.labels.auto) == 1) n.labels.auto <- rep(n.labels.auto, 3)
      df.significant <- dplyr::filter(df, !(.data$name %in% peaks.to.label))

      df.significant <- df.significant[order(abs(df.significant$lfc), decreasing=TRUE), ]
      df.highest.lfc <- utils::head(df.significant, n.labels.auto[1])
      df.rest <- utils::tail(df.significant, nrow(df.significant) - n.labels.auto[1]) |>
        dplyr::filter(abs(.data$lfc) > lfc.threshold, .data$log10pval < logp.threshold)

      df.rest <- df.rest[order(abs(df.rest$log10pval), decreasing=TRUE), ]
      df.lowest.p.vals <- utils::head(df.rest, n.labels.auto[2])
      df.rest <- utils::tail(df.rest, nrow(df.rest) - n.labels.auto[2])

      df.rest <- df.rest[order(df.rest$log2_intensity, decreasing=TRUE), ]
      df.highest.abn <- utils::head(df.rest, n.labels.auto[3])

      df.label <- rbind(df.lowest.p.vals, df.highest.lfc, df.highest.abn, df.label) |>
        dplyr::distinct(.data$name, .keep_all = TRUE)
    }
  }

  return(vp)

}

#' Create an MA plot visualising differential analysis (DA) results
#' @description This function creates an MA plot to visualise the results
#' of a DA analysis.
#' @inheritParams bulk_utils_volcano_plot
#' @param ylims a single value to create (symmetric) y-axis limits; by default
#' inferred from the data
#' @param ... parameters passed on to \code{\link{bulk_utils_ma_enhance}}
#' @keywords internal
#' @return The MA plot as a ggplot object.
bulk_utils_ma_plot <- function(
    peaks.de.results,
    pval.threshold = 0.05,
    lfc.threshold = 1,
    alpha = 0.1,
    ylims = NULL,
    add.colours = TRUE,
    add.intensity.colour.gradient = TRUE,
    add.guide.lines = TRUE,
    add.labels.auto = TRUE,
    add.labels.custom = FALSE,
    ...
){
  df = peaks.de.results |>
    dplyr::mutate(peak = .data$m_z, log10pval = log10(.data$pvalAdj), name = stringr::str_wrap(.data$display_name)) |>
    dplyr::filter(!is.na(.data$log10pval))
  log2_intensity <- NULL; lfc <- NULL
  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = .data$log2_intensity, y = .data$lfc,label=.data$name)) +
    ggplot2::theme_minimal() +
    ggplot2::xlab("Average log2(intensity)") +
    ggplot2::ylab("log2(FC)")

  if(is.null(ylims)){
    max.abs.lfc = max(abs(df$lfc))
    p <- p + ggplot2::ylim(-max.abs.lfc, max.abs.lfc)
  }else{
    p <- p + ggplot2::ylim(-abs(ylims), abs(ylims))
  }

  if(any(add.colours,
         add.intensity.colour.gradient,
         add.guide.lines,
         add.labels.auto,
         add.labels.custom)){
    p <- bulk_utils_ma_enhance(
      p = p,
      df = df,
      pval.threshold = pval.threshold,
      lfc.threshold = lfc.threshold,
      alpha = alpha,
      add.colours = add.colours,
      add.intensity.colour.gradient = add.intensity.colour.gradient,
      add.guide.lines = add.guide.lines,
      add.labels.auto = add.labels.auto,
      add.labels.custom = add.labels.custom,
      ...
    )
  }

  return(p)

}

#' Add features to an MA plot visualising differential analysis (DA) results
#' @description \code{\link{bulk_utils_ma_enhance}} is called indirectly by
#' \code{\link{bulk_utils_ma_plot}} to add extra features.
#' @param p MA plot as a ggplot object (usually passed by \code{\link{bulk_utils_ma_plot}})
#' @param df data frame of DE results for all peaks (usually passed by
#' \code{\link{bulk_utils_ma_plot}})
#' @keywords internal
#' @return The enhanced MA plot as a ggplot object.
#' @rdname bulk_utils_ma_plot
bulk_utils_ma_enhance <- function(
    p,
    df,
    pval.threshold,
    lfc.threshold,
    alpha,
    add.colours,
    point.colours = c("#bfbfbf", "orange", "red", "blue"),
    raster = FALSE,
    add.intensity.colour.gradient,
    colour.gradient.scale = list(left  = c("#99e6ff", "#000066"),
                                 right = c("#99e6ff", "#000066")),
    colour.gradient.breaks = ggplot2::waiver(),
    colour.gradient.limits = NULL,
    add.guide.lines,
    guide.line.colours = c("green", "blue"),
    add.labels.auto,
    add.labels.custom,
    annotation = NULL,
    n.labels.auto = c(5, 5, 5),
    peaks.to.label = NULL,
    seed = 0,
    label.force = 1
){

  logp.threshold = log10(pval.threshold)

  if(add.colours){
    colours = vector(length=nrow(df))
    colours[] = point.colours[1]
    colours[abs(df$lfc) > lfc.threshold] = point.colours[2]
    colours[df$log10pval < logp.threshold] = point.colours[3]
    colours[abs(df$lfc) > lfc.threshold & df$log10pval < logp.threshold] = point.colours[4]
    df$colours <- colours

    if(raster){
      p <- p + ggrastr::rasterise(ggplot2::geom_point(alpha = alpha, colour = colours))
    }else{
      p <- p + ggplot2::geom_point(alpha = alpha, colour = colours, fill = colours)
    }
  }

  if(add.intensity.colour.gradient){
    df.colour.gradient <- df |>
      dplyr::filter(abs(.data$lfc) > lfc.threshold & .data$log10pval < logp.threshold) |>
      dplyr::arrange(.data$log2_intensity)
    if(identical(colour.gradient.scale$left, colour.gradient.scale$right)){
      p <- p +
        ggplot2::geom_point(data = df.colour.gradient,
                   mapping = ggplot2::aes(x = .data$log2_intensity, y = .data$lfc, colour = .data$log2_intensity)) +
        ggplot2::scale_color_gradient(low = colour.gradient.scale$left[1],
                             high = colour.gradient.scale$left[2],
                             breaks = colour.gradient.breaks,
                             limits = colour.gradient.limits) +
        ggplot2::labs(colour = "log2(intensity)")
    }else{
      p <- p +
        ggplot2::geom_point(data = dplyr::filter(df.colour.gradient, .data$lfc < 0),
                   mapping = ggplot2::aes(x = .data$log2_intensity, y = .data$lfc, colour = .data$log2_intensity)) +
        ggplot2::scale_color_gradient(low = colour.gradient.scale$left[1],
                             high = colour.gradient.scale$left[2],
                             breaks = colour.gradient.breaks,
                             limits = colour.gradient.limits) +
        ggplot2::labs(colour = "log2(intensity)") +
        ggnewscale::new_scale_colour() +
        ggplot2::geom_point(data = dplyr::filter(df.colour.gradient, .data$lfc > 0),
                   mapping = ggplot2::aes(x = .data$log2_intensity, y = .data$lfc, colour = .data$log2_intensity)) +
        ggplot2::scale_colour_gradient(low = colour.gradient.scale$right[1],
                              high = colour.gradient.scale$right[2],
                              breaks = colour.gradient.breaks,
                              limits = colour.gradient.limits) +
        ggplot2::labs(colour = "log2(intensity)")
    }
  }

  if(add.guide.lines){
    p <- p +
      ggplot2::geom_hline(yintercept =      lfc.threshold,  colour = guide.line.colours[1]) +
      ggplot2::geom_hline(yintercept =     -lfc.threshold,  colour = guide.line.colours[1]) +
      ggplot2::geom_hline(yintercept =  2 * lfc.threshold,  colour = guide.line.colours[2]) +
      ggplot2::geom_hline(yintercept = -2 * lfc.threshold,  colour = guide.line.colours[2])
  }

  if(add.labels.auto | add.labels.custom){
    if(!is.null(annotation)){
      df <- df |>
        dplyr::mutate(
          symbol = .data$peak,
          name = ifelse(is.na(.data$symbol), .data$peak, .data$symbol)
        ) |>
        dplyr::select(-.data$symbol)
    }else{
      df <- df |> dplyr::mutate(name = .data$peak)
    }

    df.label <- tibble::tibble()
    if(add.labels.custom){
      peaks.to.rename <- peaks.to.label[names(peaks.to.label) != ""]
      peaks.to.label <- df$name[(match(peaks.to.label, c(df$name, df$peak)) - 1) %% nrow(df) + 1]
      peaks.to.label <- unique(peaks.to.label[!is.na(peaks.to.label)])
      peaks.to.rename <- peaks.to.rename[peaks.to.rename %in% peaks.to.label]
      df.label <- dplyr::filter(df, .data$name %in% peaks.to.label)
      df.label$name[match(peaks.to.rename, df.label$name)] <- names(peaks.to.rename)
      if(nrow(df.label) == 0){
        message(paste0("add.labels.custom was TRUE but no peaks specified; ",
                       "did you forget to supply peaks.to.label or annotation?"))
      }
    }

    if(add.labels.auto){
      if(length(n.labels.auto) == 1) n.labels.auto <- rep(n.labels.auto, 3)
      df.significant <- dplyr::filter(df, !(.data$name %in% peaks.to.label))
      df.significant <- df.significant[order(abs(df.significant$lfc), decreasing=TRUE), ]
      df.highest.lfc <- utils::head(df.significant, n.labels.auto[1])
      df.rest <- utils::tail(df.significant, nrow(df.significant) - n.labels.auto[1]) |>
        dplyr::filter(abs(.data$lfc) > lfc.threshold, .data$log10pval < logp.threshold)

      df.rest <- df.rest[order(abs(df.rest$log10pval), decreasing=TRUE), ]
      df.lowest.p.vals <- utils::head(df.rest, n.labels.auto[2])
      df.rest <- utils::tail(df.rest, nrow(df.rest) - n.labels.auto[2])

      df.rest <- df.rest[order(df.rest$log2_intensity, decreasing=TRUE), ]
      df.highest.abn <- utils::head(df.rest, n.labels.auto[3])

      df.label <- rbind(df.lowest.p.vals, df.highest.lfc, df.highest.abn, df.label) |>
        dplyr::distinct(.data$name, .keep_all = TRUE)
    }
  }
  return(p)
}

#' Perform GRN Inference
#' @description Performs Metabolite Peak Regulatory Network inference on a subset of the
#'   intensity matrix for a set of target metabolite peaks of interest.
#' @param intensity_matrix Numeric matrix of peak intensities with peaks as rows
#'   and samples as columns.
#' @param metadata Data frame with sample metadata; rows must match ncol(intensity_matrix).
#' @param anno Data frame with peak annotation including 'display_name' (peak names)
#'   and 'm_z' (peak IDs).
#' @param seed Numeric seed for reproducibility (default: 13).
#' @param targets Character vector of target metabolite peak names/IDs to focus GRN on.
#' @param condition Character name of metadata column to subset samples by.
#' @param samples Character vector of condition values to include in analysis.
#' @param inference_method Character specifying GRN method; currently only "GENIE3"
#'   is supported.
#' @keywords internal
#' @return Adjacency matrix with regulators in rows, targets in columns.
#'   Values represent edge weights from GENIE3 algorithm.
#' @details When sample count exceeds 1000, a random subset of 1000 samples is used
#' for computational efficiency. The random seed is set for reproducibility.
bulk_utils_infer_GRN <- function(intensity_matrix, metadata, anno, seed = 13,
                      targets, condition, samples, inference_method) {
  inference_method <- inference_method[1]
  if (targets[1] %in% anno$display_name) {
    target_ids <- anno$m_z[match(targets, anno$display_name)]
  } else {
    target_ids <- targets
  }
  # Subset samples based on condition
  sample_indices <- metadata[[condition]] %in% samples
  intensity_sampled <- intensity_matrix[, sample_indices]
  if (inference_method == "GENIE3") {
    # Limit to 1000 samples for computational efficiency
    MAX_SAMPLES <- 1000
    if (ncol(intensity_sampled) > MAX_SAMPLES) {
      set.seed(seed)
      my_sample <- sample(1:ncol(intensity_sampled), MAX_SAMPLES)
    } else {
      my_sample <- 1:ncol(intensity_sampled)
    }
    set.seed(seed)
    res <- GENIE3::GENIE3(intensity_sampled[, my_sample], targets = target_ids)
  }
  res
}

#' Convert the adjacency matrix to network links
#' @description This function converts an adjacency matrix to a data frame
#' of network links, subset to the most important ones.
#' @param weightMat the (weighted) adjacency matrix - regulators in rows,
#' targets in columns
#' @param plotConnections the number of connections to subset to
#' @keywords internal
#' @return A data frame with fields from, to and value, describing the edges
#' of the network
bulk_utils_get_link_list_rename <- function(weightMat, plotConnections){
  GENIE3::getLinkList(weightMat, plotConnections) |>
    dplyr::mutate(from = as.character(.data$regulatoryGene), 
                  to = as.character(.data$targetGene), 
                  value = .data$weight, 
                  regulatoryGene = NULL, 
                  targetGene = NULL,
                  weight = NULL)
}

#' Find recurring regulators
#' @description This function finds regulators that appear as the 
#' same network edge in more than one of the input networks.
#' @inheritParams bulk_utils_get_link_list_rename
#' @param weightMatList a list of (weighted) adjacency matrices; 
#' each list element must be an adjacency matrix with regulators in rows,
#' targets in columns
#' @keywords internal
#' @return A vector containing the names of the recurring regulators
bulk_utils_find_regulators_with_recurring_edges <- function(weightMatList, plotConnections){
  edges <- lapply(weightMatList, function(wm){
    bulk_utils_get_link_list_rename(wm, plotConnections)
  }) |>
    dplyr::bind_rows() 
  edges |>
    dplyr::group_by(.data$from, .data$to) |>
    dplyr::summarise(n = dplyr::n(), .groups = "keep") |>
    dplyr::filter(.data$n > 1) |>
    dplyr::pull(.data$from) |>
    dplyr::setdiff(edges$to)
}

#' Plot a GRN
#' @description This function creates a network plot of a GRN.
#' @inheritParams bulk_utils_get_link_list_rename
#' @param anno a data frame with peak annotation including 'display_name' (peak names)
#' and 'm_z' (peak IDs); used to highlight target peaks in the network
#' @param plotConnections the number of connections to subset to for plotting
#' @param plot_position_grid,n_networks the position of the plot in 
#' the grid (1-4) and the number of networks shown (1-4); these are
#' solely used for hiding unwanted plots in the shiny app
#' @param recurring_regulators targets to be highlighted; usually the
#' result of \code{\link{bulk_utils_find_regulators_with_recurring_edges}}
#' @keywords internal
#' @return A network plot. See visNetwork package for more details.
bulk_utils_plot_GRN <- function(weightMat, anno, plotConnections, 
                     plot_position_grid, n_networks, recurring_regulators){
  
  if(n_networks >= plot_position_grid){
    edges <- bulk_utils_get_link_list_rename(weightMat, plotConnections)
    nodes <- tibble::tibble(
      id = c(edges$to, edges$from),
      label = c(edges$to, edges$from),
      group = rep(c("target", "regulator"), each = nrow(edges)),
      color = NA
    ) |>
      dplyr::distinct(.data$id, .keep_all = TRUE)
    nodes$group[nodes$id %in% recurring_regulators] <- "recurring_regulator"
    nodes$group = paste0(nodes$group,ifelse(nodes$id %in% anno$m_z,'',paste0('_othermod')))
    
    color_target <- list("background" = '#D2E5FF')
    color_regulator <- list("background" = '#E0E0E0')
    color_recurring_regulator <- list("background" = '#ACE9B4')
    color_target_mod = list("background" = "#a2326a")
    color_regulator_mod <- list("background" = '#dcb2dc')
    color_recurring_mod = list("background" = '#580058')
    
    colors.list <- list(color_regulator, color_target, color_recurring_regulator,
                        color_target_mod, color_regulator_mod, color_recurring_mod)
    
    for(i in seq_len(nrow(nodes))){
      nodes$color[i] <- colors.list[[match(
        nodes$group[i], c("regulator", "target", "recurring_regulator",
                          "regulator_othermod", "target_othermod", "recurring_regulator_othermod")
      )]]
    }
    
    visNetwork::visNetwork(nodes, edges)
  } else {
    NULL
  }
}

#' Visualise the overlap of edges between different networks
#' @description This function creates an UpSet plot of the intersections
#' and specific differences of the edges in the input networks.
#' @inheritParams bulk_utils_get_link_list_rename
#' @inheritParams bulk_utils_find_regulators_with_recurring_edges
#' @keywords internal
#' @return An UpSet plot. See UpSetR package for more details.
bulk_utils_plot_upset <- function(weightMatList, plotConnections){
  if(length(weightMatList) > 1){
    edge.list <- lapply(weightMatList, function(wm){
      bulk_utils_get_link_list_rename(wm, plotConnections) |>
        dplyr::mutate(connection = paste0(.data$from, "-", .data$to)) |>
        dplyr::pull(.data$connection)
    })
    names(edge.list) <- paste0("Network", seq_len(length(edge.list)))
    UpSetR::upset(UpSetR::fromList(edge.list), order.by = "degree", keep.order = TRUE)
  }else{
    NULL
  }
}

#' Over-Representation Analysis (ORA) Helper Function (Bulk)
#'
#' @description Computes pathway enrichment for a set of metabolites using the hypergeometric test.
#' Filters pathways by size and returns a data frame of enrichment statistics for each pathway.
#'
#' @param pathways A named list where each element is a character vector of metabolites in a pathway.
#' @param metabolites Character vector of metabolites of interest (e.g., DE peaks).
#' @param universe Character vector of all metabolites considered as background/universe.
#' @param minSize Minimum pathway size to include (default: 1).
#' @param maxSize Maximum pathway size to include (default: length(universe) - 1).
#' @param direction Character string or vector, e.g. 'up' or 'down', indicating direction of regulation (added to output).
#' @keywords internal
#' @return Data frame with columns: pathway, total, expected, hits, Raw.p, Holm.p, FDR, metabolites, direction.
#' @details Used internally by the ORA Shiny panel. Pathways with no overlap or outside size limits are excluded.
#' @keywords internal
bulk_utils_run_ORA = function (pathways, metabolites, universe, minSize = 1, maxSize = length(universe) - 1, direction = c("up", "down")) {

  if (!is.list(pathways)) {
    stop("pathways should be a list with each element containing metabolites from the universe")
  }

  if (any(duplicated(universe))) {
    warning("There were duplicate metabolites in universe, they were collapsed")
    universe = unique(universe)
  }

  empty_ora_result = data.frame(
    "pathway" = character(),
    "total" = numeric(),
    "expected" = numeric(),
    "hits" = numeric(),
    "Raw.p" = numeric(),
    "Holm.p" = numeric(),
    "FDR" = numeric(),
    "metabolites" = character(),
    "direction" = character()
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

  metabolitesFiltered = unique(stats::na.omit(metabolites[metabolites%in%universe]))
  if (length(metabolitesFiltered) == 0) {
    warning("No metabolites from the input list belong to the universe")
    return(empty_ora_result)
  }

  overlaps = lapply(pathwaysFiltered, intersect, metabolitesFiltered)
  overlap_metabolites = lapply(overlaps, function(x) paste0(x, collapse = '; '))
  overlapsT = data.frame(
    "q" = sapply(overlaps, length),
    "m" = sapply(pathwaysFiltered, length), # set.num
    "n" = length(universe) - sapply(pathwaysFiltered, length),
    "k" = length(metabolitesFiltered)) # q.size
  pathways_pvals = with(overlapsT, stats::phyper(q - 1, m, n, k, lower.tail = FALSE))
  expected = with(overlapsT, k * (m/length(universe)))
  res = data.frame(
    "pathway" = names(pathwaysFiltered),
    "total" = overlapsT$m,
    "expected" = expected,
    "hits" = overlapsT$q,
    "Raw.p" = pathways_pvals,
    "Holm.p" = stats::p.adjust(pathways_pvals, method = "holm"),
    "FDR" = stats::p.adjust(pathways_pvals, method = "BH") ,
    "metabolites" = unlist(overlap_metabolites,use.names = F),
    "direction" = direction
  )
  res = res[order(res$Raw.p), ]
  return(res)
}



#' Execute Over-Representation Analysis (ORA) for Up/Down Peaks
#'
#' @description Runs ORA for up- and down-regulated peaks, using user-supplied pathway definitions.
#'
#' @param de_peaks Data frame of differential peaks, must include columns 'lfc' and 'metabolite_id'.
#' @param anno Data frame with annotation, must include 'metabolite_id'.
#' @param pathway_db Data frame of pathway definitions with columns PathwayID, PathwayName, and MetaboliteIDs.
#' @param min_path_size Minimum pathway size to include in ORA.
#' @param min_pathway_hits Minimum number of hits in a pathway to include in results.
#' @keywords internal
#' @return Data frame of ORA results for up- and down-regulated peaks, or NULL if no peaks are present.
#' @details Used internally by the ORA Shiny panel. Pathways with <2 or >500 entries are excluded.
bulk_utils_execute_ora <- function(
  de_peaks,
  anno,
  pathway_db,
  min_path_size = 3,
  min_pathway_hits = 2
) {

  if (!('metabolite_id' %in% colnames(de_peaks))) {
    stop("de_peaks must include a 'metabolite_id' column.")
  }
  if (!('metabolite_id' %in% colnames(anno))) {
    stop("anno must include a 'metabolite_id' column.")
  }
  required_path_cols <- c('PathwayID', 'PathwayName', 'MetaboliteIDs')
  missing_path_cols <- setdiff(required_path_cols, colnames(pathway_db))
  if (length(missing_path_cols) > 0) {
    stop("Pathway table is missing required column(s): ", paste(missing_path_cols, collapse = ', '))
  }

  parse_ids <- function(x) {
    vals <- unlist(strsplit(as.character(x), "[,;|]"))
    vals <- trimws(vals)
    vals[vals != "" & !is.na(vals)]
  }

  up_peaks = unique(parse_ids(paste(de_peaks[de_peaks$lfc > 0, ]$metabolite_id, collapse = ',')))
  down_peaks = unique(parse_ids(paste(de_peaks[de_peaks$lfc < 0, ]$metabolite_id, collapse = ',')))

  # if no peaks are up AND down regulated, return NULL and display message
  if (length(c(up_peaks,down_peaks)) == 0) {
    return(NULL)
  }

  # Prepare pathway -> metabolite ID dictionary from user-supplied pathway table.
  pathway_names <- ifelse(is.na(pathway_db$PathwayName) | pathway_db$PathwayName == "",
                          as.character(pathway_db$PathwayID),
                          as.character(pathway_db$PathwayName))
  pathway2peaks <- split(pathway_db$MetaboliteIDs, pathway_names)
  pathway2peaks <- lapply(pathway2peaks, function(vals) unique(parse_ids(paste(vals, collapse = ','))))

  # keep only pathways with < 500 and > 1 entries
  pathway2peaks = pathway2peaks[which(lapply(pathway2peaks, length) < 500 & lapply(pathway2peaks, length) > 1)]

  universe_ids <- unique(parse_ids(paste(anno$metabolite_id, collapse = ',')))

  ora_up = bulk_utils_run_ORA(
    pathways = pathway2peaks,
    metabolites = up_peaks,
    universe = universe_ids,
    minSize = min_path_size,
    maxSize = 500,
    direction = 'up'
  ) |> as.data.frame()

  ora_down = bulk_utils_run_ORA(
    pathways = pathway2peaks,
    metabolites = down_peaks,
    universe = universe_ids,
    minSize = min_path_size,
    maxSize = 500,
    direction = 'down'
  ) |> as.data.frame()

  ora_combined = rbind(ora_up,ora_down)

  if (min_pathway_hits) {
    ora_combined = ora_combined[ora_combined$hits >= min_pathway_hits,]
  }
  return(ora_combined)
}


#' Create a Volcano Plot for ORA Results
#'
#' @description Generates a volcano plot for pathway enrichment results from ORA, showing log2 fold change vs. -log10(FDR).
#'
#' @param ORA_results Data frame of ORA results (from bulk_utils_execute_ora or bulk_utils_run_ORA).
#' @param pval.threshold FDR threshold for significance (default: 0.05).
#' @param selectedPathways Character vector of pathway names to label in the plot.
#'
#' @return List with elements: 'volcano' (ggplot object) and 'data' (data frame used for plotting).
#' @details Used internally by the ORA Shiny panel. Significant pathways are colored red, others gray.
#' @keywords internal
bulk_utils_ora_volcano_plot <- function(
  ORA_results,
  pval.threshold = 0.05,
  selectedPathways
){
  df = ORA_results |>
    dplyr::mutate("log10pval" = log10(.data$FDR),
                  "lfc" = log2(.data$hits/.data$expected)) |>
    dplyr::filter(!is.na(.data$log10pval))
  df$lfc = ifelse(df$direction=='up',df$lfc,-df$lfc)
  df$significance <- ifelse(df$FDR<pval.threshold,'Significant','Non-significant')
  df.label = df[df$pathway %in% selectedPathways,]
  lfc <- NULL; log10pval <- NULL; significance <- NULL
  vp <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = .data$lfc, y = -.data$log10pval,color=.data$significance,label = .data$pathway)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::xlab("log2(FC)") +
    ggplot2::ylab("-log10(pval)") +
    ggplot2::scale_color_manual(values=c("Non-significant"="#999999", "Significant"="#FF0000"))+
    ggplot2::geom_text(data = df.label, mapping = ggplot2::aes(x = .data$lfc, y = -.data$log10pval,label = .data$pathway))
  df[,'-log10pval']=-df$log10pval

  return(list('volcano'=vp,'data'=df))
}

##' Plot PCA Results
##'
##' @description Creates a PCA plot for samples, colored by a metadata column, with optional confidence ellipses.
##'
##' @param pca.res Result of prcomp or similar PCA function.
##' @param metadata Data frame of sample metadata.
##' @param annotation.id Integer index of metadata column to use for coloring.
##' @param show.confidence.ellipses Logical; whether to show confidence ellipses (default: TRUE).
##' @param label.force Passed to ggrepel::geom_label_repel (not used here, for future extension).
##'
##' @return List with element 'plot' (ggplot object).
##' @keywords internal
bulk_utils_plot_pca <- function(
  pca.res,
  metadata,
  annotation.id,
  show.confidence.ellipses = TRUE,
  label.force = 1
){
  annotation.name <- colnames(metadata)[annotation.id]

  expr.PCA <- dplyr::mutate(
                  as.data.frame(pca.res$x),
                  "name" = factor(metadata[, 1], levels = metadata[, 1]),
                  "condition" = if(!is.factor(metadata[,annotation.id])){
                    factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
                    } else {metadata[,annotation.id]}

  )
  pca.plot <- ggplot2::ggplot(expr.PCA, ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$condition, label = .data$name)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("PC1 (proportion of variance = ", summary(pca.res)$importance[2, 1] * 100, "%)"),
         y = paste0("PC2 (proportion of variance = ", summary(pca.res)$importance[2, 2] * 100, "%)"),
         colour = annotation.name)
  if(show.confidence.ellipses){
    pca.plot <- pca.plot +
      ggplot2::stat_ellipse(geom='polygon',alpha=0.3,ggplot2::aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)

  }
  pca.plot <- pca.plot + ggplot2::geom_point()
  list('plot'=pca.plot)
}

##' Plot PCA Contributions
##'
##' @description Plots the top contributing features to a selected PCA component.
##'
##' @param pca.res Result of prcomp or similar PCA function.
##' @param comp Integer; which component to plot (default: 1).
##' @param anno Data frame with annotation, must include 'display_name' and 'm_z'.
##'
##' @return List with elements: 'plot' (ggplot object) and 'table' (data frame of top contributors).
##' @keywords internal
bulk_utils_pca_contrib <- function(pca.res,
                        comp=1,
                        anno){
  contrib = as.data.frame(pca.res$rotation)
  contrib$metab = rownames(contrib)
  colnames(contrib)[comp]='PCAComp'
  contrib = contrib[order(-abs(contrib$PCAComp)),]
  contrib = data.frame(contrib)
  contrib$display_metab = stringr::str_wrap(anno$display_name[match(contrib$metab,anno$m_z)],30)
  contrib$metab = factor(contrib$metab,levels=rev(contrib$metab))
  contrib.plot = ggplot2::ggplot(utils::head(contrib,30),ggplot2::aes(x=.data$PCAComp,y=.data$metab,fill=.data$PCAComp>0,label=.data$display_metab))+
    ggplot2::geom_bar(stat='identity') +
    ggplot2::theme_minimal()+
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab('Contribution')+
    ggplot2::ylab('')
  return(list('plot'=contrib.plot,'table'=utils::head(contrib,30)))
}

##' Perform PLS-DA
##'
##' @description Runs PLS-DA (Partial Least Squares Discriminant Analysis) on intensity matrix and metadata.
##'
##' @param intensity.matrix Numeric matrix of intensities (features × samples).
##' @param metadata Data frame of sample metadata.
##' @param separator.id Integer index of metadata column to use as group/class.
##'
##' @return PLS-DA result object from mixOmics::plsda.
##' @keywords internal
bulk_utils_perform_plsda <- function(intensity.matrix,
                          metadata,
                          separator.id
){
  return(mixOmics::plsda(t(intensity.matrix),as.vector(metadata[,separator.id]) , ncomp = 2))
}

##' Plot PLS-DA Results
##'
##' @description Creates a PLS-DA plot for samples, colored by a metadata column, with optional confidence ellipses.
##'
##' @param intensity.matrix Numeric matrix of intensities (features × samples).
##' @param metadata Data frame of sample metadata.
##' @param separator.id Integer index of metadata column to use as group/class.
##' @param annotation.id Integer index of metadata column to use for coloring.
##' @param show.confidence.ellipses Logical; whether to show confidence ellipses (default: TRUE).
##' @param label.force Passed to ggrepel::geom_label_repel (not used here, for future extension).
##'
##' @return ggplot object.
##' @keywords internal
bulk_utils_plot_plsda <- function(
  intensity.matrix,
  metadata,
  separator.id,
  annotation.id,
  show.confidence.ellipses = TRUE,
  label.force = 1
){
  annotation.name <- colnames(metadata)[annotation.id]
  my.plsda <- bulk_utils_perform_plsda(intensity.matrix,metadata,separator.id)
  coords = my.plsda$variates$X
  expr.plsda <- dplyr::mutate(
    as.data.frame(coords),
    "name" = factor(metadata[, 1], levels = metadata[, 1]),
    "condition" = if(!is.factor(metadata[,annotation.id])){
      factor(metadata[, annotation.id], levels = unique(metadata[, annotation.id]))
    } else {
      metadata[,annotation.id]
    }
  )
  plsda.plot <- ggplot2::ggplot(expr.plsda, ggplot2::aes(x = .data$comp1, y = .data$comp2, colour = .data$condition, label = .data$name)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = paste0("PLS-DA Comp1 (proportion of variance = ", round(my.plsda$prop_expl_var$X[1] * 100,digits = 1), "%)"),
         y = paste0("PLS-DA Comp2 (proportion of variance = ", round(my.plsda$prop_expl_var$X[2] * 100,digits=1), "%)"),
         colour = annotation.name)
  if(show.confidence.ellipses){
    plsda.plot <- plsda.plot +
      ggplot2::stat_ellipse(geom='polygon',alpha=0.3,ggplot2::aes(fill = .data$condition, colour = .data$condition), show.legend = FALSE)
  }
  plsda.plot <- plsda.plot + ggplot2::geom_point()
  plsda.plot
}


##' Plot PLS-DA Contributions
##'
##' @description Plots the top contributing features to a selected PLS-DA component.
##'
##' @param intensity.matrix Numeric matrix of intensities (features × samples).
##' @param metadata Data frame of sample metadata.
##' @param separator.id Integer index of metadata column to use as group/class.
##' @param comp Integer; which component to plot (default: 1).
##' @param anno Data frame with annotation, must include 'display_name' and 'm_z'.
##'
##' @return List with elements: 'plot' (ggplot object) and 'table' (data frame of top contributors).
##' @keywords internal
bulk_utils_plsda_contrib <- function(intensity.matrix,
                          metadata,
                          separator.id,
                          comp=1,
                          anno){
  my.plsda <- bulk_utils_perform_plsda(intensity.matrix,metadata,separator.id)
  contrib = data.frame(my.plsda$loadings$X)
  contrib$metab = rownames(contrib)
  colnames(contrib)[comp]='PLSDAComp'
  contrib = contrib[order(-abs(contrib$PLSDAComp)),]
  contrib = data.frame(contrib)
  contrib$display_metab = stringr::str_wrap(anno$display_name[match(contrib$metab,anno$m_z)],30)
  contrib$metab = factor(contrib$metab,levels=rev(contrib$metab))
  contrib.plot = ggplot2::ggplot(utils::head(contrib,30),ggplot2::aes(x=.data$PLSDAComp,y=.data$metab,fill=.data$PLSDAComp>0,label=.data$display_metab))+
    ggplot2::geom_bar(stat='identity') +
    ggplot2::theme_minimal()+
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab('Contribution')+
    ggplot2::ylab('')
  return(list('plot'=contrib.plot,'table'=utils::head(contrib,30)))
}

##' Plot Barplot of Peak Intensities
##'
##' @description Plots a barplot of peak intensities for each sample, grouped by condition.
##'
##' @param sub.intensity.matrix Numeric matrix or data frame of intensities (features × samples).
##' @param log.transformation Logical; whether to log2-transform intensities (default: TRUE).
##' @param condition.vector Vector of condition/group assignments for each sample.
##'
##' @return ggplot object.
##' @keywords internal
bulk_utils_peaks_barplot <- function(sub.intensity.matrix,
                          log.transformation = TRUE,
                          condition.vector){
  if (log.transformation){
    log.intensity.matrix <- data.frame(log2(as.matrix(sub.intensity.matrix) + 1))
  } else {
    log.intensity.matrix <- sub.intensity.matrix
  }
  log.intensity.matrix$peak <- stringr::word(rownames(log.intensity.matrix),sep='_',1,2)
  melted.intensity.matrix <- tidyr::pivot_longer(log.intensity.matrix,
                                                  cols = colnames(log.intensity.matrix)[1:(ncol(log.intensity.matrix)-1)])
  melted.intensity.matrix$condition <- rep(condition.vector,nrow(sub.intensity.matrix))
  p <- ggplot2::ggplot(melted.intensity.matrix, ggplot2::aes(x = .data$name,
                                            y = .data$value,
                                            fill = .data$condition)) +
    ggplot2::geom_bar(stat='identity',position='dodge') +
    ggplot2::theme_minimal() +
    ggplot2::ylab(ifelse(log.transformation,'log2 intensity','intensity')) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),legend.position = "bottom") +
    ggplot2::facet_wrap(~.data$peak,scales='free')
  p
}

##' Plot Boxplot of Peak Intensities by Metadata
##'
##' @description Plots a boxplot of peak intensities for each sample, grouped by a metadata column.
##'
##' @param sub.intensity.matrix Numeric matrix or data frame of intensities (features × samples).
##' @param metadata Data frame of sample metadata.
##' @param log.transformation Logical; whether to log2-transform intensities (default: TRUE).
##' @param metadata.column Integer index of metadata column to use for grouping (default: 2).
##'
##' @return List with elements: 'plot' (ggplot object) and 'table' (data frame used for plotting).
##' @keywords internal
bulk_utils_peaks_boxplot <- function(sub.intensity.matrix,
                          metadata,
                          log.transformation = TRUE,
                          metadata.column = 2){
  if (log.transformation){
    log.intensity.matrix <- data.frame(log2(as.matrix(sub.intensity.matrix) + 1))
  } else {
    log.intensity.matrix <- sub.intensity.matrix
  }
  log.intensity.matrix$peak <- stringr::word(rownames(log.intensity.matrix),sep='_',1,2)
  melted.intensity.matrix <- tidyr::pivot_longer(log.intensity.matrix,
                                                  cols = colnames(log.intensity.matrix)[1:(ncol(log.intensity.matrix)-1)])
  melted.intensity.matrix$metadata = rep(metadata[,metadata.column],nrow(sub.intensity.matrix))
  melted.intensity.matrix$sample = rep(metadata[,1],nrow(sub.intensity.matrix))
  p <- ggplot2::ggplot(melted.intensity.matrix, ggplot2::aes(fill = .data$metadata,
                                            x = .data$metadata,
                                            y = .data$value,label=.data$sample)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width=0.1,color='black',fill='black') +
    ggplot2::theme_minimal() +
    ggplot2::ylab(ifelse(log.transformation,'log2 intensity','intensity')) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
    ggplot2::facet_wrap(~.data$peak,scales = 'free')
  p
  return.list = list('plot'=p,'table'=melted.intensity.matrix)
}

##' Update selectize inputs for peaks with sensible defaults
##'
##' Updates a Shiny selectize input for peaks, preselecting a specified number of initial choices if available.
##'
##' @param session Shiny session object.
##' @param inputId Character; input id of the selectize control.
##' @param choices Character vector of choices.
##' @param selected_n Integer; how many initial choices to preselect (default 2).
##' @keywords internal
utils_update_peak_selectize <- function(session, inputId, choices, selected_n = 2){
  shiny::validate(shiny::need(length(choices) > 0, paste0('No choices available for ', inputId)))
  selected_vals <- utils::head(choices, selected_n)
  shiny::updateSelectizeInput(session, inputId, choices = choices, server = TRUE, selected = selected_vals)
}