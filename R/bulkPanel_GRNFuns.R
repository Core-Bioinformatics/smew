#' Helper functions for Metabolite Peak Regulatory Network (GRN) Analysis
#'
#' These functions support the GRN inference and visualization modules.
#'
#' @name GRNFuns
#' @keywords internal
NULL

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
#' @return Adjacency matrix with regulators in rows, targets in columns.
#'   Values represent edge weights from GENIE3 algorithm.
#' @details When sample count exceeds 1000, a random subset of 1000 samples is used
#' for computational efficiency. The random seed is set for reproducibility.
#' @export
infer_GRN <- function(intensity_matrix, metadata, anno, seed = 13,
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
#' @return A data frame with fields from, to and value, describing the edges
#' of the network
#' @export
#' @examples 
#' weightMat <- matrix(
#'   c(0.1, 0.4, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' get_link_list_rename(weightMat, 2)
get_link_list_rename <- function(weightMat, plotConnections){
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
#' @inheritParams get_link_list_rename
#' @param weightMatList a list of (weighted) adjacency matrices; 
#' each list element must be an adjacency matrix with regulators in rows,
#' targets in columns
#' @return A vector containing the names of the recurring regulators
#' @export
#' @examples 
#' weightMat1 <- matrix(
#'   c(0.1, 0.4, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' weightMat2 <- matrix(
#'   c(0.1, 0.2, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' find_regulators_with_recurring_edges(list(weightMat1, weightMat2), 2)
find_regulators_with_recurring_edges <- function(weightMatList, plotConnections){
  edges <- lapply(weightMatList, function(wm){
    get_link_list_rename(wm, plotConnections)
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
#' @inheritParams GRNpanel
#' @inheritParams get_link_list_rename
#' @param anno a data frame with peak annotation including 'display_name' (peak names)
#' and 'm_z' (peak IDs); used to highlight target peaks in the network
#' @param plotConnections the number of connections to subset to for plotting
#' @param plot_position_grid,n_networks the position of the plot in 
#' the grid (1-4) and the number of networks shown (1-4); these are
#' solely used for hiding unwanted plots in the shiny app
#' @param recurring_regulators targets to be highlighted; usually the
#' result of \code{\link{find_regulators_with_recurring_edges}}
#' @return A network plot. See visNetwork package for more details.
#' @export
#' @examples 
#' weightMat1 <- matrix(
#'   c(0.1, 0.4, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' weightMat2 <- matrix(
#'   c(0.1, 0.2, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' anno <- tibble::tibble(ENSEMBL = c("r1", "r2", "t1", "t2"), NAME = ENSEMBL)
#' recurring_regulators <- find_regulators_with_recurring_edges(list(weightMat1, weightMat2), 2)
#' plot_GRN(weightMat1, anno, 2, 1, 1, recurring_regulators)
#' plot_GRN(weightMat2, anno, 2, 1, 1, recurring_regulators)
plot_GRN <- function(weightMat, anno, plotConnections, 
                     plot_position_grid, n_networks, recurring_regulators){
  
  if(n_networks >= plot_position_grid){
    edges <- get_link_list_rename(weightMat, plotConnections)
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
#' @inheritParams get_link_list_rename
#' @inheritParams find_regulators_with_recurring_edges
#' @return An UpSet plot. See UpSetR package for more details.
#' @export
#' @examples 
#' weightMat1 <- matrix(
#'   c(0.1, 0.4, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' weightMat2 <- matrix(
#'   c(0.1, 0.2, 0.8, 0.3), nrow = 2, ncol = 2,
#'   dimnames = list("regulators" = c("r1", "r2"), "targets" = c("t1", "t2"))
#' )
#' plot_upset(list(weightMat1, weightMat2), 2)
plot_upset <- function(weightMatList, plotConnections){
  if(length(weightMatList) > 1){
    edge.list <- lapply(weightMatList, function(wm){
      get_link_list_rename(wm, plotConnections) |>
        dplyr::mutate(connection = paste0(.data$from, "-", .data$to)) |>
        dplyr::pull(.data$connection)
    })
    names(edge.list) <- paste0("Network", seq_len(length(edge.list)))
    UpSetR::upset(UpSetR::fromList(edge.list), order.by = "degree", keep.order = TRUE)
  }else{
    NULL
  }
}