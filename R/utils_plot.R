# Plot utility helpers for SMEW Shiny app

##' Generic downloadHandler creator for ggplot objects
##'
##' Creates a Shiny downloadHandler for exporting ggplot objects to file with dynamic filename, width, and height.
##'
##' @param plot_func Function with no arguments returning a ggplot object.
##' @param filename_func Function returning filename string.
##' @param width_func Function returning width (numeric).
##' @param height_func Function returning height (numeric).
##' @param units Character; units for ggsave (default 'in').
##' @param dpi Numeric; DPI for ggsave (default 300).
##' @return A Shiny downloadHandler object to assign to output[[id]].
##' @keywords shiny, plot, download
create_download_plot_handler <- function(plot_func,
                                         filename_func,
                                         width_func,
                                         height_func,
                                         units = 'in',
                                         dpi = 300){
  shiny::downloadHandler(
    filename = function(){ filename_func() },
    content = function(file){
      # if (!requireNamespace("showtext", quietly = TRUE)) {
      #   stop("The 'showtext' package is required for PDF export with full symbol support. Please install it via install.packages('showtext').")
      # }
      # showtext::showtext_auto(enable = TRUE)
      # font_family <- "sans"
      # # Try to register DejaVuSans if showtextdb::font_add is available
      # if (requireNamespace("showtextdb", quietly = TRUE) && "font_add" %in% getNamespaceExports("showtextdb")) {
      #   # Try to add DejaVuSans from system fonts if available
      #   try({
      #     showtextdb::font_add("DejaVuSans", regular = NULL)
      #     font_family <- "DejaVuSans"
      #   }, silent = TRUE)
      # }
#      p <- plot_func() + ggplot2::theme(text = ggplot2::element_text(family = font_family))
      suppressWarnings(ggplot2::ggsave(file, plot = plot_func(), width = width_func(), height = height_func(), units = units, dpi = dpi))
 #     showtext::showtext_auto(enable = FALSE)
    }
  )
}
