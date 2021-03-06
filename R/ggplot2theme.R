#' Tufte-inspired ggplot2 theme
#'
#' @param base_size Default font size = 14
#' @param base_family Font family. Defaults to "serif"
#' @param base_line_size defaults to base_size/22
#' @param base_rect_size defaults to base_size/22
#' @export
#' @return nothing
#' @examples
#' theme_set(theme_min())
theme_min = function (base_size = 14, base_family = "serif", base_line_size = base_size/22, 
          base_rect_size = base_size/22) 
{
  theme_bw(base_size = base_size, base_family = base_family, 
           base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace% 
    theme(axis.ticks = element_blank(), legend.background = element_blank(), 
          legend.key = element_blank(), panel.background = element_blank(), 
          panel.border = element_blank(), strip.background = element_blank(), 
          plot.background = element_blank(), panel.grid = element_blank(), complete = TRUE)
}


