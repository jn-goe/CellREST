#' theme_minimal_plots
#'
#' @param text_size Text Size.
#' @param text_col Text Color.
#' @param ... Additional parameters passed to `ggplot2::theme()`.
#' @description
#' Adapts legend ticks length and plot title position and fixes text size, text color, and a transparent background complementing `ggplot2::theme_minimal()`.
#'
#' @export
theme_minimal_plots <- function(text_size = 6, text_col = "black", ...) {
  theme_minimal(base_size = text_size) +
    theme(axis.text = element_text(color = text_col, size = text_size),
          axis.title = element_text(color = text_col, size = text_size),
          legend.text = element_text(color = text_col, size = text_size),
          legend.title = element_text(color = text_col, size = text_size),
          rect = element_rect(fill = "transparent"),
          legend.ticks.length = unit(c(.02, .02), 'in'),
          plot.title = element_text(hjust = .5), ...
    )
}


#' guide_small_legend_cont
#'
#' @description
#' Uses `ggplot2::guides()`, i.e. size and position, for plotting continuous `fill` and `color` values.
#'
#' @export
guide_small_legend_cont <- function() {
  guides(fill=guide_colorbar(keywidth=0.1,
                             keyheight=0.5,
                             default.unit="inch"),
         color=guide_colorbar(keywidth=0.1,
                              keyheight=0.5,
                              default.unit="inch")
  )
}

#' guide_small_legend_disc
#'
#' @description
#' Adapts `ggplot2::guides()`, i.e. size and position, for plotting discrete `fill` and `color` values.
#'
#' @export
guide_small_legend_disc <- function(...) {
  guides(fill=guide_legend(keywidth=0.1,
                           keyheight=0.15,
                           default.unit="inch",
                           override.aes = list(size = 1.5)),
         color=guide_legend(keywidth=0.1,
                            keyheight=0.15,
                            default.unit="inch",
                            override.aes = list(size = 1.5)),
         size=guide_legend(keywidth=0.1,
                            keyheight=0.15,
                            default.unit="inch",
                            override.aes = list(size = 1.5)),
                            ...)
}

#' numerical_to_viridis
#'
#' @param x Numerical vector.
#' @returns A vector of colors from the color palette `viridis` of the same length of `x` mimicking numerical distances between entries of `x`.
#'
#' @export
numerical_to_viridis <- function(x) {
  x <- round((100-1)*((x - min(x))/max(x - min(x)))+1)
  col_vec <- viridis::viridis(100)[x]
  return(col_vec)
}
