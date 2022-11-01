#' GGplot helper to color heatmap tiles by whether or not they are significant according to an adjusted p-value
#'
#' @param alpha Significance level
#'
#' @return
#' @export
#'
#' @examples
sig_stars <- function(alpha=0.05) {
  list(ggplot2::geom_point(ggplot2::aes(shape=ifelse(p.adj < alpha, "dot", "no_dot"))),
       ggplot2::scale_shape_manual(values=c(dot=8, no_dot=NA), guide="none"))
}

#' GGplot helper to angle x-axis labels 45 degrees
#'
#' @return
#' @export
#'
#' @examples
anglex <- function() {
  list(theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1)))
}

#' Forest plot for coefficients from tidy summary of model
#'
#' @param tb tidied summary of model
#' @param xintercept where to draw vertical reference line
#'
#' @return a ggplot object
#' @export
#'
#' @examples
forest_plot <- function(tb,xintercept=0) {
  tb %>%
    dplyr::filter(!str_detect(term,"Intercept")) %>%
    ggplot2::ggplot(ggplot2::aes(estimate,term)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(aes(xmin=conf.low,xmax=conf.high)) +
    ggplot2::geom_vline(aes(xintercept=xintercept),color="red",linetype="dashed")-> p
  p
}