#' Create a basic cell type heatmap across all neighborhoods (co-occurrence of cell types in neighborhoods)
#'
#' @param nbhds.obj Neighborhoods object, as returned by sg_to_nbhds
#'
#' @return
#' @export
#'
#' @examples
celltype_heatmap = function(nbhds.obj) {
  pheatmap::pheatmap(celltype_cor(nbhds.obj), cluster_cols=FALSE, cluster_rows=FALSE,
           col=viridis::viridis(100), breaks=seq(-1, 1, length=101))
}


#' Heatmap of neighborhood cluster centroids ("regions") and cell type composition
#'
#' @param centroids Centroids from cluster_nbhds
#'
#' @return
#' @export
#'
#' @examples
nbhd_celltype_heatmap = function(centroids) {
  pheatmap::pheatmap(centroids, cluster_cols=FALSE, cluster_rows=FALSE,
                     col=viridis::viridis(100), breaks=seq(min(centroids), max(centroids), length=101),angle_col = "45")
}

#' Title
#'
#' @param clust.out 
#' @param pat 
#'
#' @return
#' @export
#'
#' @examples
region_spatial_plot = function(clust.out,pat) {
  d = tibble::tibble(x=pat$x,y=pat$y,type=pat$marks,region=clust.out$clusters)
  p = d %>%
    ggplot2::ggplot(ggplot2::aes(x=x,y=y,colour=region)) +
    ggplot2::geom_point()
  print(p)
}

#' Title
#'
#' @param nbhds.obj 
#' @param clusters 
#' @param against 
#' @param cell_types 
#'
#' @return
#' @export
#'
#' @examples
pseudospace_plot = function(nbhds.obj, clusters, against = "Tregs",cell_types = c("Tregs","CD8+ T cells")) {
  if(!(against %in% cell_types)) {
    stop("arg against must be in cell_types!")
  }
  nbhds = nbhds.obj$scaled_nbhds
  
  p = nbhds %>%
    tibble::as_tibble() %>%
    dplyr::mutate(cluster = clusters) %>%
    dplyr::arrange(!!rlang::sym(against)) %>%
    dplyr::mutate(ord_nbhd = 1:dim(.)[1]) %>%
    tidyr::pivot_longer(cols = -c(ord_nbhd,cluster),names_to = "cell_type",values_to = "comp") %>%
    dplyr::filter(cell_type %in% cell_types) %>% {
      ggplot2::ggplot(data=.,ggplot2::aes(x=ord_nbhd,y = comp, colour = cell_type)) +
        ggplot2::geom_smooth() + 
        ggnewscale::new_scale("color") +
        ggplot2::geom_segment(ggplot2::aes(colour = cluster,y=-0.5,yend=-0.4,xend=ord_nbhd)) +
        ggplot2::scale_colour_brewer(palette = "Set1")
    }
  print(p)
}

#' Title
#'
#' @param nbhds.obj 
#' @param clusters 
#' @param columns 
#'
#' @return
#' @export
#'
#' @examples
celltype_pairsplot = function(nbhds.obj, clusters, cell_types = NULL) {
  nbhds = nbhds.obj$scaled_nbhds
  if(is.null(cell_types)) {
    cell_types = 1:ncol(nbhds)
  }
  nbhds.obj$scaled_nbhds %>%
    tibble::as_tibble() %>%
    dplyr::mutate(cluster = factor(clusters)) %>%
    GGally::ggpairs(ggplot2::aes(alpha=0.4, colour = cluster, shape = cluster),columns = cell_types)
}