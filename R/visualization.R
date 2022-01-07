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
    ggplot(aes(x=x,y=y,colour=region)) +
    geom_point()
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
    as_tibble() %>%
    mutate(cluster = clusters) %>%
    arrange(!!sym(against)) %>%
    mutate(ord_nbhd = 1:dim(.)[1]) %>%
    pivot_longer(cols = -c(ord_nbhd,cluster),names_to = "cell_type",values_to = "comp") %>%
    filter(cell_type %in% cell_types) %>% {
      ggplot(data=.,aes(x=ord_nbhd,y = comp, colour = cell_type)) +
        geom_smooth() + 
        new_scale("color") +
        geom_segment(aes(colour = cluster,y=-0.5,yend=-0.4,xend=ord_nbhd)) +
        scale_colour_brewer(palette = "Set1")
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
pairsplot = function(nbhds.obj, clusters, columns = NULL) {
  nbhds = nbhds.obj$scaled_nbhds
  if(is.null(columns)) {
    columns = 1:ncol(nbhds)
  }
  nbhds.obj$scaled_nbhds %>%
    as_tibble() %>%
    mutate(cluster = factor(clusters)) %>%
    ggpairs(aes(alpha=0.4, colour = cluster, shape = cluster),columns = columns)
}