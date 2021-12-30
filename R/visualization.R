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
