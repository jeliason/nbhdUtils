celltype_heatmap = function(nbhds.obj) {
  pheatmap(celltype_cor(nbhds.obj), cluster_cols=FALSE, cluster_rows=FALSE,
           col=viridis::viridis(100), breaks=seq(-1, 1, length=101))
}


nbhd_celltype_heatmap = function(centroids) {
  pheatmap::pheatmap(centroids, cluster_cols=FALSE, cluster_rows=FALSE,
                     col=viridis::viridis(100), breaks=seq(min(centroids), max(centroids), length=101),angle_col = "0")
}
