#' Create point pattern with certain cell types
#'
#' @param X X coordinates of cells
#' @param Y Y coordinates of cells
#' @param cell_types string of cell type
#' @param keep_types types to retain in ppp
#'
#' @return spatstat ppp object 
#' @export
create_ppp = function(X,Y,cell_types,keep_types=c('Tregs','CD8+ T cells')) {
  Xmin = min(X) - 5
  Xmax = max(X) + 5
  Ymin = min(Y) - 5
  Ymax = max(Y) + 5
  
  df = data.frame(X=X,Y=Y,cell_types=cell_types)
  
  filt = df %>%
    dplyr::filter(cell_types %in% keep_types) %>%
    dplyr::mutate(cell_types = factor(cell_types))
  
  pat = spatstat.geom::ppp(filt$X,filt$Y,c(Xmin,Xmax),c(Ymin,Ymax))
  
  spatstat.geom::marks(pat) <- filt$cell_types
  
  return(pat)
}