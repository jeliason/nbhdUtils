#' Create point pattern with certain cell types
#'
#' @param X X coordinates of cells
#' @param Y Y coordinates of cells
#' @param cell_types vector of strings of cell type
#' @param keep_types string vector of types to retain in ppp, or TRUE if should keep all types
#'
#' @return spatstat ppp object 
#' @export
create_ppp = function(X,Y,cell_types,keep_types=TRUE) {
  if(keep_types == TRUE) {
    keep_types = unique(cell_types)
  } else if(!is.vector(keep_types)) {
    stop("keep_types must be a vector of strings of cell types to keep!")
  }
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