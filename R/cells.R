celltype_cor = function(nbhds.obj) {
  nbhds = nbhds.obj$nbhds
  return(stats::cor(nbhds))
}