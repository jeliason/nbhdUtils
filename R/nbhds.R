#' Extract cell neighborhoods from statgraphs object
#'
#' @param sg statgraphs object containing the edges for each node
#' @param pat spatstat ppp object that sg was created from
#'
#' @return matrix of neighborhood counts of each cell type
#' @export
#'
#' @examples
#' X = c(1,1,5)
#' Y = c(2,4,3)
#' cell_types = c("Tregs","tumor cells","vasculature")
#' pat = create_ppp(X,Y,cell_types,keep_types=TRUE)
#' sg = spatgraphs::spatgraph(pat,type="SIG")
#' nbhds = sg_to_nbhds(sg,pat)
sg_to_nbhds = function(sg,pat) {
  edges = sg$edges
  marks = pat$marks
  cmarks = as.character(marks)
  nnodes = length(edges)
  nbhds = sapply(1:nnodes,function(i) {
    x = c(cmarks[i],cmarks[edges[[i]]])
    x = factor(x,levels=levels(marks))
    as.vector(table(x))
  })
  nbhds = t(nbhds)
  colnames(nbhds) = levels(marks)
  return(nbhds)
}

#' Scale neighborhood counts matrix
#'
#' @param nbhds neighborhoods matrix from sg_to_nbhds
#' @param scale string indicating type of scaling, can be any of 
#' "comp","global.comp","binary", or "standardize"
#' 1. comp (default) - divide by total number of cells, yields local composition of cell types
#' 2. global.comp - divide by maximum number of cells across all neighborhoods per cell type, normalizes
#' each cell type to its own maximum
#' 3. binary - 1 if cell type is present, 0 otherwise
#' 4. standardize - normalizes so that cell types have mean 0 and standard deviation of 1
#' @return list with original nbhds, scaled nbhds and type of scaling
#' @export
#'
scale_nbhds = function(nbhds,scale="comp") {
  nbhds.obj = list(nbhds=nbhds)
  if(scale == "comp") {
    sums = rowSums(nbhds)
    nbhds.obj$scaled_nbhds = sweep(nbhds,1,STATS=sums,FUN="/")
  } else if(scale == "global.comp") {
    nbhds.obj$scaled_nbhds = sweep(nbhds,2,STATS=apply(nbhds,2,max),FUN="/")
  } else if(scale == "binary") {
    nbhds.obj$scaled_nbhds = (nbhds > 0) + 0
  } else if(scale == "standardize") {
    nbhds.obj$scaled_nbhds = scale(nbhds,T,T)
  } else {
    stop("Not a recognized scaling!")
  }
  nbhds.obj$scaling = scale
  return(nbhds.obj)
}

cluster_nbhds = function(nbhds.obj, blusparam = bluster::NNGraphParam(), scaled = TRUE) {
  nbhds = if(scaled) nbhds.obj$scaled_nbhds else nbhds.obj$nbhds
  clusters = bluster::clusterRows(nbhds, BLUSPARAM = blusparam)
  nclust = length(unique(clusters))
  centroids = sapply(1:nclust, function(i) {
    colMeans(nbhds[clusters == i,])
  })
  colnames(centroids) <- paste0("Region ", as.character(1:nclust))
  centroids
}