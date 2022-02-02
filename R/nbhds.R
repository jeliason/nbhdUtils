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
#' pat = create_ppp(X,Y,cell_types,keep_types="all")
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
    sums[sums == 0] = 1
    nbhds.obj$scaled_nbhds = sweep(nbhds,1,STATS=sums,FUN="/")
  } else if(scale == "global.comp") {
    maxs = apply(nbhds,2,max)
    maxs[maxs == 0] = 1
    nbhds.obj$scaled_nbhds = sweep(nbhds,2,STATS=maxs,FUN="/")
  } else if(scale == "binary") {
    nbhds.obj$scaled_nbhds = (nbhds > 0) + 0
  } else if(scale == "standardize") {
    nbhds.obj$scaled_nbhds = scale(nbhds,T,T)
  } else if(scale == "none") {
    nbhds.obj$scaled_nbhds = nbhds
  } else {
    stop("Not a recognized scaling!")
  }
  nbhds.obj$scaling = scale
  return(nbhds.obj)
}

#' Title
#'
#' @param pat 
#' @param point 
#' @param radius 
#'
#' @return
#' @export
#'
#' @examples
get_disc_nbhd = function(pat,point,radius=50) {
  marks = pat$marks
  bound <- spatstat.geom::disc(centre=point,radius=radius)
  isin<-spatstat.geom::inside.owin(x=pat$x,y=pat$y,w=bound)
  x <- marks[isin]
  x = factor(x,levels=levels(marks))
  x = as.vector(table(x))
  names(x) = levels(marks)
  x
}

#' Title
#'
#' @param pat 
#' @param radius 
#'
#' @return
#' @export
#'
#' @examples
get_raster_nbhds = function(pat,radius=50){
  frame <- spatstat.geom::Window(pat)
  
  xrange = frame$xrange[2]-frame$xrange[1]
  xstart = frame$xrange[1] + radius/2
  xstop = frame$xrange[2] - radius/2
  
  yrange = frame$yrange[2]-frame$yrange[1]
  ystart = frame$yrange[1] + radius/2
  ystop = frame$yrange[2] - radius/2
  
  nx=floor(xrange / radius*2)
  ny=floor(yrange / radius*2)
  
  xcenters=seq(xstart,xstop,length.out=nx)
  ycenters=seq(ystart,ystop,length.out=ny)
  
  g = pracma::meshgrid(xcenters,ycenters)
  names(g) <- c("x","y")
  
  nnbhds = length(g$x)
  raster_nbhds = sapply(1:nnbhds, function(i) {
    list(get_disc_nbhd(pat,point=c(g$x[i],g$y[i]),radius=radius))
  })
  raster_nbhds = do.call(rbind,raster_nbhds)
  raster_nbhds
}

#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' Cluster neighborhoods
#'
#' @param pat 
#' @param ntype 
#' @param scale 
#' @param blusparam 
#'
#' @return
#' @export
#'
#' @examples
cluster_nbhds = function(pat=NULL,nbhds.obj=NULL,ntype="gabriel",scale="comp", par=NULL, spat.dist=NULL,
                         blusparam = bluster::NNGraphParam()) {
  if(is.null(nbhds.obj) & !is.null(pat)) {
    if(ntype == "raster") {
      nbhds = get_raster_nbhds(pat,radius=par)
    } else {
      sg = spatgraphs::spatgraph(pat,type=ntype,par=par)
      nbhds = sg_to_nbhds(sg,pat)
    }
    nbhds.obj = scale_nbhds(nbhds,scale=scale)
  } else if(is.null(pat) & is.null(nbhds.obj)) {
    stop("Must pass either pat or nbhds.obj")
  }
  nbhds = nbhds.obj$scaled_nbhds
  if(!is.null(spat.dist)) {
    spat.dist[,1] = range01(spat.dist[,1])
    spat.dist[,2] = range01(spat.dist[,2])
  }
  nbhds = cbind(nbhds,spat.dist)
  clusters = bluster::clusterRows(nbhds, BLUSPARAM = blusparam)
  nclust = length(unique(clusters))
  centroids = sapply(1:nclust, function(i) {
    colMeans(nbhds[clusters == i,])
  })
  colnames(centroids) <- paste0("Region ", as.character(1:nclust))
  list(centroids=centroids,clusters=clusters)
}



#' Sweep cluster and neighborhood parameters and compare using ARI
#'
#' @param pat 
#' @param ntypes 
#' @param scales 
#' @param clust_params 
#'
#' @return
#' @export
#'
#' @examples
nbhdClusterSweep = function(pat,ntypes=list(list(type="SIG",pars=NULL),
                                            list(type="RNG",pars=NULL),
                                            list(type="gabriel",pars=NULL)),
                            scales=c("comp","global.comp","standardize"),
                            clust_params=list(list(BLUSPARAM=bluster::SNNGraphParam(),
                                                   optParams=list(k=c(5L,10L, 15L, 20L),
                                                                  cluster.fun=c("louvain", "infomap"))),
                                              list(BLUSPARAM=bluster::SomParam(centers=1L),
                                                   optParams=list(centers=5:10))
                            )) {
  graphs.list = sapply(ntypes,function(ntype) {
    type = ntype$type
    pars = ntype$pars
    if(is.vector(pars)) {
      sapply(pars,function(par) {
        list(spatgraphs::spatgraph(pat,type=type,par=par))
      })
    } else {
      list(spatgraphs::spatgraph(pat,type=type))
    }
  })
  graphs.list = unlist(graphs.list,recursive = F)
  grid = tidyr::expand_grid(graphs.list,scales,clust_params)
  
  out=sapply(1:dim(grid)[1],function(row) {
    params = grid[row,]
    sg = params[[1]][[1]]
    scale = params[[2]]
    cl_param = params[[3]][[1]]
    blusparam = cl_param$BLUSPARAM
    optParams = cl_param$optParams
    nbhds = sg_to_nbhds(sg,pat)
    nbhds.obj = scale_nbhds(nbhds,scale=scale)
    sweepOut = bluster::clusterSweep(nbhds.obj$scaled_nbhds, BLUSPARAM = blusparam,args=optParams)
    clus = sweepOut$clusters %>% as.matrix() %>% tibble::as_tibble()
    colnames(clus)<-paste(colnames(clus), paste0("scale.",scale,"_","sg.",sg$type,".par.",sg$parameters), sep = "_")
    list(clus)
  })
  
  df = do.call(cbind,out)
  df
  aris <- rdist::pdist(t(df),metric=mclust::adjustedRandIndex)
  g <- igraph::graph.adjacency(aris, mode="undirected", weighted=TRUE)
  meta <- igraph::cluster_walktrap(g)
  med_idx <- GDAtools::medoids(aris,meta$membership)
  combos = colnames(df)
  
  
  list(clusters=df,aris=aris,meta=meta,medoids=combos[med_idx])
}