#' Create kernel density from multitype point pattern
#'
#' @param pat 
#' @param eps 
#' @param sigma 
#'
#' @return
#' @export
#'
#' @examples
pat.density = function(pat,eps,sigma,...) {
  split.pat = spatstat.geom::split.ppp(pat)
  dens.split = lapply(split.pat,function(pp) {
    tryCatch({
      spatstat.core::density.ppp(pp,sigma=sigma,diggle=T,eps=eps,...)
    },
    error=function(e) {
      spatstat.core::density.ppp(pp,diggle=T,eps=eps,...)
    })
  })
  sp = dens.split[[2]]
  Y = rep(sp$yrow,times=sp$dim[2])
  X = rep(sp$xcol,each=sp$dim[1])
  step = sp$xstep
  dens = lapply(dens.split,function(d) {
    c(d$v)
  }) %>% 
    do.call(cbind,.) %>%
    tibble::as_tibble()
  
  list(dens=dens,X=X,Y=Y,step=sp$xstep)
}

#' Helper to estimate kernel density
#'
#' @param df_raw 
#' @param spot 
#' @param eps resolution of kernel density
#' @param sigma kernel type to use
#' @param ranges X and Y ranges to use for point pattern
#'
#' @return
#' @export
#'
#' @examples
kernel.density = function(df_raw,spot,eps,sigma=NULL,ranges=NULL) {
  df = df_raw %>%
    dplyr::filter(spots == spot)
  pat = create_ppp(df$X,df$Y,cell_types = df$type, ranges=ranges)
  
  pat.density(pat,eps,sigma)
}

#' Calculate cell-type (weighted) densities/frequencies around each cell.
#'
#' @param df_raw
#' @param spot 
#' @param add_markers 
#' @param nbhd_def 
#' @param nbhd_pars 
#' @param density 
#' @param weighted 
#'
#' @return
#' @export
#'
#' @examples
cell.densities <- function(df_raw,
                           spot,
                           add_markers = FALSE,
                           nbhd_def=spdep::dnearneigh,
                           nbhd_pars=list(0,50),
                           density=TRUE,
                           weighted=FALSE) {
  df = df_raw %>%
    dplyr::filter(spots == spot)
  markers = NULL
  if(add_markers) {
    markers = get.markers(df_raw,spot) %>% dplyr::select(-CellID)
  }
  d2 <- .Cell.densities(df$X,df$Y,df$type,CellID=df$CellID,
                        type.levels = levels(df$type),
                        markers = markers,
                        nbhd_def = nbhd_def,
                        nbhd_pars = nbhd_pars,
                        density=density,
                        weighted=weighted)
  d2
}

#' Calculate cell-type (weighted) densities/frequencies around each cell.
#'
#' @param X X coordinates for each cell
#' @param Y Y coordinates for each cell
#' @param types cell type for each cell
#' @param CellID IDs for each cell
#' @param type.levels possible cell types
#' @param markers cellular biomarkers to include with densities
#' @param nbhd_def neighborhood definition, returns list of neighbors for each cell
#' @param nbhd_pars parameters for neighborhood definition
#' @param density whether to convert cell type frequencies to densities
#' @param weighted whether to weight cells in neighborhood by distance from focal cell
#' @param weighting_scheme see "style" argument for spdep::nb2listwdist
#'
#' @return
#' @export
#'
#' @examples
.Cell.densities <- function(X,Y,
                            types,
                            CellID = 1:length(X),
                            type.levels=levels(types),
                            markers=NULL,
                            nbhd_def=spdep::dnearneigh,
                            nbhd_pars=list(0,50),
                            density=TRUE,
                            weighted=FALSE,
                            weighting_scheme = "raw") {
  nb <- do.call(nbhd_def,c(list(cbind(X,Y)),nbhd_pars))
  
  nb.ar = 1
  if(density) {
    nb.ar = nbhd_pars[[2]]^2*pi
  }
  if(weighted) {
    listw=spdep::nb2listwdist(nb,sp::SpatialPoints(cbind(X,Y)),zero.policy=TRUE,alpha=1,style=weighting_scheme,type="idw")
  }
  nbhds = sapply(1:length(nb),function(i) {
    x = c(types[i],types[nb[[i]]])
    x = factor(x,levels=type.levels)
    if(length(x) == 0) {
      return(rep(0,length(type.levels)))
    }
    if(weighted) {
      wt=c(1,listw$weights[[i]])
      aggregate(wt ~ x, FUN = sum,drop=F)$wt
    } else {
      as.vector(table(x))
    }
  })
  nbhds = t(nbhds)
  colnames(nbhds) = type.levels
  
  dens = (nbhds / nb.ar) %>% tibble::as_tibble()
  if(!is.null(markers)) {
    d2 <- dplyr::bind_cols(CellID=CellID,dens,markers,X=X,Y=Y)
  } else {
    d2 <- dplyr::bind_cols(CellID=CellID,dens,X=X,Y=Y)
  }
  d2 = d2 %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0)))
  d2
}
