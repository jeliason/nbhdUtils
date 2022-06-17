#' Partial Correlations
#' Taken from ppcor R package
#' @param x data matrix
#' @param method one of "pearson", "kendall", "spearman"
#'
#' @return list with element "estimate" containing the pcor matrix
#' @export
#'
#' @examples
pcor <- function(x, method = c("pearson", "kendall", "spearman"))
{
  # correlation method
  method <- match.arg(method)
  
  # check the data
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x)) 
    stop("supply a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  
  # sample number
  n <- dim(x)[1]
  
  # given variables' number
  gp <- dim(x)[2]-2
  
  # covariance matrix
  cvx <- cov(x,method=method)
  
  # inverse covariance matrix
  # if(det(cvx) < .Machine$double.eps){
  #   warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
  #   icvx <- ginv(cvx)
  # }else
  #   icvx <- solve(cvx)
  
  icvx <- MASS::ginv(cvx)
  # partial correlation
  pcor <- -cov2cor(icvx)
  diag(pcor) <- 1
  
  list(estimate=pcor)
}

#' Semipartial Correlations
#' Taken from ppcor R package
#' @param x data matrix
#' @param method one of "pearson", "kendall", "spearman"
#'
#' @return list with element "estimate" containing the spcor matrix
#' @export
#'
#' @examples
spcor <- function(x, method = c("pearson", "kendall", "spearman"))
{
  # correlation method
  method <- match.arg(method)
  
  # check the data
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x)) 
    stop("supply a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  
  # sample number
  # n <- dim(x)[1]
  # 
  # # given variables' number
  # gp <- dim(x)[2]-2
  
  # covariance matrix
  cvx <- cov(x,method=method)
  
  # # inverse covariance matrix
  # if(det(cvx) < .Machine$double.eps){
  #   warning("The inverse of variance-covariance matrix is calculated using Moore-Penrose generalized matrix invers due to its determinant of zero.")
  #   icvx <- ginv(cvx)
  # }else
  #   icvx <- solve(cvx)
  
  icvx <- MASS::ginv(cvx)
  
  # semi-partial correaltion
  spcor <- -cov2cor(icvx)/sqrt(diag(cvx))/sqrt(abs(diag(icvx)-t(t(icvx^2)/diag(icvx))))
  diag(spcor) <- 1
  
  list(estimate=spcor)
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
  
  dens = (nbhds / nb.ar) %>% tibble::as_tibble
  if(!is.null(markers)) {
    d2 <- dplyr::bind_cols(CellID=CellID,dens,markers,X=X,Y=Y)
  } else {
    d2 <- dplyr::bind_cols(CellID=CellID,dens,X=X,Y=Y)
  }
  d2 = d2 %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0)))
  d2
  # d2 <- d2 %>%
  #   tidyr::drop_na() %>%
  #   filter(!if_any(.fns=~is.infinite(.))) %>%
  #   select(where(~sum(!is.na(.x)) > 0))
  # # colnames(d2) <- make.names(names(d2))
  # d2
}

# note: this function does NOT include the focal cell of each neighborhood!
#' Title
#'
#' @param X 
#' @param Y 
#' @param types 
#' @param CellID 
#' @param type.levels 
#' @param markers 
#' @param nbhd_def 
#' @param nbhd_pars 
#' @param weighted 
#'
#' @return
#' @export
#'
#' @examples
weighted_cells <- function(X,Y,
                            types,
                            CellID = 1:length(X),
                            type.levels=levels(types),
                            markers=NULL,
                            nbhd_def=spdep::dnearneigh,
                            nbhd_pars=list(0,50),
                            weighted=TRUE) {
  nb <- do.call(nbhd_def,c(list(cbind(X,Y)),nbhd_pars))
  
  if(weighted) {
    listw=spdep::nb2listwdist(nb,sp::SpatialPoints(cbind(X,Y)),zero.policy = T,alpha=1,style="S",type="idw")
  }
  nbhds = sapply(1:length(nb),function(i) {
    x = types[nb[[i]]]
    x = factor(x,levels=type.levels)
    if(length(x) == 0) {
      return(rep(0,length(type.levels)))
    }
    if(weighted) {
      wt=listw$weights[[i]]
      tryCatch({
        aggregate(wt~x,FUN=sum,drop=F)$wt
      },
      error = function(e) {print(i); stop(e)})
      # dplyr::count(tibble(x=x,wt=wt),x,wt=wt,.drop=F) %>% dplyr::pull(n)
    } else {
      as.vector(table(x))
    }
  })
  nbhds = t(nbhds)
  colnames(nbhds) = type.levels
  
  # colnames(nbhds) <- make.names(colnames(nbhds))
  
  tibble::as_tibble(nbhds) %>% 
    dplyr::mutate(dplyr::across(dplyr::everything(), ~tidyr::replace_na(.x, 0))) %>%
    dplyr::mutate(CellID = CellID,
           X = X,
           Y = Y)
}

#' Title
#'
#' @param x 
#' @param y 
#' @param X 
#' @param weight_fun 
#' @param weight_pars 
#' @param pcor.fun 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
.Pcors <- function(x,y,
                   X, # vars to calculate partial correlation of
                   weight_fun=mgwrsar::bisq_C,
                   weight_pars=list(200,0),
                   pcor.fun = pcor,
                   method = "kendall") {
  dists = rdist::pdist(cbind(x,y))
  weights = sapply(1:length(x), function(i) {
    wi = do.call(weight_fun,c(list(dists[i,]),weight_pars))
  })
  
  ncolX = ncol(X)
  nrowX = nrow(X)
  m.names <- colnames(X)
  pcors = array(0,dim = c(ncolX,ncolX,nrowX),dimnames = list(m.names,m.names,NULL))
  for(i in 1:nrowX) {
    wt = weights[,i]
    X.wt = X[wt > 0,]
    wt = wt[wt > 0]
    ix = which(lapply(1:ncolX,function(j) {
      var(X.wt[,j])
    }) > .Machine$double.eps)
    X.wt = X.wt[,ix]
    tryCatch({
      # pcor = corpcor::pcor.shrink(X.wt[,ix],w=wt,verbose=F)
      X.wt = sweep(X.wt,1,sqrt(wt),FUN="*")
      pcor.mat = pcor.fun(X.wt,method=method)$estimate
    }, error = function(e) {
      print(i)
      stop(e)
    })
    pcors[ix,ix,i] = pcor.mat
  }
  pcors
}

#' Title
#'
#' @param df_raw 
#' @param spot 
#' @param add_markers 
#' @param densities 
#' @param nbhd_def 
#' @param nbhd_pars 
#' @param weight_fun 
#' @param weight_pars 
#' @param pcor.fun 
#'
#' @return
#' @export
#'
#' @examples
densities.pcors <- function(df_raw,
                            spot,
                            add_markers=FALSE,
                            densities = TRUE,
                            nbhd_def=spdep::dnearneigh,
                            nbhd_pars=list(0,50),
                            weight_fun=mgwrsar::bisq_C,
                            weight_pars=list(200,0),
                            pcor.fun = pcor) {
  df = df_raw %>%
    dplyr::filter(spots == spot)
  
  d2 = cell.densities(df_raw,spot,add_markers,nbhd_def,nbhd_pars)
  pcors = .Pcors(d2$X,d2$Y,d2 %>% select(-c(CellID,X,Y)),
                 weight_fun = weight_fun,
                 weight_pars = weight_pars,
                 pcor.fun = pcor.fun)
  
  list(d2=d2,pcors=pcors,clusters=df$type)
}

#' Title
#'
#' @param pcors 
#' @param CellID 
#' @param full 
#'
#' @return
#' @export
#'
#' @examples
pcors.table <- function(pcors,CellID = 1:dim(pcors)[3],full = TRUE) {
  nfeat = dim(pcors)[1]
  ncols = nfeat*(nfeat-1)
  if (!full) {
    ncols = ncols / 2
  }
  nobs = dim(pcors)[3]
  tab = matrix(0,nrow = dim(pcors)[3],ncol = ncols)
  
  mat1 = pcors[,,1]
  if(full) {
    i1 <- lower.tri(mat1) | upper.tri(mat1)
  } else {
    i1 = lower.tri(mat1)
  }
  i2 <- which(i1, arr.ind=TRUE)
  tnames <- data.frame(Var1 = colnames(mat1)[i2[,1]], 
                       Var2 = colnames(mat1)[i2[,2]]) %>% as_tibble %>% mutate(comb = paste0(Var1,":",Var2)) %>% pull(comb)
  colnames(tab) <- tnames
  for(row in 1:nobs) {
    tab[row,] = pcors[,,row][i1]
  }
  tab = as.data.frame(tab)
  tab$CellID = CellID
  tab
}

#' Title
#'
#' @param coords 
#' @param nbhd 
#' @param neighbors 
#' @param threshold 
#'
#' @return
#' @export
#'
#' @examples
spatial_contexts = function(coords,nbhd,neighbors=50,threshold=0.9) {
  nbs = spdep::knn2nb(spdep::knearneigh(coords,neighbors))
  sapply(1:nrow(coords),function(i) {
    nb = c(i,nbs[[i]])
    tb = sort(table(nbhd[nb]),decreasing = T)
    ix = which(cumsum(tb) > threshold*(length(nb)))
    res <- tryCatch({
      paste(sort(names(tb)[1:ix[1]]),collapse = ":")
    },
    error=function(e) {
      print(e)
      print(nbhd[nb])
      print(length(nb))
      print(cumsum(tb))
      print(cumsum(tb) > threshold*(length(nb)))
      stop(e)
    })
  })
}

# d is a vector with the off-diagonal elements of pcor matrix. This function converts that vector back to a square matrix.
#' Title
#'
#' @param d 
#' @param lower.only 
#'
#' @return
#' @export
#'
#' @examples
build_matrix_pcors = function(d,diag=1) {
  d %>%
    enframe(name = "cell.types") %>%
    unnest(value) %>% 
    filter(str_detect(cell.types,":")) %>%
    separate(cell.types,c("Type1","Type2"),sep = ":") %>%
    pivot_wider(names_from = Type2,values_from = value) %>%
    column_to_rownames("Type1") %>%
    as.matrix -> m
  
  n_col = length(rownames(m))+1
  mat = matrix(0,nrow = n_col,ncol=n_col)
  mat[lower.tri(mat,diag=F)] = m[lower.tri(m,diag=T)]
  mat[upper.tri(mat,diag=F)] = t(mat)[upper.tri(mat,diag=F)]
  rownames(mat) = c(colnames(m)[1],rownames(m))
  colnames(mat) = c(colnames(m),rownames(m)[length(rownames(m))])
  diag(mat) <- diag
  mat
}

#' Title
#'
#' @param m 
#'
#' @return
#' @export
#'
#' @examples
remove_zero_cols_rows = function(m) {
  ix = which(lapply(1:ncol(m),function(j) {
    var(m[,j])
  }) > .Machine$double.eps)
  m[ix,ix]
}

plot.pcors <- function(d2,types,density.type,squared = T) {
  types.fm = make.names(types)
  types.string = paste0(types.fm[1],":",types.fm[2])
  ix = which(colnames(d2) == types.string)
  if(length(ix) == 0) {
    types.string = paste0(types.fm[2],":",types.fm[1])
    ix = which(colnames(d2) == types.string)
  }
  if(squared) {
    p <- d2 %>%
      ggplot() +
      geom_point(aes(X,Y,size = !!sym(make.names(density.type)),colour=!!sym(types.string))) +
      scale_colour_viridis(option="plasma",direction=1,limit=c(0,max(abs(d2[,ix])))) +
      labs(colour = "") + 
      ggtitle(paste0("Percentage of variation in ",types[1]," accounted for by\n variation in ",types[2]))
  }
  else {
    p <- d2 %>%
      ggplot() +
      geom_point(aes(X,Y,size = !!sym(make.names(density.type)),colour=!!sym(types.string))) +
      scale_colour_distiller(palette="PiYG",direction=1,limit=c(-max(abs(d2[,ix])),max(abs(d2[,ix])))) +
      labs(colour = "") + 
      ggtitle(paste0("Unique correlation of ",types[2]," with ",types[1]))
  }
  p
}

#' Title
#'
#' @param df_raw 
#' @param spot 
#'
#' @return
#' @export
#'
#' @examples
get.markers <- function(df_raw, spot) {
  d = df_raw %>%
    filter(spots == spot) %>%
    select(c(`CD44 - stroma`:`MMP12 - matrix metalloproteinase`,CellID))
  d
}