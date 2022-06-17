#' Create point pattern with certain cell types
#'
#' @param X X coordinates of cells
#' @param Y Y coordinates of cells
#' @param cell_types factor of cell type
#' @param keep_types string vector of types to retain in ppp, or "all" if should keep all types
#' @param ranges ranges of X and Y
#'
#' @return spatstat ppp object 
#' @export
create_ppp = function(X,Y,cell_types,keep_types="all",ranges=NULL) {
  cell_types = as.factor(cell_types)
  if("all" %in% keep_types) {
    keep_types = levels(cell_types)
  } else if(!is.vector(keep_types)) {
    stop("keep_types must be a vector of strings of cell types to keep!")
  }
  if(is.null(ranges)) {
    Xmin = min(X) - 5
    Xmax = max(X) + 5
    Ymin = min(Y) - 5
    Ymax = max(Y) + 5
  } else {
    Xmin = ranges[[1]][1]
    Xmax = ranges[[1]][2]
    Ymin = ranges[[2]][1]
    Ymax = ranges[[2]][2]
  }
  
  df = data.frame(X=X,Y=Y,cell_types=cell_types)
  
  filt = df %>%
    dplyr::filter(cell_types %in% keep_types) %>%
    dplyr::mutate(cell_types = as.factor(cell_types))
  
  pat = spatstat.geom::ppp(filt$X,filt$Y,c(Xmin,Xmax),c(Ymin,Ymax))
  
  spatstat.geom::marks(pat) <- filt$cell_types
  
  return(pat)
}

#' Split groups into list with names of each group
#'
#' @param .tbl tibble to perform group_by
#' @param ... additional arguments to dplyr::group_by
#'
#' @return list of groups with names
#' @export
#'
#' @examples
named_group_split <- function(.tbl, keep = FALSE, ...) {
  grouped <- dplyr::group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!dplyr::group_keys(grouped), sep = " / "))
  
  grouped %>% 
    dplyr::group_split(., .keep = keep) %>% 
    rlang::set_names(names)
}

#' Check a BRMS model using the DHARMa simulated residuals
#'
#' @param model brms model
#' @param integer integer response? (TRUE/FALSE)
#' @param plot make plot?
#' @param ... further arguments for DHARMa::plotResiduals 
#'
#' @return a DHARMa object
#' @export
#'
#' @examples
check_brms <- function(model,             
                       integer = TRUE,   
                       plot = TRUE,       
                       ...                
) {
  
  mdata <- brms::standata(model)
  if (!"Y" %in% names(mdata))
    stop("Cannot extract the required information from this brms model")
  
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model, ndraws = 1000, re.form = NULL)),
    observedResponse = mdata$Y, 
    fittedPredictedResponse = apply(
      t(brms::posterior_epred(model, ndraws = 1000, re.form = NULL)),
      1,
      mean),
    integerResponse = integer)
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  
  invisible(dharma.obj)
}


#' Make geometric neighborhoods dataframe
#'
#' @param df cleaned data frame from Schurch et al
#' @param spot.labels spots for which to create neighborhoods
#' @param keep_types cell types to keep
#' @param radius radius of neighborhood
#'
#' @return dataframe of neighborhood composition
#' @export
#'
#' @examples
make_spot_nbhds <- function(df,spot.labels,keep_types,radius=50) {

  nbhds.list = sapply(spot.labels, function(spot) {
    filt.df = df %>%
      dplyr::filter(spots == spot) %>%
      dplyr::select(-`...1`)
    pat = create_ppp(filt.df$X,filt.df$Y,filt.df$ClusterName,
                     keep_types = keep_types)
    sg = spatgraphs::spatgraph(pat,type="geometric",par=radius)
    
    nbhds = sg_to_nbhds(sg,pat) %>% tibble::as_tibble() %>% dplyr::mutate(spot = spot) %>%
      dplyr::mutate(X = filt.df$X, Y = filt.df$Y)
    
    list(nbhds)
  })
  sapply(nbhds.list,function(nbhds) length(colnames(nbhds)))
  nbhds = dplyr::bind_rows(nbhds.list) %>% tibble::as_tibble()
  nbhds[is.na(nbhds)] <- 0
  names(nbhds) <- make.names(colnames(nbhds))
  nbhds
}

#' Fit or load model along with fitting time
#'
#' @param fit.model model fitting function
#' @param file filename
#' @param ... arguments to pass to fit.model
#'
#' @return fitted model, from disk if already exists, otherwise will fit
#' @export
#'
#' @examples
run_time_model <- function(fit.model,file,...) {
  # file = paste0(file,".rds")
  if(file.exists(file)) {
    cat("Model already fit, reading from disk.\n")
    res = readRDS(file)
    model = res$m
  } else {
    cat("Model fitting\n")
    start = Sys.time()
    model = fit.model(...)
    end = Sys.time()
    total_time = end - start
    cat("Total time fitting: ", total_time,"\n")
    res = list(m=model,time=total_time)
    saveRDS(res,file)
  }
  model
}

#' Estimate kernel density in multitype point pattern
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
  
  split.pat = spatstat.geom::split.ppp(pat)
  dens.split = lapply(split.pat,function(pp) {
    tryCatch({
      spatstat.core::density.ppp(pp,sigma=sigma,diggle=T,eps=eps)
      # density(pp,sigma=bw.scott,eps=eps,diggle=T)
    },
    error=function(e) {
      spatstat.core::density.ppp(pp,diggle=T,eps=eps)
    })
  })
  # cat(paste0(unlist(lapply(dens.split,function(dens) length(dens$xcol)*length(dens$yrow))),collapse = "\n"),file = paste0("dens_size_",spot,".txt"))
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