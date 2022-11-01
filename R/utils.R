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
named_group_split <- function(.tbl, ...,keep = FALSE) {
  grouped <- dplyr::group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!dplyr::group_keys(grouped), sep = " / "))
  
  grouped %>% 
    dplyr::group_split(., .keep = keep) %>% 
    rlang::set_names(names)
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

#' Saves object as an RDS file to RESULTS_PATH subfolder
#'
#' @param obj R object
#' @param filename 
#'
#' @return
#' @export
#'
#' @examples
saveObj <- function(obj,filename,overwrite=TRUE) {
  full = paste0(RESULTS_PATH,filename)
  if(!file.exists(full) | overwrite) {
    saveRDS(obj,full)
  }
}

#' Reads object from RESULTS_PATH or other subfolder
#'
#' @param filename 
#' @param root if not NULL, will read from RESULTS_PATH
#'
#' @return R object
#' @export
#'
#' @examples
readObj <- function(filename,root=NULL) {
  if(is.null(root)) {
    obj <- readRDS(paste0(RESULTS_PATH,filename))
  } else {
    obj <- readRDS(paste0(root,filename))
  }
}

#' Standardize vector to [0,1]
#'
#' @param x vector to standardize
#'
#' @return standardized vector
#' @export
#'
#' @examples
range01 <- function(x){(x-min(x))/(max(x)-min(x))}