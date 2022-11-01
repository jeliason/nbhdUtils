#' Read a df of pcors (as produced in the sensitivity analysis nb) and convert to a list of pcor matrices
#'
#' @param file filename to be read from NETWORK_PATH subfolder
#' @param NETWORK_PATH path where networks are stored
#' @param split_on variable to split the df on
#'
#' @return list of pcor matrices
#' @export
#'
#' @examples
read_pcors_df_to_list <- function(file,NETWORK_PATH,split_on = "Spot") {
  pcors <- readRDS(paste0(NETWORK_PATH,file))
  pcors %>%
    named_group_split(!!sym(split_on)) -> pcors
  
  lapply(pcors,\(x) {
    x %>%
      dplyr::pivot_wider(names_from = "type2", values_from = "value") -> pcor
    
    rnames <- pcor$type1
    pcor %>%
      dplyr::select(-type1) %>%
      as.matrix() -> pcor
    pcor <- pcor[,c(ncol(pcor),1:(ncol(pcor)-1))]
    rownames(pcor) <- rnames
    diag(pcor) <- 1
    pcor
  }) -> pcors
  pcors
}

#' Read a df of pat cross functions and convert to a list of spotwise matrices
#'
#' @param file filename to be read from NETWORK_PATH subfolder
#' @param NETWORK_PATH path where networks are stored
#' @param distance distance at which to evaluate pat cross function
#' @param value name of method used
#' @param split_on variable to split on
#'
#' @return list of "split_on"-wise matrices
#' @export
#'
#' @examples
read_pat_cross_df_to_list <- function(file,NETWORK_PATH,distance=20,value = "gx",split_on = "spot") {
  pat_cross <- readRDS(paste0(NETWORK_PATH,file))
  pat_cross %>%
    named_group_split(!!sym(split_on)) -> pat_cross
  
  lapply(pat_cross,\(x) {
    x %>%
      dplyr::filter(dist == distance) %>%
      dplyr::select(!!sym(value),type1,type2) %>%
      dplyr::pivot_wider(names_from = "type1", values_from = all_of(value)) -> mat
    
    rnames <- mat$type2
    mat %>%
      dplyr::select(-type2) %>%
      as.matrix() -> mat
    rownames(mat) <- rnames
    mat
  }) -> mats
  
  mats
}

#' Helper to read in pcors df (removes lower triangle)
#'
#' @param file 
#' @param NETWORK_PATH 
#' @param shave logical, remove lower triangle?
#' @param split_on variable to split on
#'
#' @return tibble of pcors
#' @export
#'
#' @examples
read_pcors_df <- function(file,NETWORK_PATH,shave=TRUE,split_on = "Spot") {
  if(shave) {
    pcors_list <- read_pcors_df_to_list(file,NETWORK_PATH,split_on = split_on)
    lapply(1:length(pcors_list),\(i) {
      pcor <- pcors_list[[i]]
      nm <- names(pcors_list)[[i]]
      df  = pcor %>%
        corrr::as_cordf() %>%
        corrr::shave() %>%
        corrr::stretch() %>%
        dplyr::rename(type1=x,type2=y,value=r)-> df
      df %>%
        dplyr::mutate(!!sym(split_on) := nm) -> df
      df
    }) %>%
      dplyr::bind_rows() -> pcors
  } else {
    pcors <- readRDS(paste0(NETWORK_PATH,file))
  }
  pcors
}