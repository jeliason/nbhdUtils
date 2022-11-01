#' Calculate median nearest neighbor distance between all cell types in patient
#'
#' @param df table containing a column of patients and a column of spots
#' @param pats list of spatstat point patterns (as made by create_ppp)
#'
#' @return tibble of average distance between all cell types
#' @export
#'
#' @examples
# mean_nn_dist <- function(df,pats) {
#   patients <- unique(df$Patient)
#   lapply(1:length(patients),function(i) {
#     patient = patients[[i]]
#     
#     sp <- df %>% dplyr::filter(Patient == patient) %>% pull(Spot) %>% unique(.)
#     pat_list <- pats[which(names(pats) %in% sp)]
#     if(length(pat_list) == 0) {return()}
#     # get avg nn dist within each spot and bind together per patient
#     lapply(1:length(pat_list),\(j) {
#       pat = pat_list[[j]]
#       spot = names(pat_list)[[j]]
#       nnda = nndist(pat,by=marks(pat),k=1)
#       
#       mean_nn = aggregate(nnda, by=list(from=marks(pat)), median) %>%
#         pivot_longer(-from, names_to = "to",values_to = "average_distance")
#       
#       mean_nn$Spot = spot
#       
#       mean_nn
#     }) %>% bind_rows() -> mean_nn
#     
#     # average over nn dist within each patient, across spots
#     mean_nn %>%
#       mutate(average_distance = ifelse(is.infinite(average_distance),NA,average_distance)) %>%
#       group_by(from,to) %>%
#       summarise(average_distance = mean(average_distance,na.rm=T)) %>%
#       mutate(average_distance = ifelse(is.nan(average_distance),NA,average_distance)) %>%
#       mutate(Patient = patient) -> mean_nn
#     
#   }) %>%
#     bind_rows() %>%
#     ungroup() -> mean_nn
# }