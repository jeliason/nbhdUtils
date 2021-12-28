create_ppp = function(df,cell_types=c('Tregs','CD8+ T cells')) {
  Xmin = min(df$X) - 5
  Xmax = max(df$X) + 5
  Ymin = min(df$Y) - 5
  Ymax = max(df$Y) + 5
  
  filt = df %>%
    filter(ClusterName %in% cell_types) %>%
    mutate(ClusterName = factor(ClusterName))
  
  pat = ppp(filt$X,filt$Y,c(Xmin,Xmax),c(Ymin,Ymax))
  
  marks(pat) <- filt$ClusterName
  
  return(pat)
}