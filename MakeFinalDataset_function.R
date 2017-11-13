MakeFinalDataset <- function(dta, n_neigh = NULL, hierarchical_method = 'ward.D2',
                             coord_names, trt_name, out_name) {
  
  dta <- copy(subdta)
  setnames(dta, trt_name, 'Trt')
  
  memb <- GetHierClusters(dta, hierarchical_method = hierarchical_method,
                          n_neigh = n_neigh, coord_names = coord_names)
  dta[, cluster := memb]
  
  # What are the observed percentage of treated?
  obs_alpha <- dta[, sum(Trt) / length(Trt), by = cluster]
  hist(obs_alpha$V1, breaks = 100)
  print(paste(sum(obs_alpha$V1 %in% c(0, 1)), 'all treated/control'))
  
  keep_clusters <- obs_alpha$cluster[obs_alpha$V1 > 0 & obs_alpha$V1 < 1]
  dta <- subset(dta, cluster %in% keep_clusters)
  dta[, neigh := as.numeric(as.factor(cluster))]
  
  return(list(data = dta, obs_alpha = obs_alpha))
}