#' Power plant clustering.
#' 
#' Performing clustering of the power plants and acquires proportion of treated
#' in each cluster.
#' 
#' @param dta The data including coordinates.
#' @param n_neigh The number of clusters we want.
#' @param hierarchical_method The name of the hierarchical method to be used.
#' Options correspond to options in hclust(). Defaults to 'ward.D2'.
#' @param The names of the columns including the coordinates. The names should
#' be for longitude first, and then for latitude.
#' @param trt_name The name of the binary treatment column.
#' @param subset_clusters Logical. Whether clusters with proportion of treated
#' equal to 0 or 1 should be dropped. Defaults to TRUE.
MakeFinalDataset <- function(dta, n_neigh = NULL,
                             hierarchical_method = 'ward.D2',
                             coord_names, trt_name, subset_clusters = TRUE) {
  
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
  if (subset_clusters) {
    dta <- subset(dta, cluster %in% keep_clusters)
  }
  dta[, neigh := as.numeric(as.factor(cluster))]
  
  return(list(data = dta, obs_alpha = obs_alpha))
}