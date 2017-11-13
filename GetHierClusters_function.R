GetHierClusters <- function(dta, n_neigh = NULL, hierarchical_method = 'ward.D2',
                            coord_names) {
  
  d <- as.dist(fields::rdist.earth(dta[, coord_names, with = FALSE]))
  cluster_tree <- hclust(d, method = hierarchical_method)
  memb <- cutree(cluster_tree, k = n_neigh)  
  
  tab <- table(memb)
  print(paste(length(tab), 'clusters'))
  print(paste('The range of cluster size is:', min(tab), '-', max(tab)))
  print(paste(sum(tab == 1), 'clusters have only 1 observation'))
  
  return(memb)
}