MakeFinalDataset <- function(dta, clustering = c('hierarchical', 'MSA', 'CBSA'),
                             n_neigh = NULL, hierarchical_method = 'ward.D2',
                             coord_names, trt_name, out_name,
                             msa_path = NULL, cbsa_path = NULL) {
  
  clustering <- match.arg(clustering)
  if (clustering == 'MSA' & is.null(msa_path)) {
    msa_path <- "~/Dropbox/arepa_linkage/shapefile link/MSA/msa.shp"
  }
  if (clustering == 'CBSA' & is.null(cbsa_path)) {
    cbsa_path <- paste0('~/Dropbox/arepa_linkage/shapefile link/',
                        'cb_2016_us_cbsa_500k/cb_2016_us_cbsa_500k.shp')
  }
  
  
  dta <- copy(subdta)
  setnames(dta, trt_name, 'Trt')
  
  memb <- GetClusters(dta, clustering = clustering, n_neigh = n_neigh,
                      hierarchical_method = hierarchical_method,
                      coord_names = coord_names, msa_path = msa_path,
                      cbsa_path = cbsa_path)
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