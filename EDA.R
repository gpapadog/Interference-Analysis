setwd('~/Documents/Interference/Application/SnCR_Gas_plants/')
load('~/Dropbox/DATAverse/subdta.dat')
source('~/Documents/Interference/Application/functions/GetClusters_function.R')
source('~/Documents/Interference/Application/functions/IndirectEffectPlot_function.R')

library(rgdal)
library(raster)
library(proj4)
library(Interference)
library(lme4)
library(data.table)
library(gplots)
library(gridExtra)
library(ggplot2)

msa_path <- "~/Dropbox/arepa_linkage/shapefile link/MSA/msa.shp"
cbsa_path <- paste0('~/Dropbox/arepa_linkage/shapefile link/',
                    'cb_2016_us_cbsa_500k/cb_2016_us_cbsa_500k.shp')

clustering <- 'hierarchical'
n_neigh <- 50
hierarchical_method <- 'ward.D2'
coord_names <- c('Fac.Longitude', 'Fac.Latitude')
trt_name <- 'SnCR'
out_name <- 'mean4maxOzone'

dta <- copy(subdta)
setnames(dta, trt_name, 'Trt')

memb <- GetClusters(dta, clustering = clustering, n_neigh = n_neigh,
                    hierarchical_method = hierarchical_method,
                    coord_names = coord_names, msa_path = msa_path,
                    cbsa_path = cbsa_path)
dta[, cluster := memb]


# What are the observed percentage of treated?
obs_alpha <- dta[, sum(Trt) / length(Trt), by = cluster]





library(grDevices)
library(plyr)
library(RColorBrewer)
library(ggmap)

us.dat <- map_data("state")

cols <- gray.colors(n_neigh, start = 0.1, end = 0.1, gamma = 2.2, alpha = NULL)

x <- dta[, list(Fac.Latitude, Fac.Longitude, cluster, Trt)]
x <- as.data.frame(x)
x$cluster <- as.factor(x$cluster)
x$Trt <- as.factor(x$Trt)

df <- x
find_hull <- function(x) {
  x[chull(x$Fac.Longitude, x$Fac.Latitude), ]
}
hulls <- ddply(df, "cluster", find_hull)

ggplot() + 
  geom_polygon(aes(long, lat, group = group), color = 'grey55', fill = 'grey90',
               data = us.dat) +
  geom_polygon(aes(Fac.Longitude, Fac.Latitude, fill = cluster), colour = 'grey55',
               data = hulls, alpha=.2) + 
  geom_point(data = x, aes(Fac.Longitude, Fac.Latitude, colour = cluster,
                           fill = cluster, shape = Trt), size = 2) +
  scale_color_manual('', values = cols) +
  guides(colour = 'none', fill = 'none') +
  scale_fill_manual("", values = cols) +
  scale_shape_manual(values=c(17, 10), labels = c("Other", "SCR/SNCR")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank()) +
  xlab('') + ylab('') +
  guides(shape = guide_legend(title = "")) +
  theme(legend.text = element_text(size = 13))










