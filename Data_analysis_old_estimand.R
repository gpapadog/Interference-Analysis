setwd('~/Github/Interference-Analysis/')
load('~/Dropbox/DATAverse/subdta.dat')
source('GetHierClusters_function.R')
source('MakeFinalDataset_function.R')
source('Old_estimands_IPW_function.R')

library(Interference)
library(rgdal)
library(raster)
library(proj4)
library(lme4)
library(data.table)
library(gplots)
library(gridExtra)
library(ggplot2)

n_neigh <- 50
hierarchical_method <- 'ward.D2'
coord_names <- c('Fac.Longitude', 'Fac.Latitude')
trt_name <- 'SnCR'
out_name <- 'mean4maxOzone'

dta <- MakeFinalDataset(dta = subdta, hierarchical_method = hierarchical_method,
                        n_neigh = n_neigh, coord_names = coord_names,
                        trt_name = trt_name, out_name = out_name)
obs_alpha <- dta$obs_alpha
dta <- dta$data

# ----------- Fitting the propensity score model ------------ #s

dta[, MedianHHInc := scale(MedianHHInc)]
dta[, MedianHValue := scale(MedianHValue)]
dta[, logHeatInput := scale(logHeatInput)]
dta[, mean4MaxTemp := scale(mean4MaxTemp)]
dta[, logPopPerSQM := scale(logPopPerSQM)]

cov_names <- c('pctCapacity_byHI', 'logHeatInput', 'Phase2', 'mostlyGas',
               'small_nunits', 'med_nunits', 'mean4MaxTemp', 'PctUrban',
               'PctWhite', 'PctBlack', 'PctHisp', 'PctHighSchool', 
               'MedianHHInc', 'PctPoor', 'PctOccupied', 'PctMovedIn5',
               'MedianHValue', 'logPopPerSQM')
cov_cols <- which(names(dta) %in% cov_names)
cov_names <- names(dta)[cov_cols]

glm_form <- paste('Trt ~ (1 | neigh) +', paste(cov_names, collapse = ' + '))
glmod <- glmer(as.formula(glm_form), data = dta, family = 'binomial',
               control = glmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5)))

phi_hat <- list(coefs = summary(glmod)$coef[, 1],
                re_var = as.numeric(summary(glmod)$varcor))


# ----------- Calculating the IPW ----------- #

setnames(dta, 'Trt', 'A')
setnames(dta, 'mean4maxOzone', 'Y')
dta <- as.data.frame(dta)

trt_col <- which(names(dta) == 'Trt')
out_col <- which(names(dta) == out_name)
alpha_range <- quantile(obs_alpha$V1[! (obs_alpha$V1 %in% c(0, 1))], probs = c(0.2, 0.8))
alpha <- seq(alpha_range[1], alpha_range[2], length.out = 40)


n_neigh <- length(unique(dta$neigh))
neigh_ind <- sapply(1 : n_neigh, function(x) which(dta$neigh == x))
yhat_group <- Old_IPW(alpha, dta, neigh_ind, phi_hat, cov_cols)

scores <- CalcScore(dta = dta, neigh_ind = NULL, phi_hat = phi_hat,
                    cov_cols = cov_cols, trt_name = 'A')
ypop <- Ypop(ygroup = yhat_group, ps = 'estimated', scores = scores)

yhat_pop <- ypop$ypop
yhat_pop_var <- ypop$ypop_var


# --------- Direct effect ----------- #
de <- DE(ypop = yhat_pop, ypop_var = yhat_pop_var, alpha = alpha)

# --------- Indirect effect ----------- #
ie <- IE(ygroup = yhat_group[, 1, ], ps = 'estimated', scores = scores)





# ---------------------------------------- #
# ------------ PAPER PLOTTING ------------ #
# ---------------------------------------- #

de_plot <- data.frame(alpha = alpha, de = de[1, ], low = de[3, ],
                      high = de[4, ])
g_de <- ggplot() + geom_line(aes(alpha, de), data = de_plot) +
  geom_ribbon(data = de_plot, aes(x = alpha, ymin = low, ymax = high),
              alpha=0.3) +
  geom_abline(intercept = 0, slope = 0, linetype = 2) +
  ylab(expression(DE(alpha))) + xlab(expression(alpha))


a1 <- c(1, 10, 20, 30, 40)
g_ie <- NULL

for (ii in a1) {
  
  ie_plot <- data.frame(alpha = alpha, ie = ie[1, ii, ])
  ie_plot$low <- ie_plot$ie - 1.96 * sqrt(ie[2, ii, ])
  ie_plot$high <- ie_plot$ie + 1.96 * sqrt(ie[2, ii, ])
  
  g_ie1 <- ggplot() + geom_line(aes(alpha, ie), data = ie_plot) +
    geom_ribbon(data = ie_plot, aes(x = alpha, ymin = low, ymax = high),
                alpha = 0.3) +
    geom_abline(intercept = 0, slope = 0, linetype = 2) +
    ylab(paste0('IE(', round(alpha[ii], 3), ',', expression(alpha), ')')) +
    xlab(expression(alpha))
  
  g_ie[[length(g_ie) + 1]] <- g_ie1
} 


grid.arrange(grobs = append(list(g_de), g_ie[c(2, 4)]), nrow = 1)
