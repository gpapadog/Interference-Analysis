setwd('~/Documents/Interference/Application/SnCR_Gas_plants/')
load('~/Dropbox/DATAverse/subdta.dat')
source('~/Documents/Interference/Application/functions/GetClusters_function.R')
source('~/Documents/Interference/Application/functions/IndirectEffectPlot_function.R')
source('MakeFinalDataset_function.R')

library(rgdal)
library(raster)
library(proj4)
library(Interference)
library(lme4)
library(data.table)
library(gplots)
library(gridExtra)
library(ggplot2)


clustering <- 'hierarchical'
n_neigh <- c(30, 70, 50)
hierarchical_method <- c('ward.D2', 'ward.D2', 'complete')
method_name <- c('Ward 30', 'Ward 70', 'Complete 50')
coord_names <- c('Fac.Longitude', 'Fac.Latitude')
trt_name <- 'SnCR'
out_name <- 'mean4maxOzone'

cov_names <- c('pctCapacity_byHI', 'logHeatInput', 'Phase2', 'mostlyGas',
               'small_nunits', 'med_nunits', 'mean4MaxTemp', 'PctUrban',
               'PctWhite', 'PctBlack', 'PctHisp', 'PctHighSchool', 
               'MedianHHInc', 'PctPoor', 'PctOccupied', 'PctMovedIn5',
               'MedianHValue', 'logPopPerSQM')

extra_alphas <- c(0.226, 0.414)

plot_data <- NULL

for (cc in 1 : length(n_neigh)) {
  
  dta <- MakeFinalDataset(dta = subdta, clustering = clustering,
                          n_neigh = n_neigh[cc],
                          hierarchical_method = hierarchical_method[cc],
                          coord_names = coord_names, trt_name = trt_name,
                          out_name = out_name)
  obs_alpha <- dta$obs_alpha
  dta <- dta$data

  # ----------- Fitting the propensity score model ------------ #s
  
  dta[, MedianHHInc := scale(MedianHHInc)]
  dta[, MedianHValue := scale(MedianHValue)]
  dta[, logHeatInput := scale(logHeatInput)]
  dta[, mean4MaxTemp := scale(mean4MaxTemp)]
  dta[, logPopPerSQM := scale(logPopPerSQM)]
  
  cov_cols <- which(names(dta) %in% cov_names)
  cov_names <- names(dta)[cov_cols]
  
  glm_form <- paste('Trt ~ (1 | neigh) +', paste(cov_names, collapse = ' + '))
  glmod <- glmer(as.formula(glm_form), data = dta, family = 'binomial',
                 control = glmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 2e5)))
  
  phi_hat <- list(coefs = summary(glmod)$coef[, 1],
                  re_var = as.numeric(summary(glmod)$varcor))
  
  # ----------- Calculating the IPW ----------- #
  
  trt_col <- which(names(dta) == 'Trt')
  out_col <- which(names(dta) == out_name)
  
  alpha_range <- quantile(obs_alpha$V1[! (obs_alpha$V1 %in% c(0, 1))], probs = c(0.2, 0.8))
  alpha <- seq(alpha_range[1], alpha_range[2], length.out = 40)
  alpha <- sort(c(alpha, extra_alphas))
  
  yhat_group <- GroupIPW(dta = dta, cov_cols = cov_cols, phi_hat = phi_hat,
                         alpha = alpha, trt_col = trt_col, out_col = out_col)
  yhat_group <- yhat_group$yhat_group
  
  ypop_trueps <- YpopTruePS(yhat_group, alpha)
  yhat_pop <- ypop_trueps$ypop
  
  scores <- CalcScore(dta = dta, neigh_ind = NULL, phi_hat = phi_hat,
                      cov_cols = cov_cols, trt_name = 'Trt')
  yhat_pop_var <- VarEstPS(dta = dta, yhat_group = yhat_group,
                           yhat_pop = yhat_pop, neigh_ind = NULL,
                           phi_hat = phi_hat, cov_cols = cov_cols,
                           var_true = ypop_trueps$ypop_var, scores = scores)

  de <- DE(ypop = yhat_pop, ypop_var = yhat_pop_var, alpha = alpha)
  de <- rbind(de, low_int = de[1, ] - 1.96 * sqrt(de[2, ]))
  de <- rbind(de, high_int = de[1, ] + 1.96 * sqrt(de[2, ]))

  ie_var <- IEvar(ygroup = yhat_group[, 1, ], alpha = alpha, ps = 'estimated',
                  scores = scores)
  ie <- IE(ypop = yhat_pop[1, ], ypop_var = ie_var, alpha = alpha)
  
  
  a1 <- which(alpha == extra_alphas[1])
  a2 <- which(alpha == extra_alphas[2])
  
  all_quant <- rep(c('DE', 'IE1', 'IE2'), each = length(alpha))
  x <- data.frame(alpha = rep(alpha, 3), quant = all_quant,
                  method = method_name[cc],
                  est = NA, LB = NA, UB = NA)
  x[1 : length(alpha), 4 : 6] <- t(de[c(1, 3, 4), ])
  x[(length(alpha) + 1) : (2 * length(alpha)), 4 : 6] <- t(ie[c(1, 3, 4), a1, ])
  x[(2 * length(alpha) + 1) : (3 * length(alpha)), 4 : 6] <- t(ie[c(1, 3, 4), a2, ])
  
  plot_data <- rbind(plot_data, x)
}




# ---------------------------------------- #
# ------------ PAPER PLOTTING ------------ #
# ---------------------------------------- #

g <- NULL
nrow <- 1
ncol <- 3

g[[1]] <- ggplot(aes(x = alpha, y = est, ymin = LB, ymax = UB),
                 data = subset(plot_data, quant == 'DE')) +
  geom_line() +
  facet_wrap(~ method, nrow = nrow, ncol = ncol) +
  ylab(expression(DE(alpha))) +
  geom_ribbon(alpha=0.3) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

g[[2]] <- ggplot(aes(x = alpha, y = est, ymin = LB, ymax = UB),
                 data = subset(plot_data, quant == 'IE1')) +
  geom_line() +
  facet_wrap(~ method, nrow = nrow, ncol = ncol) +
  ylab(expression(IE(0.226,alpha))) +
  geom_ribbon(alpha=0.3) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

g[[3]] <- ggplot(aes(x = alpha, y = est, ymin = LB, ymax = UB),
                 data = subset(plot_data, quant == 'IE2')) +
  geom_line() +
  facet_wrap(~ method, nrow = nrow, ncol = ncol) +
  ylab(expression(IE(0.414,alpha))) +
  geom_ribbon(alpha=0.3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(expression(alpha))



grid.arrange(grobs = g, nrow = ncol, ncol = nrow)
