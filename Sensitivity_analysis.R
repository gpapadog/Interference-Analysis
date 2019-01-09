setwd('~/Github/Interference-Analysis/')
load('~/Dropbox/DATAverse/analysis_dat.dat')
source('GetHierClusters_function.R')
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
n_neigh <- c(70, 50)
hierarchical_method <- c('ward.D2', 'complete')
method_name <- c('Ward 70', 'Complete 50')
coord_names <- c('Fac.Longitude', 'Fac.Latitude')
trt_name <- 'SnCR'
out_name <- 'mean4maxOzone'
ps_with_re <- TRUE
B <- 3
alpha_level <- 0.05
num_alphas <- 5

# Renaming the analysis data to subdta and excluding the NOx emissions column.
subdta <- data.table::copy(analysis_dat)
subdta[, totNOxemissions := NULL]

cov_names <- c('pctCapacity_byHI', 'logHeatInput', 'Phase2', 'mostlyGas',
               'small_nunits', 'med_nunits', 'mean4MaxTemp', 'PctUrban',
               'PctWhite', 'PctBlack', 'PctHisp', 'PctHighSchool', 
               'MedianHHInc', 'PctPoor', 'PctOccupied', 'PctMovedIn5',
               'MedianHValue', 'logPopPerSQM')

extra_alphas <- c(0.1, 0.4)

plot_data <- NULL

for (cc in 1 : length(n_neigh)) {
  
  set.seed(1234)
  
  dta <- MakeFinalDataset(dta = subdta, n_neigh = n_neigh[cc],
                          hierarchical_method = hierarchical_method[cc],
                          coord_names = coord_names, trt_name = trt_name,
                          subset_clusters = FALSE)
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
  
  alpha_range <- quantile(obs_alpha$V1, probs = c(0.2, 0.8))
  alpha <- seq(alpha_range[1], alpha_range[2], length.out = num_alphas)
  alpha <- sort(c(alpha, extra_alphas))
  alpha[1] <- ifelse(alpha[1] == 0, 0.001, alpha[1])
  
  yhat_group <- GroupIPW(dta = dta, cov_cols = cov_cols, phi_hat = phi_hat,
                         alpha = alpha, trt_col = trt_col, out_col = out_col)
  yhat_group <- yhat_group$yhat_group
  
  scores <- CalcScore(dta = dta, neigh_ind = NULL, phi_hat = phi_hat,
                      cov_cols = cov_cols, trt_name = 'Trt')
  ypop <- Ypop(ygroup = yhat_group, ps = 'estimated', scores = scores,
               dta = dta)
  
  yhat_pop <- ypop$ypop
  yhat_pop_var <- ypop$ypop_var
  
  ps_info_est <- list(glm_form = glm_form, ps_with_re = ps_with_re,
                      gamma_numer = phi_hat[[1]], use_control = TRUE)
  boots_est <- BootVar(dta = dta, B = B, alpha = alpha, ps = 'est',
                       cov_cols = cov_cols, ps_info_est = ps_info_est,
                       verbose = TRUE, trt_col = trt_col, out_col = out_col,
                       return_everything = FALSE)
  
  de <- DE(ypop = yhat_pop, ypop_var = yhat_pop_var, boots = boots_est,
           alpha = alpha)
  ie <- IE(ygroup = yhat_group[, 1, ], ps = 'estimated', scores = scores,
           boots = boots_est)
  
  
  a1 <- which(alpha == extra_alphas[1])
  a2 <- which(alpha == extra_alphas[2])
  
  all_quant <- rep(c('DE', 'IE1', 'IE2'), each = length(alpha))
  x <- data.frame(alpha = rep(alpha, 3), quant = all_quant,
                  method = method_name[cc],
                  est = NA, LB = NA, UB = NA)
  x[1 : length(alpha), 4 : 6] <- t(de[c(1, 8, 9), ])
  x[(length(alpha) + 1) : (2 * length(alpha)), 4 : 6] <- t(ie[c(1, 8, 9), a1, ])
  x[(2 * length(alpha) + 1) : (3 * length(alpha)), 4 : 6] <- t(ie[c(1, 8, 9), a2, ])
  
  plot_data <- rbind(plot_data, x)
}

# load('~/Documents/Research/Interference/Revisions/Application/results_sens_1.dat')
# plot_data <- res$plot_data

# ---------------------------------------- #
# ------------ PAPER PLOTTING ------------ #
# ---------------------------------------- #

g <- NULL
nrow <- 1
ncol <- 2

g[[1]] <- ggplot(aes(x = alpha, y = est, ymin = LB, ymax = UB),
                 data = subset(plot_data, quant == 'DE')) +
  geom_line() +
  facet_wrap(~ method, nrow = nrow, ncol = ncol, scales = 'free_x') +
  ylab(expression(DE(alpha))) +
  geom_ribbon(alpha=0.3) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = 12))

g[[2]] <- ggplot(aes(x = alpha, y = est, ymin = LB, ymax = UB),
                 data = subset(plot_data, quant == 'IE1')) +
  geom_line() +
  facet_wrap(~ method, nrow = nrow, ncol = ncol, scales = 'free_x') +
  ylab(expression(IE(0.1,alpha))) +
  geom_ribbon(alpha=0.3) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())

g[[3]] <- ggplot(aes(x = alpha, y = est, ymin = LB, ymax = UB),
                 data = subset(plot_data, quant == 'IE2')) +
  geom_line() +
  facet_wrap(~ method, nrow = nrow, ncol = ncol, scales = 'free_x') +
  ylab(expression(IE(0.4,alpha))) +
  geom_ribbon(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab(expression(alpha)) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())


grid.arrange(grobs = g, nrow = 3, ncol = 1, heights = c(1, 0.9, 1))

