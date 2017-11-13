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
n_neigh <- 50
hierarchical_method <- 'ward.D2'
coord_names <- c('Fac.Longitude', 'Fac.Latitude')
trt_name <- 'SnCR'
out_name <- 'mean4maxOzone'

dta <- MakeFinalDataset(dta = subdta, clustering = clustering,
                        n_neigh = n_neigh,
                        hierarchical_method = hierarchical_method,
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


# What counterfactual F alpha are we considering?
alpha_range <- quantile(obs_alpha$V1[! (obs_alpha$V1 %in% c(0, 1))], probs = c(0.2, 0.8))

nonzero_alpha <- obs_alpha$V1[obs_alpha$V1 != 0]
probs <- data.frame(alpha = unique(sort(nonzero_alpha)))

probs$n_clust <- sapply(probs$alpha, function(x) sum(nonzero_alpha == x))
probs$obs_prob <- probs$n_clust / sum(probs$n_clust)

probs$in_range <- ifelse(probs$alpha > alpha_range[1] & probs$alpha < alpha_range[2], 1, 0)

probs$obs_prob_in <- 0
wh <- which(probs$in_range == 1)
probs$obs_prob_in[wh] <- probs$obs_prob[wh]
probs$obs_prob_in <- probs$obs_prob_in / sum(probs$obs_prob_in)

probs$hypoth_prob <- 0
wh <- which(probs$in_range == 1 & probs$alpha > median(nonzero_alpha))
probs$hypoth_prob[wh] <- probs$obs_prob[wh] / sum(probs$obs_prob[wh])




# ---------------------------------------- #
# ------------ PAPER PLOTTING ------------ #
# ---------------------------------------- #

all_probs <- c(probs$obs_prob, probs$obs_prob_in, probs$hypoth_prob)
all_dist <- c('Observed', 'Observed-restricted', 'Counterfactual')
plot_Falphas <- data.frame(alpha = rep(probs$alpha, 3), prob = all_probs,
                           Distribution = rep(all_dist, each = nrow(probs)))
plot_Falphas$Distribution <- factor(plot_Falphas$Distribution,
                                    levels = all_dist)

ggplot(aes(x = alpha, y = prob), data = plot_Falphas) +
  geom_segment(aes(xend = alpha), yend = 0, size = 0.7) +
  geom_vline(xintercept = alpha_range, linetype = 2) +
  xlab(expression(alpha)) +
  ylab('Probability mass function') +
  facet_wrap(~ Distribution) + coord_flip() +
  theme(panel.spacing = unit(2, "lines"))


alpha <- probs$alpha[probs$in_range == 1]
p1 <- probs$obs_prob_in[probs$in_range == 1]
p2 <- probs$hypoth_prob[probs$in_range == 1]




# ----------- Calculating the IPW ----------- #

trt_col <- which(names(dta) == 'Trt')
out_col <- which(names(dta) == out_name)

yhat_group <- GroupIPW(dta = dta, cov_cols = cov_cols, phi_hat = phi_hat,
                       alpha = alpha, trt_col = trt_col,
                       out_col = out_col)$yhat_group

ypop_trueps <- YpopTruePS(yhat_group, alpha)
yhat_pop <- ypop_trueps$ypop

scores <- CalcScore(dta = dta, neigh_ind = NULL, phi_hat = phi_hat,
                    cov_cols = cov_cols, trt_name = 'Trt')

# --------- Indirect effect ----------- #

ie_var <- IEvar(ygroup = yhat_group[, 1, ], alpha = alpha, ps = 'estimated',
                scores = scores)

prob_diff <- p2 - p1
ie_Falpha <- c(est = sum(prob_diff * yhat_pop[1, ]),
               var = delta_method(ie_var, vec = prob_diff))
rep(as.numeric(ie_Falpha[1]), 3) + c(- 1, 0, 1) * 1.96 * sqrt(ie_Falpha[2])




