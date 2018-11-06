setwd('~/Github/Interference-Analysis/')
load('~/Dropbox/DATAverse/subdta.dat')
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
set.seed(1234)

n_neigh <- 50
hierarchical_method <- 'ward.D2'
coord_names <- c('Fac.Longitude', 'Fac.Latitude')
trt_name <- 'SnCR'
out_name <- 'mean4maxOzone'
ps_with_re <- TRUE
B <- 500
alpha_level <- 0.05
save_results <- TRUE
save_path <- '~/Documents/Research/Interference/Revisions/Application/'

dta <- MakeFinalDataset(dta = subdta, hierarchical_method = hierarchical_method,
                        n_neigh = n_neigh, coord_names = coord_names,
                        trt_name = trt_name, subset_clusters = FALSE)
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
alpha_range <- quantile(obs_alpha$V1, probs = c(0.2, 0.8))
probs <- data.frame(alpha = unique(sort(obs_alpha$V1)))

probs$n_clust <- sapply(probs$alpha, function(x) sum(obs_alpha$V1 == x))
probs$obs_prob <- probs$n_clust / sum(probs$n_clust)

probs$in_range <- ifelse(probs$alpha > alpha_range[1] & probs$alpha < alpha_range[2], 1, 0)

probs$obs_prob_in <- 0
wh <- which(probs$in_range == 1)
probs$obs_prob_in[wh] <- probs$obs_prob[wh]
probs$obs_prob_in <- probs$obs_prob_in / sum(probs$obs_prob_in)

probs$hypoth_prob <- 0
wh <- which(probs$in_range == 1 & probs$alpha > median(obs_alpha$V1))
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
  theme(panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12))



alpha <- probs$alpha[probs$in_range == 1]
p1 <- probs$obs_prob_in[probs$in_range == 1]
p2 <- probs$hypoth_prob[probs$in_range == 1]




# ----------- Calculating the IPW ----------- #

trt_col <- which(names(dta) == 'Trt')
out_col <- which(names(dta) == out_name)


yhat_group <- GroupIPW(dta = dta, cov_cols = cov_cols, phi_hat = phi_hat,
                       alpha = alpha, trt_col = trt_col,
                       out_col = out_col)$yhat_group
yhat_pop <- Ypop(ygroup = yhat_group)$ypop


ps_info_est <- list(glm_form = glm_form, ps_with_re = ps_with_re,
                    gamma_numer = phi_hat[[1]], use_control = TRUE)
boots <- BootVar(dta = dta, B = B, alpha = alpha, ps = 'est',
                 cov_cols = cov_cols, ps_info_est = ps_info_est,
                 verbose = TRUE, trt_col = trt_col, out_col = out_col,
                 return_everything = FALSE)



# --------- Asymptotic variance  ----------- #

scores <- CalcScore(dta = dta, neigh_ind = NULL, phi_hat = phi_hat,
                    cov_cols = cov_cols, trt_name = 'Trt')
ie_var <- IEvar(ygroup = yhat_group[, 1, ], ps = 'estimated', scores = scores)

prob_diff <- p2 - p1
ie_Falpha <- c(est = sum(prob_diff * yhat_pop[1, ]),
               var = delta_method(ie_var, vec = prob_diff))
rep(as.numeric(ie_Falpha[1]), 3) +
  c(- 1, 0, 1) * qnorm(1 - alpha_level / 2) * sqrt(ie_Falpha[2])


# ------ Bootstrap confidence interval ------------ #

boots_ie <- apply(boots[1, , ], 2, function(x) sum(x * prob_diff))
c(ie_Falpha[1], quantile(boots_ie, probs = c(0, 1) + c(1, - 1) * alpha_level / 2))



if (save_results) {
  res <- list(yhat_group = yhat_group, boots = boots, scores = scores,
              ie_var = ie_var, B = B)
  save(res, file = paste0(save_path, 'results1_Falpha.dat'))
}

