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


n_neigh <- 50
hierarchical_method <- 'ward.D2'
coord_names <- c('Fac.Longitude', 'Fac.Latitude')
trt_name <- 'SnCR'
out_name <- 'mean4maxOzone'
estimand <- '1'  # Set to '1' for the estimand depending on covariates.
                 # Set to '2' for pi-estimand.

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

# ----------- Calculating the IPW ----------- #

trt_col <- which(names(dta) == 'Trt')
out_col <- which(names(dta) == out_name)
alpha_range <- quantile(obs_alpha$V1, probs = c(0.2, 0.8))
alpha <- seq(alpha_range[1], alpha_range[2], length.out = 40)
alpha <- sort(c(alpha, 0.1, 0.4))

yhat_group <- GroupIPW(dta = dta, cov_cols = cov_cols, phi_hat = phi_hat,
                       alpha = alpha, trt_col = trt_col, out_col = out_col,
                       estimand = estimand)
yhat_group <- yhat_group$yhat_group

scores <- CalcScore(dta = dta, neigh_ind = NULL, phi_hat = phi_hat,
                    cov_cols = cov_cols, trt_name = 'Trt')
ypop <- Ypop(ygroup = yhat_group, ps = 'estimated', scores = scores, dta = dta)

yhat_pop <- ypop$ypop
yhat_pop_var <- ypop$ypop_var

# --------- Direct effect ----------- #
de <- DE(ypop = yhat_pop, ypop_var = yhat_pop_var, alpha = alpha)


# --------- Indirect effect ----------- #
ie <- IE(ygroup = yhat_group[, 1, ], ps = 'estimated', scores = scores)



# ------------  PLOTTING ------------ #

de_plot <- data.frame(alpha = alpha, de = de[1, ], low = de[3, ],
                      high = de[4, ])
a1 <- which(alpha == 0.1)
ie1_plot <- as.data.frame(t(ie[c(1, 3, 4), a1, ]))
a1 <- which(alpha == 0.4)
ie2_plot <- as.data.frame(t(ie[c(1, 3, 4), a1, ]))

res_array <- array(NA, dim = c(length(alpha), 3, 3))
dimnames(res_array) <- list(alpha = alpha, quant = c('DE', 'IE1', 'IE2'),
                            stat = c('est', 'LB', 'HB'))
res_array[, 1, ] <- as.matrix(de_plot[, - 1])
res_array[, 2, ] <- as.matrix(ie1_plot)
res_array[, 3, ] <- as.matrix(ie2_plot)

res_df <- plyr::adply(res_array[, , 1], 1 : 2)
res_df$LB <- c(de_plot$low, ie1_plot$LB, ie2_plot$LB)
res_df$UB <- c(de_plot$high, ie1_plot$UB, ie2_plot$UB)

f_names <- list('DE' = expression(DE(alpha)),
                'IE1' = expression(IE(0.1,alpha)),
                'IE2' = expression(IE(0.4,alpha)))
f_labeller <- function(variable, value){
  return(f_names[value])
}
res_df$alpha <- as.numeric(levels(res_df$alpha))[res_df$alpha]


ggplot(data = res_df, aes(x = alpha, y = V1, group = quant)) +  geom_line() +
  facet_wrap(~ quant, nrow = 1, labeller = f_labeller) +
  geom_ribbon(data = res_df, aes(ymin = LB, ymax = UB, group = quant), alpha = 0.3) +
  xlab(expression(alpha)) + ylab('') +
  theme(axis.title = element_text(size = 12),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0.1, 0.4, by = 0.1)) +
  geom_hline(yintercept = 0, linetype = 2)



library(plot3D)
par(mar = c(1, 1, 1, 3))
persp3D(alpha, alpha, ie[1, , ], theta=50, phi=50, axes=TRUE,
        nticks=5, ticktype = "detailed", xlab='α1', ylab='α2', zlab='',
        colkey = list(length = 0.5, shift = -0.1))



