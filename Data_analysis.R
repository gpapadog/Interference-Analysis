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

# ----------- Calculating the IPW ----------- #

trt_col <- which(names(dta) == 'Trt')
out_col <- which(names(dta) == out_name)
alpha_range <- quantile(obs_alpha$V1[! (obs_alpha$V1 %in% c(0, 1))], probs = c(0.2, 0.8))
alpha <- seq(alpha_range[1], alpha_range[2], length.out = 40)

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


# --------- Direct effect ----------- #

de <- DE(ypop = yhat_pop, ypop_var = yhat_pop_var, alpha = alpha)
de <- rbind(de, low_int = de[1, ] - 1.96 * sqrt(de[2, ]))
de <- rbind(de, high_int = de[1, ] + 1.96 * sqrt(de[2, ]))

plot(alpha, de[1, ], ylim = range(de[c(1, 3, 4), ]))
lines(alpha, de[3, ])
lines(alpha, de[4, ])
abline(h = 0, col = 'red')



# --------- Indirect effect ----------- #

ie_var <- IEvar(ygroup = yhat_group[, 1, ], alpha = alpha, ps = 'estimated',
                scores = scores)
ie <- IE(ypop = yhat_pop[1, ], ypop_var = ie_var, alpha = alpha)

par(mfrow = c(3, 3), mar = rep(2, 4))
for (ii in 1 : length(alpha)) {
  print(sum(ie[3, ii, ] < 0 & ie[4, ii, ] > 0))
  plot(alpha, ie[1, ii, ], ylim = range(ie[c(1, 3, 4), ii, ]),
       type = 'l', lwd = 2, main = alpha[ii])
  lines(alpha, ie[3, ii, ])
  lines(alpha, ie[4, ii, ])
  abline(h = 0, col = 'red')
}



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
  
  ie_plot <- as.data.frame(t(ie[c(1, 3, 4), ii, ]))
  names(ie_plot) <- c('ie', 'low', 'high')

  g_ie1 <- ggplot() + geom_line(aes(alpha, ie), data = ie_plot) +
    geom_ribbon(data = ie_plot, aes(x = alpha, ymin = low, ymax = high),
                alpha = 0.3) +
    geom_abline(intercept = 0, slope = 0, linetype = 2) +
    ylab(paste0('IE(', round(alpha[ii], 3), ',', expression(alpha), ')')) +
    xlab(expression(alpha))
  
  g_ie[[length(g_ie) + 1]] <- g_ie1
} 


grid.arrange(grobs = append(list(g_de), g_ie[c(2, 4)]), ncol = 3)
