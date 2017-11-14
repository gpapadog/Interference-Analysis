Old_IPW <- function(alpha, dta, neigh_ind, phi_hat, cov_cols) {
  
  n_neigh <- length(neigh_ind)
  
  yhat_group <- array(NA, dim = c(n_neigh, 2, length(alpha)))
  dimnames(yhat_group) <- list(neigh = 1 : n_neigh, it = c(0, 1), alpha = alpha)
  
  for (aa in 1 : length(alpha)) {
    curr_alpha <- alpha[[aa]]
    print(paste0('alpha = ', curr_alpha))
    
    for (it in c(0, 1)) {
      curr_it <- it
      
      for (nn in 1 : n_neigh) {
        
        # We only calculate group average if there are individuals with that treatment.
        if (any(dta$A[neigh_ind[[nn]]] == curr_it)) {
          y_curr <- 0
          
          bern_prob <- curr_alpha ^ curr_it * (1 - curr_alpha) ^ (1 - curr_it)
          
          for (ind in neigh_ind[[nn]]) {
            y_curr <- y_curr + (dta$A[ind] == curr_it) * dta$Y[ind] / bern_prob
          }
          
          # TRUE PS.
          denom <- DenomIntegral(A = dta$A[neigh_ind[[nn]]],
                                 X = dta[neigh_ind[[nn]], cov_cols],
                                 phi_hat = phi_hat, alpha = curr_alpha)
          denom <- length(neigh_ind[[nn]]) * denom$value
          yhat_group[nn, it + 1, aa] <- y_curr / denom 
          
        }
      }
    }
  }
  
  return(yhat_group = yhat_group)
}

