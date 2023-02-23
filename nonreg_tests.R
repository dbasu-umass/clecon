# Function to compute deviation between
# vector of relative price of production
# and vector of relative value using
# non-regression-based measures

nonreg_tests <- function(x,y,w,w_avg,mev,Q){
  
  # ---- Inputs to the function
  # x (1xn): price vector
  # y (1xn): value vector 
  # w (1xn): vector of nominal wage rate
  # w_avg: (scalar) nominal wage rate
  # mev: (scalar) monetary expression of value (using gross output)
  # Q (1xn): vector of gross output
  
  # Remove observations corresponding to zero prices
  mydat <- data.frame(
    x=as.vector(x),
    y=as.vector(y),
    w=as.vector(w),
    Q=as.vector(Q)
    )
  mydat1 <- mydat[mydat$x!=0, ]

  # ---------- Relative vectors (All Combinations) -------- #
  
  # All possible relative prices
  x2 <- combn(mydat1$x, 2)
  relp_all <- x2[1,]/x2[2,]
  
  # All possible relative values
  y2 <- combn(mydat1$y, 2)
  relv_all <- y2[1,]/y2[2,]
  
  # ------------- Measures ---------------------- #
  # --- RMSE%
  rmse_rel_all <- sqrt(mean(((relp_all/relv_all)-1)^2))
  
  # --- Minimum Absolute Distance
  mad_rel_all <- mean(abs((relp_all/relv_all)-1))
  
  # --- Classical distance measure (CDM)
  # Relative wage vector
  w_rel <- mydat1$w/w_avg
  w_comb <- combn(w_rel, 2)
  rel_w <- w_comb[1,]/w_comb[2,]
  
  # d vector
  d_j <- w_rel * mydat1$y
  d_comb <- combn(d_j, 2)
  rel_d <- d_comb[1,]/d_comb[2,]
  
  # Vector of weights
  omega_j <- mydat1$Q/sum(mydat1$Q)
  omega_2 <- combn(omega_j, 2)
  rel_omega <- omega_2[1,]*omega_2[2,]
  
  # CDM
  cdm_rel_all <- sum( abs((relp_all/rel_d)-1) * rel_omega )
  
  # --- Mean Absolute Weighted Distance
  mawd_rel_all <- sum(abs((relp_all/relv_all)-1)*rel_omega)
  
  # --- Angle in degrees
  z <- relp_all/relv_all
  tan_alpha <- (sd(z)/mean(z))
  alpha_rel_all <- (atan(tan_alpha))*(180/pi)
  
  # -- Distance using angle
  d_rel_all <- sqrt(2*(1-cos(alpha_rel_all*(pi/180))))
  
  
  # ---- Results ------- #
  # Return result
  return(
    list(
      "RMSE" = rmse_rel_all, 
      "MAD" = mad_rel_all,
      "MAWD" = mawd_rel_all, 
      "CDM" = cdm_rel_all,
      "Angle" = alpha_rel_all,
      "Distance (using angle)" = d_rel_all
    )
  )
  
}
