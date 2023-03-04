# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the circulating capital 
# model using the Standard Interpretation
# of Marx's labor theory of value
# with uniform wage rates and allowing for
# unproductive industries

ppstdint1 <- function(A, Ap, l, b, Q, Qp, lp_simple){
  
  # -- Inputs to the function
  # A (nxn): input output matrix
  # Ap (mxm): input output matrix for productive industries
  # l (1xn): direct labor vector (not adjusted for complexity)
  # l_simple (1xn): direct labor vector (adjusted for complexity)
  # lp_simple (1xn): direct labor vector for productive industries (adjusted for complexity)
  # b (nx1): real wage vector
  # Q (nx1): gross output vector
  # Qp (mx1): gross output vector for productive industries
  
  # ---- M
  M <- A + b%*%l
  
  # Is M nonnegative?
  nn_M <- ifelse(min(M)>=0,1,0)
  # Is M irreducible?
  require(popdemo)
  ir_M <- ifelse(popdemo::isIrreducible(M),1,0)
  
  # ---- Uniform rate of profit
  maxEigenv <- max(Mod(eigen(M)$values))
  r <- (1/maxEigenv)-1
  
  # -- Maximal rate of profit (when b is the 0 vector)
  R <- 1/(max(Mod(eigen(A)$values)))-1
  
  # ---- Relative price of production vector
  # First column of eigenvector matrix of M
  # The vector has all real elements (of the same sign)
  p_rel <- Re(eigen(M)$vectors[,1])
  
  # ---- Vector of values (for productive sectors only)
  # Note: we use the labor input adjusted for complexity
  Ip <- diag(ncol(Ap))
  lambda_p <- lp_simple%*%solve(Ip - Ap)
  colnames(lambda_p) <- colnames(lp_simple)
  
  # Normalization using gross output
  mev <- (p_rel%*%Q)/(lambda_p%*%Qp)
  
  # ----- Absolute price of production vector
  p_abs <- mev[1,1]*matrix(p_rel,nrow=1)
  colnames(p_abs) <- colnames(l)
  
  
  # Monetary expression of value (using gross output)
  mev_gross <- (p_abs%*%Q)/(lambda_p%*%Qp)
  
  # ----- Results as a list
  return(list("Max Eigen Value (M)" = maxEigenv,
              "Uniform Rate of Profit" = r,
              "Maximal Rate of Profit" = R,
              "Prices of Production (Absolute)" = p_abs,
              "Prices of Production (Relative)" = p_rel,
              "Values (Productive Industries)" = lambda_p,
              "Monetary Expression of Value (Gross)" = mev_gross[1,1],
              "M: Nonnegative (1=Y,0=N)" = nn_M,
              "M: Irreducible (1=Y,0=N)" = ir_M
  )
  )
}

