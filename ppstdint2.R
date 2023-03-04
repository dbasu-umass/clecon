# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the model with capital 
# stock using the Standard Interpretation
# with a uniform wage rate and no
# unproductive industries

ppstdint2 <- function(A, l, b, Q, D, K, t, l_simple){
  
  # -- Inputs for the function
  # A (nxn): input output matrix
  # l (1Xn): direct labor vector (not adjusted for complexity)
  # l_simple (1Xn): direct labor vector (adjusted for complexity)
  # b (nx1): vector of consumption
  # Q (nx1): gross output vector
  # D (nxn): depreciation matrix
  # K (nxn): capital stock coefficient matrix
  # t (nxn): turnover diagonal matrix
  
  # Identity matrix 
  I <- diag(ncol(A))
  
  # -- Maximum eigenvalue of N
  N <- (K + (A+b%*%l)%*%t)%*%solve(I-A-D-b%*%l)
  maxEigenv <- max(Mod(eigen(N)$values))
  
  # Is N nonnegative?
  nn_N <- ifelse(min(N)>=0,1,0)
  # Is N irreducible?
  require(popdemo)
  ir_N <- ifelse(popdemo::isIrreducible(N),1,0)
  
  # -- Uniform rate of profit
  r <- (1/maxEigenv)
  
  # -- Maximal rate of profit (when b is the 0 vector)
  M <- (K + (A)%*%t)%*%solve(I-A-D)
  R <- 1/(max(Mod(eigen(M)$values)))
  
  # ----- Solve for price of production vector
  # Rel Price = First column of eigen vector matrix of N
  # The vector has all real elements (of the same sign)
  # If all elements <0, multiply with -1
  p_rel_neg <- (-1)*Re(eigen(N)$vectors[,1])
  p_rel_pos <- Re(eigen(N)$vectors[,1])
  if (Re(eigen(N)$vectors[1,1])<0) {
    p_rel <- p_rel_neg
  }else{
    p_rel <- p_rel_pos
  }
  
  # Vector of values 
  # Note: we use the labor input adjusted for complexity
  lambda <- l_simple%*%solve(I - A -D)
  colnames(lambda) <- colnames(l_simple)
  
  # Normalization using gross output
  mev_num <- (matrix(p_rel,nrow=1)%*%matrix(Q,ncol=1))
  mev_den <- (matrix(lambda,nrow=1)%*%matrix(Q,ncol=1))
  mev <- mev_num/mev_den
  
  # ----- Absolute price of production vector
  p_abs <- mev[1,1]*matrix(p_rel,nrow=1)
  colnames(p_abs) <- colnames(l)
  
  # Results as a list
  return(list("Max Eigen Value (N)" = maxEigenv,
              "Maximal Rate of Profit" = R,
              "Uniform Rate of Profit" = r,
              "Prices of Production (Absolute)" = p_abs,
              "Values" = lambda,
              "Monetary Expression of Value (Gross)" = mev[1,1],
              "N: Nonnegative (1=Y,0=N)" = nn_N,
              "N: Irreducible (1=Y,0=N)" = ir_N
              )
        )
  
}
