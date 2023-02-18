# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the model with capital 
# stock using the New Interpretation

ppnewint5 <- function(A, l, w, v, Q, D, K, t){
  
  # -- Inputs for the function
  # A (nxn): input output matrix
  # l (1Xn): direct labor vector
  # w: average wage rate (scalar)
  # v: value of labor power (scalar)
  # Q (nx1): gross output vector
  # D (nxn): depreciation matrix
  # K (nxn): capital stock coefficient matrix
  # t (nxn): turnover diagonal matrix
  
  # Identity matrix 
  I <- diag(ncol(A))
  
  # Net output
  y <- (I-A-D)%*%Q
  
  # -- Maximum eigenvalue of N
  N <- (K + A%*%t)%*%solve(I-A-D)
  maxEigenv <- max(Mod(eigen(N)$values))
  
  # -- Maximum rate of profit
  R <- (1/maxEigenv)
  
  
  # ----- Solve for uniform rate of profit
  
  # -- Define Univariate function of rate of profit
  myfunc <- function(r2){
    
    B1 <- solve(I - A - D - r2*K - r2*A%*%t)
    C1 <- (w*l + r2*w*l%*%t)
    
    return(
      (C1%*%B1%*%y) - ((w*l)%*%Q)/v
    )
  }
  
  # Find root to get uniform rate of profit
  # Note: upper bound should be kept less than
  # R because the function blows up at R
  r <- uniroot(myfunc,c(0,R-0.01))$root
  
  # ----- Solve for price of production vector
  p_abs <- (w*l + r*w*l%*%t)%*%solve(I - A - D -r*K - r*A%*%t)
  
  # Vector of values 
  lambda <- l%*%solve(I - A -D)
  
  # MEV
  mev <- (p_abs%*%y)/(l %*%Q)
  
  
  return(list("Max Eigen Value (N)" = maxEigenv,
              "Maximal Rate of Profit" = R,
              "Uniform Rate of Profit" = r,
              "Prices of Production (Absolute)" = p_abs,
              "Values" = lambda,
              "Monetary Expression of Value" = mev[1,1])
  )
  
}

