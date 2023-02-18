# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the circulating capital 
# model using the New Interpretation
# Marx's labor theory of value and
# allowing for wage differential across
# industries

ppnewint2 <- function(A, l, w, v, Q){
  
  # -- Inputs
  # A (nxn): input output matrix
  # l (1Xn): direct labor vector
  # w (1xn): nominal wage rate vector
  # v: value of labor power (scalar)
  # Q (nx1): gross output vector
  
  # Identity matrix 
  I <- diag(ncol(A))
  
  # Net output
  y <- (I-A)%*%Q
  
  # -- Maximum eigenvalue of A
  maxEigenv <- max(Mod(eigen(A)$values))
  
  # -- Maximum rate of profit
  R <- (1/maxEigenv)-1
  
  
  # ----- Solve for uniform rate of profit
  
  # -- Define Univariate Function of rate of profit
  myfunc <- function(r2){
    
    return(
      (1+r2)*(l%*%diag(w))%*%solve(I-(1+r2)*A)%*%y - ((l%*%diag(w))%*%Q)/v
    )
  }
  
  # Find root to get uniform rate of profit
  # Note: upper bound should be kept less than
  # R because the function blows up at R
  r <- uniroot(myfunc,c(0,R-0.01))$root
  
  # ----- Solve for price of production vector
  p_abs <- (1+r)*(l%*%diag(w))%*%solve(I-(1+r)*A)
  colnames(p_abs) <- colnames(l)
  
  # Vector of values 
  lambda <- l%*%solve(I - A)
  colnames(lambda) <- colnames(l)
  
  # MEV
  mev <- (p_abs%*%y)/(l %*%Q)
  
  # Monetary expression of value (using gross output)
  mev_gross <- (p_abs%*%Q)/(lambda%*%Q)
  
  
  return(list("Max Eigen Value (A)" = maxEigenv,
              "Maximal Rate of Profit" = R,
              "Uniform Rate of Profit" = r,
              "Prices of Production (Absolute)" = p_abs,
              "Prices of Production (Relative)" = p_abs/p_abs[1],
              "Values" = lambda,
              "Monetary Expression of Value" = mev[1,1],
              "Monetary Expression of Value (Gross)" = mev_gross[1,1])
  )
  
}

