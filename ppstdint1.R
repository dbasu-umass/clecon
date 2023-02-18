# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the circulating capital 
#  model using the Standard Interpretation
# Marx's labor theory of value

ppstdint1 <- function(A, l, b, Q){
  
  # -- Inputs to the function
  # A (nxn): input output matrix
  # l (1xn): direct labor vector
  # b (nx1): real wage vector
  # Q (nx1): gross output vector
  
  # ---- M
  M <- A + b%*%l
  
  # ---- Uniform rate of profit
  maxEigenv <- max(Mod(eigen(M)$values))
  r <- (1/maxEigenv)-1
  
  # ---- Relative price of production vector
  # First column of eigenvector matrix of M
  # The vector has all real elements (of the same sign)
  p_rel <- Re(eigen(M)$vectors[,1])
  
  # ---- Vector of values 
  I <- diag(ncol(A))
  lambda <- l%*%solve(I - A)
  colnames(lambda) <- colnames(l)
  
  # Normalization using gross output
  mev <- (p_rel%*%Q)/(lambda%*%Q)
  
  # ----- Absolute price of production vector
  p_abs <- mev[1,1]*matrix(p_rel,nrow=1)
  colnames(p_abs) <- colnames(l)
  
  
  # Monetary expression of value (using gross output)
  mev_gross <- (p_abs%*%Q)/(lambda%*%Q)
  
  # ----- Output of function
  return(list("Max Eigen Value (M)" = maxEigenv,
              "Uniform Rate of Profit" = r,
              "Prices of Production (Absolute)" = p_abs,
              "Prices of Production (Relative)" = p_rel,
              "Values" = lambda,
              "Monetary Expression of Value (Gross)" = mev_gross[1,1])
  )
}

