# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the circulating capital 
#  model using the Sraffian approach

ppsraffa1 <- function(A, l, Q, w, pshare){
  
  # -- Inputs to the function
  # A (nxn): input output matrix
  # l (1xn): direct labor vector
  # Q (nx1): gross output vector
  # w: average wage rate (scalar)
  # pshare : profit share (scalar)
  
  # ---- M
  I <- diag(ncol(A))
  H <- A%*%solve(I-A)
  
  # ---- Maximal rate of profit
  maxEigenv <- max(Mod(eigen(H)$values))
  R <- (1/maxEigenv)
  
  # ---- Vector of values 
  lambda <- l%*%solve(I - A)
  colnames(lambda) <- colnames(l)
  
  # ---- Price of production vector
  p_abs <- w*lambda%*%solve(I - R*pshare*H)
  colnames(p_abs) <- colnames(lambda)
  
  # Monetary expression of value (using gross output)
  mev_gross <- (p_abs%*%Q)/(lambda%*%Q)
  
  # ----- Output of function
  return(list("Max Eigen Value (H)" = maxEigenv,
              "Maximal Rate of Profit" = R,
              "Prices of Production (Absolute)" = p_abs,
              "Values" = lambda,
              "Monetary Expression of Value (Gross)" = mev_gross[1,1])
  )
}
