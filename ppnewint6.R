# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the model with capital 
# stock using the New Interpretation
# with wage differentials but no
# unproductive sectors

ppnewint6 <- function(A, l, w, v, Q, D, K, t, l_simple){
  
  # -- Inputs for the function
  # A (nxn): input output matrix
  # l (1Xn): direct labor vector (not adjusted for complexity)
  # l_simple (1xn): direct labor inputs (adjusted for complexity)
  # w (1xn): vector of wage rates
  # v: value of labor power (scalar)
  # Q (nx1): gross output vector
  # D (nxn): depreciation matrix
  # K (nxn): capital stock coefficient matrix
  # t (nxn): turnover diagonal matrix
  
  # Necessary conditon (v<1)
  if(v>=(l_simple%*%diag(w)%*%Q)/(l%*%diag(w)%*%Q)){
    stop("Uniform rate of profit cannot be computed")
  } else{
    
    # Identity matrix 
    I <- diag(ncol(A))
    
    # Net output
    y <- (I-A-D)%*%Q
    
    # -- Maximum eigenvalue of N
    N <- (K + A%*%t)%*%solve(I-A-D)
    maxEigenv <- max(Mod(eigen(N)$values))
    
    # Is N nonnegative?
    nn_N <- ifelse(min(N)>=0,1,0)
    # Is N irreducible?
    require(popdemo)
    ir_N <- ifelse(popdemo::isIrreducible(N),1,0)
    
    # -- Maximum rate of profit
    R <- (1/maxEigenv)
    
    
    # ----- Solve for uniform rate of profit
    
    # -- Define Univariate function of rate of profit
    myfunc <- function(r2){
      
      B1 <- solve(I - A - D - r2*K - r2*A%*%t)
      C1 <- (l%*%diag(w) + r2*l%*%diag(w)%*%t)
      
      return(
        (C1%*%B1%*%y) - ((l_simple%*%diag(w))%*%Q)/v
      )
    }
    
    # Find root to get uniform rate of profit
    # Note: upper bound should be kept less than
    # R because the function blows up at R
    r <- uniroot(myfunc,c(0,R-0.01))$root
    
    # ----- Solve for price of production vector
    p_abs <- (l%*%diag(w) + r*l%*%diag(w)%*%t)%*%solve(I - A - D -r*K - r*A%*%t)
    
    # Vector of values 
    # Note: we use the labor input adjusted for complexity
    lambda <- l_simple%*%solve(I - A -D)
    colnames(lambda) <- colnames(l_simple)
    
    # MEV
    mev <- (p_abs%*%y)/(l %*%Q)
    
    # Results as a list
    return(list("Max Eigen Value (N)" = maxEigenv,
                "Maximal Rate of Profit" = R,
                "Uniform Rate of Profit" = r,
                "Prices of Production (Absolute)" = p_abs,
                "Values" = lambda,
                "Monetary Expression of Value" = mev[1,1],
                "N: Nonnegative (1=Y,0=N)" = nn_N,
                "N: Irreducible (1=Y,0=N)" = ir_N
    )
    )
  }

}
