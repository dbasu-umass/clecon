# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the model with capital 
# stock using the New Interpretation
# with a uniform wage rate and allowing for
# unproductive sectors

ppnewint7 <- function(A, Ap, l, lp, w, v, Q, Qp, D, Dp, K, t, lp_simple){
  
  # -- Inputs for the function
  # A (nxn): input output matrix
  # Ap (mxm): input output matrix for productive industries
  # l (1Xn): direct labor vector (not adjusted for complexity)
  # lp (1Xm): direct labor vector for productive industries (not adjusted for complexity)
  # lp_simple (1Xm): direct labor vector for productive industries (adjusted for complexity)
  # w: average wage rate (scalar)
  # v: value of labor power (scalar)
  # Q (nx1): gross output vector
  # Qp (mx1): gross output vector for productive industries
  # D (nxn): depreciation matrix
  # Dp (mxm): depreciation matrix for productive industries
  # K (nxn): capital stock coefficient matrix
  # t (nxn): turnover diagonal matrix
  
  # Necessary condition for solutions
  if(v>=(lp_simple%*%Qp)/(l%*%Q)){
    stop("Necessary condition violated. Uniform rate of profit cannot be computed")
  } else{
    # Identity matrix 
    I <- diag(ncol(A))
    Ip <- diag(ncol(Ap))
    
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
      C1 <- (w*l + r2*w*l%*%t)
      
      return(
        (C1%*%B1%*%y) - ((w*lp_simple)%*%Qp)/v
      )
    }
    
    # Find root to get uniform rate of profit
    # Note: upper bound should be kept less than
    # R because the function blows up at R
    r <- uniroot(myfunc,c(0,(R-0.00001)))$root
    
    # ----- Solve for price of production vector
    p_abs <- (w*l + r*w*l%*%t)%*%solve(I - A - D -r*K - r*A%*%t)
    
    # Vector of values 
    # Note: we use the labor input adjusted for complexity
    lambda <- lp_simple%*%solve(Ip - Ap -Dp)
    colnames(lambda) <- colnames(lp_simple)
    
    # MEV
    mev <- (p_abs%*%y)/(lp %*%Qp)
    
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
    
  } # end of if else statement
  
  
}
