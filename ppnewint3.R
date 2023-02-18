# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the circulating capital 
# model using the New Interpretation
# Marx's labor theory of value and
# allowing for some unproductive industries

ppnewint3 <- function(A, Ap, l, lp, w, v, Q, Qp){
  
  # -- Inputs to the function
  # A (nxn): input output matrix
  # Ap (mxm): input output matrix for productive sectors
  # l (1Xn): direct labor vector
  # lp (1Xm): direct labor vector for productive sectors
  # w: average wage rate (scalar)
  # v: value of labor power (scalar)
  # Q (nx1): gross output vector
  # Qp (mx1): gross output vector for productive sectors
  

  if(v>=(lp%*%Qp)/(l%*%Q)){
    stop("Uniform rate of profit cannot be computed")
  } else{
    
    # Identity matrix 
    I <- diag(ncol(A))
    Ip <- diag(ncol(Ap))
    
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
        (1+r2)*l%*%(solve(I-(1+r2)*A))%*%y - (lp%*%Qp)/v
      )
    }
    
    # Find root to get uniform rate of profit
    # Note: upper bound should be kept less than
    # R because the function blows up at R
    r <- uniroot(myfunc,c(0,R-0.01))$root
    
    # ----- Solve for price of production vector
    p_abs <- (1+r)*(w*l)%*%solve(I-(1+r)*A)
    
    # Vector of values 
    lambda <- lp%*%solve(Ip - Ap)
    
    # MEV
    mev <- (p_abs%*%y)/(lp %*%Qp)
    
    
    return(list("Max Eigen Value (A)" = maxEigenv,
                "Maximal Rate of Profit" = R,
                "Uniform Rate of Profit" = r,
                "Prices of Production (Absolute)" = p_abs,
                "Prices of Production (Relative)" = p_abs/p_abs[1],
                "Values" = lambda,
                "Monetary Expression of Value" = mev[1,1])
    )
  }
  
  
}

