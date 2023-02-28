# Function to compute the vector of 
# labor values, the vector of 
# price of production and the uniform 
# rate of profit for the circulating capital 
# model using the New Interpretation of
# Marx's labor theory of value and
# allowing for wage differential across
# industries and unproductive sectors

ppnewint4 <- function(A, Ap, l, lp, w, wp, v, Q, Qp, lp_simple){
  
  # -- Inputs to the function
  # A (nxn): input output matrix
  # Ap (mxm): input output matrix for productive sectors
  # l (1Xn): direct labor vector
  # lp (1Xm): direct labor vector for productive sectors (not adjusted for complexity)
  # lp_simple (1Xm): direct labor vector for productive sectors (adjusted for complexity)
  # w (1xn): vector of wage rates
  # wp (1xm): vector of wage rates for productive sectors
  # v: value of labor power (scalar)
  # Q (nx1): gross output vector
  # Qp (mx1): gross output vector for productive sectors
  
  # Necessary condition for existence of solution
  if(v>=(lp%*%(diag(wp))%*%Qp)/(l%*%(diag(w))%*%Q)){
    stop("Uniform rate of profit cannot be computed")
  } else{
    
    # Identity matrix 
    I <- diag(ncol(A))
    Ip <- diag(ncol(Ap))
    
    # Net output
    y <- (I-A)%*%Q
    
    # -- Maximum eigenvalue of A
    maxEigenv <- max(Mod(eigen(A)$values))
    
    # Is A nonnegative?
    nn_A <- ifelse(min(A)>=0,1,0)
    # Is A irreducible?
    require(popdemo)
    ir_A <- ifelse(popdemo::isIrreducible(A),1,0)
    
    # -- Maximum rate of profit
    R <- (1/maxEigenv)-1
    
    
    # ----- Solve for uniform rate of profit
    
    # -- Define Univariate Function of rate of profit
    myfunc <- function(r2){
      
      return(
        (1+r2)*(l%*%diag(w))%*%(solve(I-(1+r2)*A))%*%y - ((lp%*%diag(wp))%*%Qp)/v
      )
    }
    
    # Find root to get uniform rate of profit
    # Note: upper bound should be kept less than
    # R because the function blows up at R
    r <- uniroot(myfunc,c(0,R-0.01))$root
    
    # ----- Solve for price of production vector
    p_abs <- (1+r)*(l%*%diag(w))%*%solve(I-(1+r)*A)
    
    # Vector of values 
    # Note: we use the labor input adjusted for complexity
    lambda <- lp_simple%*%solve(Ip - Ap)
    colnames(lambda) <- colnames(lp_simple)
    
    # MEV
    mev <- (p_abs%*%y)/(lp %*%Qp)
    
    # Monetary expression of value (using gross output)
    mev_gross <- (p_abs%*%Q)/(lambda%*%Qp)
    
    # Results as a list
    return(list("Max Eigen Value (A)" = maxEigenv,
                "Maximal Rate of Profit" = R,
                "Uniform Rate of Profit" = r,
                "Prices of Production (Absolute)" = p_abs,
                "Prices of Production (Relative)" = p_abs/p_abs[1],
                "Values" = lambda,
                "Monetary Expression of Value" = mev[1,1],
                "Monetary Expression of Value (Gross)" = mev_gross[1,1],
                "A: Nonnegative (1=Y,0=N)" = nn_A,
                "A: Irreducible (1=Y,0=N)" = ir_A
                )
          )
    
  }  
  
  
}
