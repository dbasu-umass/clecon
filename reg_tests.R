# Function to compute deviation between
# vector of relative price of production
# and vector of relative value using
# regression-based measures

reg_tests <- function(x,y){
  
  # -- Input variables
  # x: price vector
  # y: value vector
  
  # --- Testing
  # x <- res_si1$`Prices of Production (Absolute)`
  # y <- res_si1$Values
  
  # Remove any zero prices
  mydat <- data.frame(x=as.vector(x),y=as.vector(y))
  mydat1 <- mydat[mydat$x!=0, ]
  
  # All possible relative prices
  x2 <- combn(mydat1$x, 2)
  relp <- x2[1,]/x2[2,]
  
  # All possible relative values
  y2 <- combn(mydat1$y, 2)
  relv <- y2[1,]/y2[2,]
  
  # Run regression
  dev1 <- lm(log(relp) ~ log(relv))
  dev2 <- lm(relp ~ relv)
  
  # Conduct F-test
  l1 <- linearHypothesis(dev1, c("(Intercept)=0","log(relv)=1"))
  l2 <- linearHypothesis(dev2, c("(Intercept)=0","relv=1"))
  
  # Result
  reg_result <- c(
    summary(dev1)$r.squared,
    l1$F[2],
    round(l1$`Pr(>F)`[2], digits=4),
    summary(dev2)$r.squared,
    l2$F[2],
    round(l2$`Pr(>F)`[2], digits=4)
  )
  
  # Return result
  return(reg_result)
}
