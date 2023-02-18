# Function to compute deviation between
# vector of relative price of production
# and vector of relative value using
# regression-based measures

reg_tests <- function(x,y){
  
  # -- Input variables
  # x: price vector
  # y: value vector
  
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
  
   # Return result
  return(
    list(
      "RSquared (log-log)" = summary(dev1)$r.squared,
      "FStat (log-log)" = l1$F[2],
      "Pvalue (log-log)" = round(l1$`Pr(>F)`[2], digits=4),
      "RSquared (level-level)" = summary(dev2)$r.squared,
      "FStat (level-level)" = l2$F[2],
      "Pvalue (level-level)" = round(l2$`Pr(>F)`[2], digits=4)
    )
  )
}
