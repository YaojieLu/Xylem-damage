
################################ Permanent ################################
############################ Vulnerable plants ############################
options(digits=20)
source("Permanent/Vulnerable plants/Functions.r")

# Minimize f1 and f2
f3 <- function(wLr){
  res1 <- optbelowf(wLr)$maximum
  res2 <- optabovef(wLr)$maximum
  res <- abs(res1-wLr)+abs(res2-wLr)
  return(res)
}
# Results
ca <- 400
k <- 0.05
MAP <- 222

res <- optimize(f3, c(0.215, 0.225), tol=.Machine$double.eps)
res
optbelowf(res$minimum)
optabovef(res$minimum)
