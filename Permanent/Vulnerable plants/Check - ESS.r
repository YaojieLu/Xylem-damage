
################################ Permanent ################################
############################ Vulnerable plants ############################
options(digits=20)
source("Permanent/Vulnerable plants/Functions.r")

# Auxiliary function
# Optimization among wLi < wLr
f1 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(0.001, wLr), tol=.Machine$double.eps, maximum=T)
  return(res)
}
# Optimization among wLi > wLr
f2 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(wLr, 0.999), tol=.Machine$double.eps, maximum=T)
  return(res)
}

# Minimize f1 and f2
f3 <- function(wLr){
  res1 <- f1(wLr)$maximum
  res2 <- f2(wLr)$maximum
  res <- abs(res1-wLr)+abs(res2-wLr)
  return(res)
}
# Results
ca <- 400
k <- 0.05
MAP <- 222

res <- optimize(f3, c(0.215, 0.225), tol=.Machine$double.eps)
res
f1(res$minimum)
f2(res$minimum)
