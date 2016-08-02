
################################ Permanent ################################
############################ Vulnerable plants ############################
options(digits=20)
source("Permanent/Vulnerable plants/Functions.r")

f1 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  resL <- optimize(averBif1, c(0.1, wLr), tol=.Machine$double.eps, maximum=T)
  resH <- optimize(averBif1, c(wLr, 0.3), tol=.Machine$double.eps, maximum=T)
  res <- resH$maximum-resL$maximum
  message(wLr)
  return(res)
}
f2 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  resL <- optimize(averBif1, c(0.1, wLr), tol=.Machine$double.eps, maximum=T)
  resH <- optimize(averBif1, c(wLr, 0.3), tol=.Machine$double.eps, maximum=T)
  return(list(resL=resL, resH=resH))
}
# Results
ca <- 400
k <- 0.1
MAP <- 500

ESS <- optimize(f1, c(0.207, 0.23), tol=.Machine$double.eps)
ESS
res1 <- f2(ESS$minimum)
res1
res2 <- averBif(ESS$minimum, ESS$minimum)
res2
res1$resL$objective-res2
res1$resH$objective-res2
