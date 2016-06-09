
source("Partial/gswL(w).r")

# averBi
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05,
                    gamma=1/((MAP/365/k)/1000)*nZ){

  wLLr <- wLLf(wLr)$root
  wLLi <- wLLf(wLi)$root
  gswLfr <- Vectorize(function(w)gswLf(w, wLr))
  gswLfi <- Vectorize(function(w)ifelse(w>wLLi, gswLf(w, wLi), 0))
  Bfmi <- Vectorize(function(w)ifelse(gswLfi(w)>0, Bfm(w, gswLfi(w), wLi), 0))
  
  Ef <- function(w){h*VPD*gswLfr(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.4)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLLr)/(1-wLLr)-k*integralrEf(w))
  f1 <- function(w)Bfmi(w)*fnoc(w)
  
  res1 <- integrate(f1, wLLr, 1, rel.tol=.Machine$double.eps^0.4)
  fLnoc <- 1/k*exp(-k*integrate(rEf, wLLr, 1, rel.tol=.Machine$double.eps^0.4)$value)
  res2 <- fLnoc*Bfmi(wLLr)
  res <- res1$value+res2
  message("Resident; Invader; averBI: ", wLr, " ", wLi, " ", res)
  return(res)
}

f1 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res1 <- optimize(averBif1, c(0.135, 0.14), tol=.Machine$double.eps, maximum=T)$maximum
  res <- (res1-wLr)^2
  return(res)
}
f2 <- Vectorize(f1)

# Result
pkx <- 0.5
ca <- 400
k <- 0.05
MAP <- 1825

begin <- proc.time()
optimize(f2, c(0.135, 0.14), tol=.Machine$double.eps)
end <- proc.time()
message(sprintf("completed in %.2f hours", (end[3]-begin[3])/3600))
