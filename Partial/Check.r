
############################# Partial #############################
options(digits=20)
source("Partial/gswL(w).r")

# averBi
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  wLLr <- wLLf(wLr)$root
  wLLi <- wLLf(wLi)$root
  gswLfr <- Vectorize(function(w)ifelse(w>wLLr, gswLf(w, wLr), 0))
  gswLfi <- Vectorize(function(w)ifelse(w>wLLi, gswLf(w, wLi), 0))
  
  Ef <- function(w){h*VPD*gswLfr(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.25)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLLr)/(1-wLLr)-k*integralrEf(w)*1/(1-wLLr))*1/(1-wLLr)
  f1 <- Vectorize(function(w)Bfm(w, gswLfi(w), wLi)*fnoc(w))
  fLnoc <- 1/k*exp(-k*integrate(rEf, wLLr, 1, rel.tol=.Machine$double.eps^0.25)$value*1/(1-wLLr))
  res1 <- integrate(f1, wLLr, 1, rel.tol=.Machine$double.eps^0.25)
  res2 <- ifelse(wLLi<wLLr, fLnoc*Bfm(wLLr, gswLfi(wLLr), wLi), 0)
  res <- res1$value+res2
  browser()
  return(res)
}

# Results
ca <- 400
k <- 0.025
MAP <- 365
pkx <- 0.5
wLr <- 0.135
averBif(wLr-1e-10, wLr)
