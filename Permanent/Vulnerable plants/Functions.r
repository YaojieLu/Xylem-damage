
################################ Permanent ################################
############################ Vulnerable plants ############################
library("pracma")

# xylem conductance function
kxf <- function(px, kxmax=5, c=9.53, d=1.28)kxmax*exp(-(-px/d)^c)
#kxf <- function(px, kxmax=5, c=2.64, d=3.54)kxmax*exp(-(-px/d)^c)

# PLC(px)
PLCf <- function(px)1-kxf(px)/kxmax

# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# pxmin(w)
pxminf <- function(w){
  ps <- psf(w)
  f1 <- function(x)(ps-x)*kxf(x)#*h2
  res <- optimize(f1, c(-20, ps), tol=.Machine$double.eps, maximum=TRUE)
  return(res$maximum)
}

# original gsmax
gsmaxf <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000){
  pxmin <- pxminf(w)
  kxmin <- kxf(pxmin)
  res <- (psf(w)-pxmin)*h2*kxmin/(h*VPD)
  return(res)
}

# family ESS
gswLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000){
  pxmin <- psf(wL)
  kxmin <- kxf(pxmin)
  res <- (psf(w)-pxmin)*h2*kxmin/(h*VPD)
  return(res)
}

# family ESS as px(w)
pxwLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000){
  ps <- psf(w)
  pxmin <- psf(wL)
  kxmin <- kxf(pxmin)
  res <- ps-h*VPD*gswLf(w, wL)/kxmin/h2
  return(res)
}

# family ESS as px(ps)
pxpswLf <- function(ps, wL, pe=-1.58*10^-3, b=4.38){
  w <- (ps/pe)^(-1/b)
  res <- pxwLf(w, wL)
  return(res)
}

# A
Af <- function(gs, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# averBi
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1,
                    h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswLfr <- Vectorize(function(w)ifelse(w<wLr, 0, gswLf(w, wLr)))
  gswLfi <- Vectorize(function(w)ifelse(w<wLi, 0, gswLf(w, wLi)))
  
  Ef <- function(w)h*VPD*gswLfr(w)
  rEf <- function(w)1/Ef(w)
  integralrEf <- Vectorize(function(w)quadinf(rEf, w, 1)$Q)
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLr)/(1-wLr)-k*1/(1-wLr)*integralrEf(w))*1/(1-wLr)
  fA <- function(w)Af(gswLfi(w))*cPDF*fnoc(w)
  
  fLnoc <- 1/k*exp(-k*1/(1-wLr)*integralrEf(wLr))
  
  cPDF <- 1/(fLnoc+quadinf(fnoc, wLr, 1)$Q)
  fL <- cPDF*fLnoc
  res1 <- quadinf(fA, wLr, 1)
  res2 <- fL*Af(gswLfi(wLr))
  res <- res1$Q+res2
  return(res)
}

# [0, wLr]
optbelowf <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(0.207, wLr), tol=.Machine$double.eps, maximum=T)
  return(res)
}

# [wLr, 1]
optabovef <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(wLr, 0.999), tol=.Machine$double.eps, maximum=T)
  return(res)
}

# fL
fLf <- function(wL,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1,
                h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000,
                gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswf <- Vectorize(function(w)gswLf(w, wL))
  Ef <- function(w)h*VPD*gswf(w)
  rEf <- function(w)1/Ef(w)
  integralrEf <- Vectorize(function(w)quadinf(rEf, w, 1)$Q)
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wL)/(1-wL)-k*1/(1-wL)*integralrEf(w))*1/(1-wL)
  
  fLnoc <- 1/k*exp(-k*1/(1-wL)*integralrEf(wL))
  cPDF <- 1/(fLnoc+quadinf(fnoc, wL, 1)$Q)
  fL <- cPDF*fLnoc
  browser()
  return(fL)
}
