
############################# Partial #############################
############################# Normal #############################
# ps(w)
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)
# the original PLC(px)
PLCf <- function(px, c=2.64, d=3.54)1-exp(-(-px/d)^c)
# PLC(w)
PLCwf <- function(w)PLCf(psf(w))*100
# modified gsmax
gsmaxfm <- function(w, wL,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54){
  
  pxL <- psf(wL)
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  ps <- psf(w)
  f1 <- function(x)(ps-x)*h2*kxfm(x)/(h*VPD)
  
  res <- ifelse(pxL<ps, optimize(f1, c(pxL,ps), tol=.Machine$double.eps, maximum=T)$objective, 0)
  return(res)
}

# Af(gs)
Af <- function(gs, ca, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# modified mf(w, gs)
mfm <- function(w, gs, wL,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                h2=l*LAI/nZ*p/1000, kxmax=5, c=2.64, d=3.54, h3=10){
  
  pxL <- psf(wL)
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- psf(w)
    f1 <- function(x)(ps-x)*h2*kxfm(x)
    res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, ps)
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- psf(w)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxfm(x)-h*VPD*gs)^2
    res <- ifelse(pxmin<ps, optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum, ps)
    return(res)
  }
  
  px <- pxf(w, gs)
  kx <- kxfm(px)
  PLC <- 1-kx/kxmax
  res <- h3*(PLC-PLCfm(0))
  return(res)
}

# modified B(w, gs)
Bfm <- function(w, gs, wL, ca)Af(gs, ca)-mfm(w, gs, wL)

# family ESS
gswLf <- function(w, wL, ca){
  Bfm1 <- function(gs)Bfm(w, gs, wL, ca)
  gsmaxfm1 <- function(w)gsmaxfm(w, wL)
  res <- ifelse(0<gsmaxfm1(w), optimize(Bfm1, c(0, gsmaxfm1(w)), tol=.Machine$double.eps, maximum=T)$maximum, 0)
  return(res)
}

# LHS endpoint
wLLf <- function(wL, ca){
  Bfm1 <- function(w)Bfm(w, gswLf(w, wL, ca), wL, ca)
  res <- uniroot(Bfm1, c(wL, 1), tol=.Machine$double.eps)
  return(res)
}

# averB for invader
averBif <- function(wLi, wLr, ca,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  wLLr <- wLLf(wLr, ca)$root
  wLLi <- wLLf(wLi, ca)$root
  gswLfr <- Vectorize(function(w)ifelse(w<wLLr, 0, gswLf(w, wLr, ca)))
  gswLfi <- Vectorize(function(w)ifelse(w<wLLi, 0, gswLf(w, wLi, ca)))
  
  Ef <- function(w){h*VPD*gswLfr(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.25)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLLr)/(1-wLLr)-k*integralrEf(w)*1/(1-wLLr))*1/(1-wLLr)
  f1 <- Vectorize(function(w)Bfm(w, gswLfi(w), wLi, ca)*fnoc(w))
  fLnoc <- 1/k*exp(-k*integrate(rEf, wLLr, 1, rel.tol=.Machine$double.eps^0.25)$value*1/(1-wLLr))
  res1 <- integrate(f1, wLLr, 1, rel.tol=.Machine$double.eps^0.25)
  res2 <- ifelse(wLLi<wLLr, fLnoc*Bfm(wLLr, gswLfi(wLLr), wLi, ca), 0)
  res <- res1$value+res2
  message("Resident;Invader;averBI: ", wLr, " ", wLi, " ", res)
  return(res)
}
