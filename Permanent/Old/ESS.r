
# A
Af <- function(gs, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# gswL(w)
gswLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05){
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  pxmin <- pe*wL^(-b)
  res <- (pe*w^-b-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}

# averBi
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswLfr <- function(w)gswLf(w, wLr)
  gswLfi <- function(w)ifelse(w>wLi, gswLf(w, wLi), 0)
  
  Ef <- function(w){h*VPD*gswLfr(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1)$value})#, rel.tol=.Machine$double.eps^0.25
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLr)/(1-wLr)-k*integralrEf(w))
  fLnoc <- 1/k*exp(-k*integrate(rEf, wLr+1e-14, 1)$value)
  f1 <- function(w)Af(gswLfi(w))*fnoc(w)

  res1 <- integrate(f1, wLr, 1)
  res2 <- fLnoc*Af(gswLfi(wLr))
  res <- res1$value+res2
  return(res)
}

f1 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res1 <- optimize(averBif1, c(0.125, 0.132), tol=.Machine$double.eps, maximum=T)$maximum
  return((res1-wLr)^2)
}
f2 <- Vectorize(f1)
