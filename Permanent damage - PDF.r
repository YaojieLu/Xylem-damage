
options (digits=20)

# Function
# f0
f0f <- function(w0,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, h2=h/1000, kmax=5, b=4.38, c=5.71, d=10.05,
                gamma=1/((MAP/365/k)/1000)*nZ){
  
  gsf1 <- function(w, w0)h2/(h*VPD)*pe*kmax*(w^(-b)-w0^(-b))*exp(-(-pe*w0^(-b)/d)^c)
  gsf2 <- Vectorize(function(w, w0)ifelse(w>w0, gsf1(w, w0), 0))
  gsf <- function(w)gsf2(w, w0=w0)
  
  Ev <- function(w){h*VPD*gsf(w)}
  rEv <- function(w){1/Ev(w)}
  integralrEv <- Vectorize(function(w){integrate(rEv, w, 1, rel.tol = .Machine$double.eps^0.5)$value})
  fnoc <- function(w){rEv(w)*exp(-gamma*w-k*integralrEv(w))}
  integralfnocpdf <- integrate(fnoc, w0, 1, rel.tol = .Machine$double.eps^0.5)$value
  cpdf <- 1/(integralfnocpdf+1/k)
  return(cpdf/k)
  #f <- function(w)cpdf*fnoc(w)
  #return(f(w0+1e-11))
}

# averA
averAf <- function(w0,
                   a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, h2=h/1000, kmax=5, b=4.38, c=5.71, d=10.05,
                   Vcmax=50, cp=30, Km=703, Rd=1,
                   gamma=1/((MAP/365/k)/1000)*nZ){
  
  gsf1 <- function(w, w0)h2/(h*VPD)*pe*kmax*(w^(-b)-w0^(-b))*exp(-(-pe*w0^(-b)/d)^c)
  gsf2 <- Vectorize(function(w, w0)ifelse(w>w0, gsf1(w, w0), 0))
  gsf <- function(w)gsf2(w, w0=w0)
  
  Ev <- function(w){h*VPD*gsf(w)}
  rEv <- function(w){1/Ev(w)}
  integralrEv <- Vectorize(function(w){integrate(rEv, w, 1, rel.tol = .Machine$double.eps^0.5)$value})
  fnoc <- function(w){rEv(w)*exp(-gamma*w-k*integralrEv(w))}
  integralfnocpdf <- integrate(fnoc, w0, 1, rel.tol = .Machine$double.eps^0.5)$value
  cpdf <- 1/(integralfnocpdf+1/k)
  f <- function(w)cpdf*fnoc(w)
  
  Af <- function(w)LAI*1/2*(Vcmax+(Km+ca)*gsf(w)-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gsf(w)+((ca+Km)*gsf(w)+Rd)^2-2*Rd*Vcmax)^(1/2))
  Aff <- function(w)Af(w)*f(w)
  averA <- integrate(Aff, w0, 1)$value#+cpdf/k*Af(w0)
  return(averA)
}

# averwf
averwf <- function(w0,
                   a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, h2=h/1000, kmax=5, b=4.38, c=5.71, d=10.05,
                   gamma=1/((MAP/365/k)/1000)*nZ){
  
  gsf1 <- function(w, w0)h2/(h*VPD)*pe*kmax*(w^(-b)-w0^(-b))*exp(-(-pe*w0^(-b)/d)^c)
  gsf2 <- Vectorize(function(w, w0)ifelse(w>w0, gsf1(w, w0), 0))
  gsf <- function(w)gsf2(w, w0=w0)
  
  Ev <- function(w){h*VPD*gsf(w)}
  rEv <- function(w){1/Ev(w)}
  integralrEv <- Vectorize(function(w){integrate(rEv, w, 1, rel.tol = .Machine$double.eps^0.5)$value})
  fnoc <- function(w){rEv(w)*exp(-gamma*w-k*integralrEv(w))}
  integralfnocpdf <- integrate(fnoc, w0, 1, rel.tol = .Machine$double.eps^0.5)$value
  cpdf <- 1/(integralfnocpdf+1/k)
  f <- function(w)cpdf*fnoc(w)
  
  wf <- function(w)w*f(w)
  res <- integrate(wf, w0, 1, rel.tol = .Machine$double.eps^0.5)$value+cpdf/k*w0
  return(res)
}

# Result
ca <- 400
k <- 0.05
MAP <- 1825
w0 <- 0.20007244034968721
f0f(w0)
averAf(w0)
averwf(w0)
