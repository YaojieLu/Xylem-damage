
options (digits=20)

# Function
ESSw0Rf1 <- function(ca, k, MAP,
                     a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                     pe=-1.58*10^-3, h2=h/1000, kmax=5, b=4.38, c=5.71, d=10.05,
                     Vcmax=50, cp=30, Km=703, Rd=1,
                     gamma=1/((MAP/365/k)/1000)*nZ){
  
  w0Rf1 <- function(w0R){
    
    w0If1 <- function(w0I){
      
      gsf1 <- function(w, w0)h2/(h*VPD)*pe*kmax*(w^(-b)-w0^(-b))*exp(-(-pe*w0^(-b)/d)^c)
      gsf2 <- Vectorize(function(w, w0)ifelse(w>w0, gsf1(w, w0), 0))
      gsRf <- function(w)gsf2(w, w0=w0R)
      gsIf <- function(w)gsf2(w, w0=w0I)
      
      Ev <- function(w){h*VPD*gsRf(w)}
      rEv <- function(w){1/Ev(w)}
      integralrEv <- Vectorize(function(w){integrate(rEv, w, 1, rel.tol = .Machine$double.eps^0.5)$value})
      fnoc <- function(w){rEv(w)*exp(-gamma*w-k*integralrEv(w))}
      integralfnocpdf <- integrate(fnoc, w0R, 1, rel.tol = .Machine$double.eps^0.5)$value
      cpdf <- 1/(integralfnocpdf+1/k)
      f <- function(w)cpdf*fnoc(w)
      
      AIf <- function(w)LAI*1/2*(Vcmax+(Km+ca)*gsIf(w)-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gsIf(w)+((ca+Km)*gsIf(w)+Rd)^2-2*Rd*Vcmax)^(1/2))
      AIff <- function(w)AIf(w)*f(w)
      averAI <- integrate(AIff, w0R, 1)$value+cpdf/k*AIf(w0R)
      return(averAI)
    }
    w0If2 <- function(w0I)-w0If1(w0I)
    w0If3 <- Vectorize(w0If2)
    optw0I <- optimize(w0If3, c(0.2, 0.999))$minimum
    res <- optw0I-w0R
    return(res^2)
  }
  w0Rf2 <- Vectorize(w0Rf1)
  res <- optimize(w0Rf2, c(0.2, 0.999))
  return(res)
}
ESSw0Rf2 <- Vectorize(ESSw0Rf1)

# Result
res1 <- ESSw0Rf2(400, 0.025, 365)#0.20007244034968721 -10
res2 <- ESSw0Rf2(400, 0.025, 3650)#0.20007244034968721 -10
res3 <- ESSw0Rf2(400, 0.05, 1825)#0.20007244034968721 -10
res4 <- ESSw0Rf2(400, 0.1, 365)#0.20004285115680265 -12
res5 <- ESSw0Rf2(400, 0.1, 3650)#0.20007244034968721 -10