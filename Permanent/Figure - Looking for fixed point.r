
############################# Permanent #############################
options(digits=20)
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
  
  gswLfr <- function(w)ifelse(w<wLr+1e-10, 0, gswLf(w, wLr))
  gswLfi <- function(w)ifelse(w<wLi+1e-10, 0, gswLf(w, wLi))
  
  Ef <- function(w){h*VPD*gswLfr(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.25)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLr-1e-10)/(1-wLr-1e-10)-k*integralrEf(w)*1/(1-wLr-1e-10))*1/(1-wLr-1e-10)
  f1 <- function(w)Af(gswLfi(w))*fnoc(w)
  fLnoc <- 1/k*exp(-k*integrate(rEf, wLr+1e-10, 1, rel.tol=.Machine$double.eps^0.25)$value*1/(1-wLr-1e-10))
  res1 <- integrate(f1, wLr+1e-10, 1, rel.tol=.Machine$double.eps^0.25)
  res2 <- fLnoc*Af(gswLfi(wLr+1e-10))
  res <- res1$value+res2
  #browser()
  return(res)
}

f1 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  resL <- optimize(averBif1, c(0.1, wLr), tol=.Machine$double.eps, maximum=T)
  resH <- optimize(averBif1, c(wLr, 0.2), tol=.Machine$double.eps, maximum=T)
  res <- ifelse(resL$objective>resH$objective, resL$maximum, resH$maximum)
  return(res)
}
f2 <- Vectorize(f1)

# Results
ca <- 400
k <- 0.025
MAP <- 365

x <- seq(0.12, 0.145, by=0.0001)
res <- f2(x)

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(4, 4, 2, 2), mfrow=c(1,1))
plot(x, res, xlim=c(0.12, 0.145), ylim=c(0.12, 0.15), type="p",
     xlab=expression(italic(w[Lr])), ylab=expression(Best~italic(w[Li])),
     cex.lab=1.3, col="blue", lwd=2)
abline(a=0, b=1)
