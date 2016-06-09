
options(digits=20)
# A
Af <- function(gs, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# gswL(w)
gswLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=9.53, d=1.28){
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  pxmin <- pe*wL^(-b)
  res <- (pe*w^-b-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}

# averBi
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=9.53, d=1.28,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswLfr <- function(w)ifelse(w<wLr+1e-10, 0, gswLf(w, wLr))
  gswLfi <- function(w)ifelse(w<wLi+1e-10, 0, gswLf(w, wLi))
  
  Ef <- function(w){h*VPD*gswLfr(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.25)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLr-1e-10)/(1-wLr-1e-10)-k*integralrEf(w)*1/(1-wLr-1e-10))*1/(1-wLr-1e-10)
  fLnoc <- 1/k*exp(-k*integrate(rEf, wLr+1e-10, 1, rel.tol=.Machine$double.eps^0.25)$value*1/(1-wLr-1e-10))
  cPDF <- 1/(fLnoc+integrate(fnoc, wLr+1e-10, 1, rel.tol=.Machine$double.eps^0.25)$value)
  
  f1 <- function(w)cPDF*fnoc(w)
  f2 <- function(w)Af(gswLfi(w))*f1(w)
  fL <- cPDF*fLnoc
  res1 <- integrate(f2, wLr+1e-10, 1, rel.tol=.Machine$double.eps^0.25)
  res2 <- fL*Af(gswLfi(wLr+1e-10))
  res <- res1$value+res2
  return(res)
}

# Auxiliary function
# Optimization among wLi < wLr
f1 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(0.1, wLr), tol=.Machine$double.eps, maximum=T)
  return(res$maximum)
}
# Optimization among wLi > wLr
f2 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(wLr, 0.3), tol=.Machine$double.eps, maximum=T)
  return(res$maximum)
}
f11 <- Vectorize(f1)
f21 <- Vectorize(f2)

# Results
ca <- 400
k <- 0.1
MAP <- 300

x <- seq(0.207, 0.23, by=0.0001)
resL <- f11(x)
resH <- f21(x)

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(4, 4, 3, 2), mfrow=c(1,1))
plot(x, resL, xlim=c(0.207, 0.23), lwd=2, ylim=c(0.207, 0.23), type="l",
     xlab=expression(italic(w[Lr])), ylab=expression(italic(w[Li])),
     main="c=9.53, d=1.28, k=0.1, MAP=300",
     cex.lab=1.3, col="blue")
points(x, resH, type="l", cex.lab=1.3, col="red", lwd=2)
abline(a=0, b=1)
