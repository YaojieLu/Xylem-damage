
options(digits=4)
# A
Af <- function(gs, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# gswL(w)
gswLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54){
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  pxmin <- pe*wL^(-b)
  res <- (pe*w^-b-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}

# averBi
averBif <- function(wLi, wLr,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54,
                    gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswLfr <- function(w)ifelse(w<wLr+1e-10, 0, gswLf(w, wLr))
  gswLfi <- function(w)ifelse(w<wLi+1e-10, 0, gswLf(w, wLi))
  
  Ef <- function(w){h*VPD*gswLfr(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.4)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLr-1e-12)/(1-wLr-1e-12)-k*integralrEf(w)*1/(1-wLr-1e-12))*1/(1-wLr-1e-12)
  f1 <- function(w)Af(gswLfi(w))*fnoc(w)
  fLnoc <- 1/k*exp(-k*integrate(rEf, wLr+1e-12, 1, rel.tol=.Machine$double.eps^0.4)$value*1/(1-wLr-1e-12))
  res1 <- integrate(f1, wLr+1e-12, 1, rel.tol=.Machine$double.eps^0.4)
  res2 <- fLnoc*Af(gswLfi(wLr+1e-12))
  res <- res1$value+res2
  return(res)
}

# Results
ca <- 400
k <- 0.05
MAP <- 1000

wLr <- 0.1627
f <- Vectorize(function(wLi)averBif(wLi, wLr)*100)
x <- seq(0.16, 0.17, by=0.0001)
y <- f(x)

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 3, 2, 2), mfrow=c(1,1))
plot(x, y, xlim=c(0.16, 0.17), ylim=c(140, 145),
     xlab=expression(italic(w[Li])), ylab=NA, yaxt="n",
     main=as.expression(bquote(italic(w[Lr])==.(wLr))),
     type="l", cex.lab=1.3)
points(wLr, f(wLr), pch=16, cex=1.5)
mtext(expression("Relative"~bar(italic(B))),side=2,line=1, cex=1.3)
