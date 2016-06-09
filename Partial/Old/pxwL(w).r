
source("Partial/gswL(w).r")

pxwLf <- function(w, wL,
                  a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5){
  
  gswLf1 <- Vectorize(function(w)gswLf(w, wL))
  pxL <- pe*wL^(-b)
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)(ps-x)*h2*kxfm(x)
    res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, ps)
    return(res)
  }
  
  ps <- pe*w^(-b)
  pxmin <- pxminf(w)
  gs <- gswLf1(w)
  f1 <- function(x)((ps-x)*h2*kxfm(x)-h*VPD*gs)^2
  res <- ifelse(pxmin<ps, optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum, ps)
  return(res)
}

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
pkx <- 0.5
ca <- 400

with(data.frame(wL=c(0.13)),
     {pxwLf1 <- Vectorize(function(w)pxwLf(w, wL))
     wLL <- wLLf(wL)$root
     curve(pxwLf1, wLL, 1, xlim=c(0, 1), ylim=c(-15, 0), cex.lab=1.3,
           xlab=expression(italic(w)),
           ylab=expression(italic(psi[x])~(MPa)))
     abline(h=-1.58*10^-3*wL^(-4.38))})
with(data.frame(wL=c(0.1)),
     {pxwLf1 <- Vectorize(function(w)pxwLf(w, wL))
     wLL <- wLLf(wL)$root
     curve(pxwLf1, wLL, 1, add=T, col="red")
     abline(h=-1.58*10^-3*wL^(-4.38), col="red")})
legend("topleft", c("0.13", "0.1"), lty=c(1, 1), col=c("black", "red"))
