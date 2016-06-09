
# the original PLC(px)
PLCf <- function(px, c=5.71, d=10.05)1-exp(-(-px/d)^c)

# modified gsmax
gsmaxfm <- function(w, wL,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05){
  
  pxL <- pe*wL^(-b)
  # modified PLC
  PLCfm <- function(x)PLCf(pxL)-(PLCf(pxL)-PLCf(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  ps <- pe*w^(-b)
  f1 <- function(x)(ps-x)*h2*kxfm(x)/(h*VPD)
  
  res <- ifelse(pxL<ps, optimize(f1, c(pxL,ps), tol=.Machine$double.eps, maximum=T)$objective, 0)
  return(res)
}

# gsmax
gsmaxf <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05){
  
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  ps <- pe*w^(-b)
  f1 <- function(x)(ps-x)*h2*kxf(x)/(h*VPD)
  res <- optimize(f1, c(-20,pe), tol=.Machine$double.eps, maximum=T)
  return(res$objective)
}
gsmaxf1 <- Vectorize(gsmaxf)

#pkx <- 0.5
#wL <- 0.13
#gsmaxfm1 <- Vectorize(function(w)gsmaxfm(w, wL))
#curve(gsmaxf1, wL, wL+0.01, xlim=c(wL, wL+0.01))
#curve(gsmaxfm1, wL, wL+0.01, add=T, col="red")
#

# Figure
pkx <- 0.5
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
with(data.frame(wL=c(0.13)),
     {gsmaxfm1 <- Vectorize(function(w)gsmaxfm(w, wL))
       curve(gsmaxf1, wL, 1, xlim=c(0, 1), ylim=c(0, 2))
       curve(gsmaxfm1, wL, 1, add=T, col="red")})
