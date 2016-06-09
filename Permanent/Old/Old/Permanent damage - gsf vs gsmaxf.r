
gsf1 <- function(w, w0,
                 a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 pe=-1.58*10^-3, h2=h/1000, kmax=5, b=4.38, c=5.71, d=10.05)h2/(h*VPD)*pe*kmax*(w^(-b)-w0^(-b))*exp(-(-pe*w0^(-b)/d)^c)
gsf2 <- Vectorize(function(w, w0)ifelse(w>w0, gsf1(w, w0), 0))
gsf3 <- function(w)gsf2(w, w0=w0)

gsmaxf1 <- function(w,
                    a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, h2=h/1000, kmax=5, b=4.38, c=5.71, d=10.05){
  ps <- pe*w^(-b)
  kf <- function(x)kmax*exp(-(-x/d)^c)
  f1 <- function(x)-((ps-x)*h2*kf(x))
  pxmin <- optimize(f1, c(-20,0), tol = .Machine$double.eps)$minimum
  gsmax <- (ps-pxmin)*h2*kf(pxmin)/(h*VPD)
  return(gsmax)
}
gsmaxf2 <- Vectorize(gsmaxf1)

f <- function(w0)gsmaxf2(w0+1e-5)-gsf2(w0+1e-5, w0)
# Figure
w0 <- 0.61439263589695858
curve(gsf3, 0, 1)
curve(gsmaxf2, 0, 1, col="red")
