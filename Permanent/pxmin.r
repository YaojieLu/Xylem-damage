
# ps
psf <- function(w, pe=-1.58*10^-3, b=4.38)pe*w^(-b)

# pxmin
pxminf <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, h2=h/1000, kxmax=5, c=5.71, d=10.05){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  ps <- psf(w)
  f1 <- function(x)(ps-x)*h2*kxf(x)
  res <- optimize(f1, c(-20000000,pe), tol=.Machine$double.eps, maximum=T)$maximum
  return(res)
}

# Figure
# gswL(w)
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1,1))
pxminf1 <- Vectorize(pxminf)
curve(pxminf1, 0.11, 1, xlim=c(0, 1), ylim=c(-20, 0), cex.lab=1.3,
      xlab=expression(italic(w)),
      ylab=expression(psi[x]))
curve(psf, 0.11, 1, add=T, col="red")