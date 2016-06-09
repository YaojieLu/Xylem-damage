
# gsmax
gsmaxf <- function(w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)(ps-x)*h2*kxf(x)
    res <- optimize(f1, c(-20000000,pe), tol=.Machine$double.eps, maximum=T)$maximum
    return(res)
  }
  
  ps <- pe*w^(-b)
  pxmin <- pxminf(w)
  res <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}
gsmaxf1 <- Vectorize(gsmaxf)

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

# Figure
# gswL(w)
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
curve(gsmaxf1, 0, 1, xlim=c(0, 1), ylim=c(0, 0.5), cex.lab=1.3,
      main="Permanent xylem damage - family ESS",
      xlab=expression(italic(w)),
      ylab=expression(italic(g[s])~(mol~m^-2~s^-1)))
with(data.frame(wL=c(0.2)),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col="darkgreen")})
with(data.frame(wL=c(0.19)),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col="blue")})
with(data.frame(wL=c(0.18)),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col="red")})
with(data.frame(wL=c(0.15)),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col="orange")})
legend("topleft", c("0.2", "0.19", "0.18", "0.14"), lty=c(1, 1, 1, 1), col=c("darkgreen", "blue", "red", "orange"), title=expression(italic(w[L])))
