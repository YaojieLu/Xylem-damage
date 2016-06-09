
# gsmax
gswLwf <- function(wL, w,
                   a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                   pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  px <- pe*wL^(-b)
  ps <- pe*w^(-b)
  res <- (ps-px)*h2*kxf(px)/(h*VPD)
  return(res)
}


# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1,1))
with(data.frame(w=c(1)),
     {
       gswLf <- Vectorize(function(wL)gswLwf(wL, w))
       label <- optimize(gswLf, c(0, 1), maximum=T)$maximum
       curve(gswLf, 0, 1, n=10001,
             xlab=expression(italic(w[L])),
             ylab=expression(italic(g[s])~(mol~m^-2~s^-1)),
             xlim=c(0, 1), ylim=c(0, 2), cex.lab=1.3, col="black")
       abline(v=label)
       text(x=0.27, y=1.7, labels=c(eval(substitute(expression(italic(w[L])==x), list(x=label)))), cex=1.3)
     }
)
with(data.frame(w=c(0.2)),
     {
       gswLf <- Vectorize(function(wL)gswLwf(wL, w))
       curve(gswLf, 0, 1, n=10001, col="red", add=T)
     }
)
with(data.frame(w=c(0.15)),
     {
       gswLf <- Vectorize(function(wL)gswLwf(wL, w))
       curve(gswLf, 0, n=10001, 1, col="blue", add=T)
     }
)
legend("topright", c("1", "0.2", "0.15"), lty=c(1, 1, 1), col=c("black", "red", "blue"), title=expression(italic(w)))
