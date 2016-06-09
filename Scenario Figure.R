
# PLC(px)
PLCf <- function(px, kxmax=5, c=5.71, d=10.05){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  kx <- kxf(px)
  res <- (1-kx/kxmax)*100
  return(res)
}

# Figure
windows(18, 6)
par(mgp=c(2.2, 1, 0), lwd=2, mar=c(1.5, 2, 1.5, 0), mfrow=c(1,3))
# Scenario 1
curve(x^2, -15, 0, type="n",
      xaxt="n", yaxt="n",
      xlim=c(-15, 0), ylim=c(0, 100),
      main="Scenario 1",
      xlab=NA, ylab=NA,
      cex.main=2.5, bty="n")
axis(1, xlim=c(-15, 0), at=c(-15, 0), pos=0, lwd=2, cex.axis=1.5)
mtext(expression(psi[x]),side=1,line=0.3, cex=1.5)
axis(2, xlim=c(0, 100), at=c(0, 100), pos=-15, lwd=2, cex.axis=1.5)
mtext("PLC (%)",side=2,line=0.2, cex=1.5)
curve(PLCf, -15, -1e-3, add=T)
points(-8.5, PLCf(-8.5), pch=16, cex=2)
# Scenario 2
curve(x^2, -15, 0, type="n",
      xaxt="n", yaxt="n",
      xlim=c(-15, 0), ylim=c(0, 100),
      main="Scenario 2",
      xlab=NA, ylab=NA,
      cex.main=2.5, bty="n")
axis(1, xlim=c(-15, 0), at=c(-15, 0), pos=0, lwd=2, cex.axis=1.5)
mtext(expression(psi[x]),side=1,line=0.3, cex=1.5)
axis(2, xlim=c(0, 100), at=c(0, 100), pos=-15, lwd=2, cex.axis=1.5)
mtext("PLC (%)",side=2,line=0.2, cex=1.5)
curve(PLCf, -15, -8.5, add=T)
curve(PLCf, -8.5, -1e-3, lty=2, add=T)
curve(PLCf(-8.5)*x^0, -8.5, -1e-3, add=T)
points(-8.5, PLCf(-8.5), pch=16, cex=2)
# Scenario 3
curve(x^2, -15, 0, type="n",
      xaxt="n", yaxt="n",
      xlim=c(-15, 0), ylim=c(0, 100),
      main="Scenario 3",
      xlab=NA, ylab=NA,
      cex.main=2.5, bty="n")
axis(1, xlim=c(-15, 0), at=c(-15, 0), pos=0, lwd=2, cex.axis=1.5)
mtext(expression(psi[x]),side=1,line=0.3, cex=1.5)
axis(2, xlim=c(0, 100), at=c(0, 100), pos=-15, lwd=2, cex.axis=1.5)
mtext("PLC (%)",side=2,line=0.2, cex=1.5)
curve(PLCf, -15, -8.5, add=T)
curve(PLCf, -8.5, -1e-3, lty=2, add=T)
curve(PLCf(-8.5)*x^0, -8.5, -1e-3, lty=2, add=T)
curve(PLCf(-8.5)-(PLCf(-8.5)-PLCf(x))/2.3, -8.5, -1e-3, add=T)
points(-8.5, PLCf(-8.5), pch=16, cex=2)
