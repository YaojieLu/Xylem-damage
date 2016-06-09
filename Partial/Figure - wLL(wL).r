
source("Partial/gswL(w).r")

pkx <- 0.5
ca <- 400
f <- Vectorize(function(wL)wLLf(wL)$root)

# Figure
# Low - wLL increases with decreasing wL
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
curve(f, 0.12, 0.13, xlim=c(0.12, 0.13), ylim=c(0.135, 0.14), cex.lab=1.3,
      main="Partial xylem damage - family ESS",
      xlab=expression(italic(w[L])),
      ylab=expression(italic(w[LL])))

# Medium
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
pkx <- 0.5
ca <- 400
curve(f, 0.13, 0.14, xlim=c(0.13, 0.14), ylim=c(0.13, 0.16), cex.lab=1.3,
      main="Partial xylem damage - family ESS",
      xlab=expression(italic(w[L])),
      ylab=expression(italic(w[LL])))

# High - wLL increases with increasing wL
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
pkx <- 0.5
ca <- 400
curve(f, 0.14, 0.16, xlim=c(0.14, 0.16), ylim=c(0.14, 0.16), cex.lab=1.3,
      main="Partial xylem damage - family ESS",
      xlab=expression(italic(w[L])),
      ylab=expression(italic(w[LL])))
