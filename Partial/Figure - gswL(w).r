
source("Partial/gswL(w).r")

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
pkx <- 0.5
ca <- 400

with(data.frame(wL=c(0.3)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL))
     wLL <- wLLf(wL)$root
     curve(gswLf1, wLL, 1, xlim=c(0, 1), ylim=c(0, 2), cex.lab=1.3,
           main="Partial xylem damage - family ESS",
           xlab=expression(italic(w)),
           ylab=expression(italic(g[s])~(mol~m^-2~s^-1)))})
with(data.frame(wL=c(0.161)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL))
     wLL <- wLLf(wL)$root
     curve(gswLf1, wLL, 1, add=T, col="blue")})
with(data.frame(wL=c(0.16)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL))
     wLL <- wLLf(wL)$root
     curve(gswLf1, wLL, 1, add=T, col="darkgreen")})
with(data.frame(wL=c(0.13)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL))
     wLL <- wLLf(wL)$root
     curve(gswLf1, wLL, 1, add=T, col="orange")})
with(data.frame(wL=c(0.12)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL))
     wLL <- wLLf(wL)$root
     curve(gswLf1, wLL, 1, add=T, col="red")})
legend("topleft", c("0.3", "0.161", "0.16", "0.13", "0.12"), lty=c(1, 1, 1, 1, 1), col=c("black", "blue", "darkgreen", "orange", "red"))
