
############################# Partial #############################
############################# Normal #############################
source("Partial/Normal plants/Functions.r")

# Initializing
ca <- 400
pkx <- 0.9

# Figure
windows(24, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4.5, 1, 0.5), mfrow=c(1, 3))
# Family ESS
with(data.frame(wL=c(0.15)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, xlim=c(0, 1), ylim=c(0, 0.3), cex.lab=1.5,
           xlab=expression(italic(w)),
           ylab=expression(italic(g[s])~(mol~m^-2~s^-1)))})
with(data.frame(wL=c(0.16)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, add=T, col="blue")})
with(data.frame(wL=c(0.21)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, add=T, col="darkgreen")})
with(data.frame(wL=c(0.22)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, add=T, col="orange")})
with(data.frame(wL=c(0.23)),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, add=T, col="red")})
legend("bottomright", title=expression(italic(w[L])*"="), c("0.14", "0.15", "0.21", "0.22", "0.23"), lty=c(1, 1, 1, 1, 1), col=c("black", "blue", "darkgreen", "orange", "red"))
text(0.05, 0.3*0.95, "a", cex=1.5)

# Soil water retention curve
curve(PLCwf, 0, 1, xlim=c(0, 1), ylim=c(0, 100), cex.lab=1.5,
      xlab=expression(italic(w)),
      ylab="PLC (%)")
text(0.95, 95, "b", cex=1.5)
