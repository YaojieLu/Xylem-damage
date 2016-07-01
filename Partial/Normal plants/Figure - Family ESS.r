
############################# Partial #############################
############################# Normal #############################
source("Partial/Normal plants/Functions.r")

# Initializing
ca <- 400
pkx <- 0.9
wL <- 0.18
wLL <- wLLf(wL, ca)

# Figure
wL <- c(0.22, 0.21, 0.20, 0.15, 0.1)
Cols <- rainbow(5)
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="r", lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(1, 1))
# Family ESS
with(list(wL=wL[1], Col=Cols[1]),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, xlim=c(0, 1), ylim=c(0, 0.3), cex.lab=1.5,
           xlab=expression(italic(w)),
           ylab=expression(italic(g[s])~(mol~m^-2~s^-1)), col=Col)
     segments(0, 0, wLL, 0, col=Col, lty=2)
     })
with(list(wL=wL[2], Col=Cols[2]),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, add=T, col=Col)
     segments(0, 0, wLL, 0, col=Col, lty=2)
     })
with(list(wL=wL[3], Col=Cols[3]),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, add=T, col=Col)
     segments(0, 0, wLL, 0, col=Col, lty=2)
     })
with(list(wL=wL[4], Col=Cols[4]),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, add=T, col=Col)
     segments(0, 0, wLL, 0, col=Col, lty=2)
     })
with(list(wL=wL[5], Col=Cols[5]),
     {gswLf1 <- Vectorize(function(w)gswLf(w, wL, ca))
     wLL <- wLLf(wL, ca)$root
     curve(gswLf1, wLL, 1, add=T, col=Col)
     segments(0, 0, wLL, 0, col=Col, lty=2)
     })

legend("bottomright", title=expression(italic(w[L])), as.character(wL), lty=c(1, 1, 1, 1, 1), col=Cols)
box()

dev.copy2pdf(file = "Figures/Figure 10.pdf")
