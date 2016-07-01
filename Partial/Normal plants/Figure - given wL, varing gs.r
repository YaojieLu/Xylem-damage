
############################# Partial #############################
############################# Normal #############################
source("Partial/Normal plants/Functions.r")

# Initializing
ca <- 400
pkx <- 0.9
wL <- 0.18
wLL <- wLLf(wL, ca)

inversegsmaxfm <- function(gs)optimize(function(w)(gsmaxfm(w, wL)-gs)^2, c(wLL$root, 1), tol=.Machine$double.eps)$minimum
Bgsmaxf <- Vectorize(function(gs)Bfm(inversegsmaxfm(gs), gs, wL, ca))
inversegswLf <- function(gs)optimize(function(w)(gswLf(w, wL, ca)-gs)^2, c(wLL$root, 1), tol=.Machine$double.eps)$minimum
BgswLf <- Vectorize(function(gs)Bfm(inversegswLf(gs), gs, wL, ca))

# Figure
w <- c(1, 0.2, wLL$root)
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="r", yaxs="r", lwd=2, mar=c(3.5, 3.5, 0.5, 0.5), mfrow=c(1, 1))
with(data.frame(w=w[1]),
     {
       Bgsf <- Vectorize(function(gs)Bfm(w, gs, wL, ca))
       maxgs <- gsmaxfm(w, wL)
       ESSgs <- gswLf(w, wL, ca)
       curve(Bgsf, 0, maxgs,
             xlab=NA, ylab=NA, xaxt="n",
             xlim=c(0, 0.3), ylim=c(-10, 15), cex.lab=1.5, lwd=1)
       text(maxgs+0.015, Bgsf(maxgs)-0.5, as.expression(bquote(italic(w)==.(w))), cex=1.3)
     }
)
with(data.frame(w=w[2]),
     {
       Bgsf <- Vectorize(function(gs)Bfm(w, gs, wL, ca))
       maxgs <- gsmaxfm(w, wL)
       curve(Bgsf, 0, maxgs, lwd=1, add=T)
       text(maxgs+0.0175, Bgsf(maxgs)-0.7, as.expression(bquote(italic(w)==.(w))), cex=1.3)
     }
)
with(data.frame(w=w[3]),
     {
       Bgsf <- Vectorize(function(gs)Bfm(w, gs, wL, ca))
       maxgs <- gsmaxfm(w, wL)
       curve(Bgsf, 0, maxgs, lwd=1, add=T)
       text(maxgs+0.02, Bgsf(maxgs)-0.7, expression(italic(w==w[LL])), cex=1.3)
     }
)

# B(gs) at given w with gs is indicated by gs=0, gsmaxf and ESS
segments(0, Bfm(wLL$root, 0, wL, ca), 0, Bfm(1, 0, wL, ca), col="red")
curve(Bgsmaxf, gsmaxfm(wLL$root, wL), gsmaxfm(1, wL), col="orange", add=T)
segments(-0.5*0.04, 0, gswLf(wLL$root, wL, ca), 0, lwd=1, lty=3)
curve(BgswLf, gswLf(wLL$root, wL, ca), gswLf(1, wL, ca), col="blue", lty=2, add=T)

# legend and labels
axis(1, pos=-10-25*0.04, lwd=2, at=c(0, 0.1, 0.2, 0.3))
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)), side=1, line=2.6, cex=1.3)
mtext(expression(italic(B)~(mu*mol~m^-2~s^-1)), side=2, line=1.85, cex=1.3)
legend("bottomright", expression(family~italic(g[smax]), italic(g[s]==0), family~ESS~italic(g[s])),
       col = c("orange","red", "blue"), lty=c(1, 1, 2), lwd=c(2, 2, 2))
box()

dev.copy2pdf(file = "Figures/Figure 9.pdf")
