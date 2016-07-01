
############################ Normal plants ############################
############################## Permanent ##############################
source("Permanent/Normal plants/Functions.r")

# Initializing
pxminf1 <- Vectorize(pxminf)
gsmaxf1 <- Vectorize(gsmaxf)

f <- function(w)psf(w)-pxminf1(1)
w1 <- uniroot(f, c(0.1, 1), tol=.Machine$double.eps) #0.18686045477181246

# Figure
Cols <- rainbow(4)
windows(8, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1), mfrow=c(2, 1))
# family ESS
with(list(wL=c(0.2), Col=Cols[1]),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1,
           xlim=c(0, 1), ylim=c(0, 0.3), yaxt="n",
           xlab=expression(italic(w)), ylab=NA,
           cex.lab=1.3, col=Col)})
with(list(wL=c(0.19), Col=Cols[2]),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col=Col)})
with(list(wL=c(0.18), Col=Cols[3]),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col=Col)})
with(list(wL=c(0.15), Col=Cols[4]),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col=Col)})

axis(2, pos=0, lwd=2, at=c(0, 0.1, 0.2, 0.3))
mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.9, cex=1.3)
legend("topleft", title=expression(italic(w[L])), c("0.2", "0.19", "0.18", "0.15"), bg="white", lty=c(1, 1, 1, 1), col=Cols)
text(0.95, 0.3*0.95, "a", cex=1.5)
box()

# family ESS as px(ps)
with(list(wL=c(0.2), Col=Cols[1]),
     {pxpswLf1 <- function(ps)pxpswLf(ps, wL)
     curve(pxpswLf1, psf(wL), 1,
           xlim=c(-8, 0), ylim=c(-8, 0), xlab=NA, ylab=NA,
           cex.lab=1.3, col=Col)})
with(list(wL=c(0.19), Col=Cols[2]),
     {pxpswLf1 <- function(ps)pxpswLf(ps, wL)
     curve(pxpswLf1, psf(wL), 1, add=T, col=Col)})
with(list(wL=c(0.18), Col=Cols[3]),
     {pxpswLf1 <- function(ps)pxpswLf(ps, wL)
     curve(pxpswLf1, psf(wL), 1, add=T, col=Col)})
with(list(wL=c(0.15), Col=Cols[4]),
     {pxpswLf1 <- function(ps)pxpswLf(ps, wL)
     curve(pxpswLf1, psf(wL), 1, add=T, col=Col)})
abline(a=0, b=1, lty=2, lwd=1)

mtext(expression(italic(psi[s])~(MPa)),side=1,line=2.3, cex=1.3)
mtext(expression(italic(psi[x])~(MPa)),side=2,line=1.7, cex=1.3)
text(-8*0.95, -8*0.05, "b", cex=1.5)
box()

dev.copy2pdf(file = "Figures/Figure 7.pdf")

## water use envelope and soil water retention curve
#curve(pxminf1, 0.1, 1, xlim=c(0, 1), ylim=c(-8, 0), cex.lab=1.3,
#      xlab=expression(italic(w)),
#      ylab=expression(italic(psi)~(MPa)), col="darkgreen")
#curve(psf, 0, 1, add=T, col="red")
#arrows(1, pxminf1(1), w1$root, pxminf1(1), lty=2, lwd=1, col="darkgreen")
#arrows(w1$root, pxminf1(1), w1$root, -8, lty=2, lwd=1, col="red")
#legend("bottomright", expression(italic(psi[xmin](w)), italic(psi[xmin](1)), italic(psi[s])(w), italic(w[1])), lwd=c(2, 1, 2, 1), lty=c(1, 2, 1, 2), col=c("darkgreen", "darkgreen", "red", "red"))
#text(0.05, -8*0.06, "b", cex=1.5)
#
## family ESS as px(w)
#with(data.frame(wL=c(0.2)),
#     {pxwLf1 <- function(w)pxwLf(w, wL)
#     curve(pxwLf1, wL, 1,
#           xlim=c(0, 1), ylim=c(-8, 0),
#           xlab=expression(italic(w)), ylab=expression(italic(psi[x])~(MPa)),
#           cex.lab=1.3, col="darkgreen")
#     segments(wL, -8, wL, pxwLf1(wL), lty=2, lwd=1, col="darkgreen")
#     })
#with(data.frame(wL=c(0.19)),
#     {pxwLf1 <- function(w)pxwLf(w, wL)
#     curve(pxwLf1, wL, 1, add=T, col="blue")
#     segments(wL, -8, wL, pxwLf1(wL), lty=2, lwd=1, col="blue")
#     })
#with(data.frame(wL=c(0.18)),
#     {pxwLf1 <- function(w)pxwLf(w, wL)
#     curve(pxwLf1, wL, 1, add=T, col="red")
#     segments(wL, -8, wL, pxwLf1(wL), lty=2, lwd=1, col="red")
#     })
#with(data.frame(wL=c(0.15)),
#     {pxwLf1 <- function(w)pxwLf(w, wL)
#     curve(pxwLf1, wL, 1, add=T, col="orange")
#     segments(wL, -8, wL, pxwLf1(wL), lty=2, lwd=1, col="orange")
#     })
#
#text(0.95, -8*0.06, "c", cex=1.5)
#
