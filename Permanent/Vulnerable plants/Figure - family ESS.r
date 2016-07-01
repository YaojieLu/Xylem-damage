
################################ Permanent ################################
############################ Vulnerable plants ############################
source("Permanent/Vulnerable plants/Functions.r")

# Initializing
pxminf1 <- Vectorize(pxminf)
gsmaxf1 <- Vectorize(gsmaxf)

f <- function(w)psf(w)-pxminf1(1)
w1 <- uniroot(f, c(0.1, 1), tol=.Machine$double.eps) #0.18686045477181246

# Figure
Cols <- rainbow(4)
wL <- c(0.24, 0.23, 0.22, 0.207)
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1), mfrow=c(1, 1))
with(list(wL=wL[1], Col=Cols[1]),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1,
           xlim=c(0, 1), ylim=c(0, 0.25),
           xlab=expression(italic(w)), ylab=NA,
           cex.lab=1.3, col=Col)})
with(list(wL=wL[2], Col=Cols[2]),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col=Col)})
with(list(wL=wL[3], Col=Cols[3]),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col=Col)})
with(list(wL=wL[4], Col=Cols[4]),
     {gswLf1 <- function(w)gswLf(w, wL)
     curve(gswLf1, wL, 1, add=T, col=Col)})

mtext(expression(italic(g[s])~(mol~m^-2~s^-1)),side=2,line=1.7, cex=1.3)
legend("topleft", as.character(wL), bg="white", lty=c(1, 1, 1, 1), col=Cols)
