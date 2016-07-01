
source("Permanent/Normal plants/Functions.r")

# Initialize
pe <- -1.58*10^-3
psL <- -3.5
pkx <- 0.5

PLCf1 <- function(px)PLCf(px)*100
PLCmax <- PLCf1(psL)
PLCmf <- function(px)PLCmax-(PLCmax-PLCf1(px))*pkx

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", lwd=2, mar=c(3.5, 3.5, 0.5, 0.5), mfrow=c(1, 1))
curve(PLCf1, -10, psL,
      xlim=c(-10, 0), ylim=c(0, 100),
      xlab=expression(psi[x]~(MPa)), ylab="PLC (%)",
      cex.lab=1.3, lty=2)

curve(PLCf1, psL, pe, col="orange", add=T)
segments(psL, PLCmax, pe, PLCmax, col="red")
curve(PLCmf, psL, pe, col="blue", add=T)

arrows(-1, PLCmf(-1)*1.1, -1, PLCmax, lwd=1, col="blue")
arrows(-1, PLCmf(-1)*0.9, -1, PLCf1(-1), lwd=1, col="blue")
arrows(psL, PLCmax, psL, -4, lty=2, lwd=1)
text(psL-0.3, PLCmax/2, expression(psi[xL]), cex=1.3)
arrows(psL, PLCmax, -10, PLCmax, lty=2, lwd=1)
text((-10-psL)/2+psL, PLCmax+3, expression(PLC[italic(max)]), cex=1.3)

legend("topright", c("Scenario 1", "Scenario 2", "Scenario 3"), lty=c(1, 1, 1), col=c("orange", "red", "blue"))
box()

dev.copy2pdf(file = "Figures/Figure 6.pdf")
