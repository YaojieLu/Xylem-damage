
options(digits=3)
# PLC(px)
PLCf <- function(px, c=5.71, d=10.05)1-exp(-(-px/d)^c)
# PLCwL(wL)
PLCwLf <- function(wL, pe=-1.58*10^-3, b=4.38, c=5.71, d=10.05)1-exp(-(-(pe*wL^(-b))/d)^c)

# Figure
# PLC - px
windows(16, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1,2))
curve(PLCf, -15, -1e-3,
      xlim=c(-15, 0), ylim=c(0, 1),
      xlab=expression(psi[x]),
      ylab=expression(PLC[m]),
      cex.lab=1.3, col="red")
# pk=0.5
with(data.frame(pxL=c(-10)),
     {curve(PLCf(pxL)-(PLCf(pxL)-PLCf(x))*0.5, pxL, -1e-3, add=T)
       points(pxL, PLCf(pxL), pch=16, cex=1.2)})
text(-8.7, 0.628, expression(italic(psi[x])==-10), cex=1.3)
text(-2, 0.35, expression(italic(p[k])==0.5), cex=1.3)
text(-0.4, 0.96, "a", cex=1.3)
# pk=0.8
with(data.frame(pxL=c(-10)),
     {curve(PLCf(pxL)-(PLCf(pxL)-PLCf(x))*0.8, pxL, -1e-3, add=T)
       points(pxL, PLCf(pxL), pch=16, cex=1.2)})
text(-2, 0.16, expression(italic(p[k])==0.8), cex=1.3)

# PLC - wL
curve(PLCwLf, 0.1, 1,
      xlim=c(0, 1), ylim=c(0, 1),
      xlab=expression(italic(w)),
      ylab=expression(PLC[m]),
      cex.lab=1.3, col="red")
# pk=0.5
with(data.frame(wL=c(0.13555435977796931)),
     {curve(PLCwLf(wL)-(PLCwLf(wL)-PLCwLf(x))*0.5, wL, 1, add=T)
       points(wL, PLCwLf(wL), pch=16, cex=1.2)})
text(0.225, 0.628, expression(italic(w)==0.135), cex=1.3)
text(0.867, 0.35, expression(italic(p[k])==0.5), cex=1.3)
text(0.973, 0.96, "b", cex=1.3)
# pk=0.8
with(data.frame(wL=c(0.13555435977796931)),
     {curve(PLCwLf(wL)-(PLCwLf(wL)-PLCwLf(x))*0.8, wL, 1, add=T)
       points(wL, PLCwLf(wL), pch=16, cex=1.2)})
text(0.867, 0.16, expression(italic(p[k])==0.8), cex=1.3)
