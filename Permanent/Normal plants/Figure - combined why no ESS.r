
############################ Normal plants ############################
############################## Permanent ##############################
options(digits=3)
source("Permanent/Normal plants/Functions.r")
data <- read.csv("Permanent/Normal plants/LHS maxmum & RHS maximum (0.025, 3000).csv")

pxminf1 <- Vectorize(pxminf)

ca <- 400
k <- 0.025
MAP <- 3000

# Figure
wLr <- c(0.16, 0.171, 0.18)
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(2, 2))
# LHS max vs. RHS max
plot(data$wLr, data$resLmax, xlim=c(0.14, 0.24), lwd=2, ylim=c(0.14, 0.24), type="l",
     xlab=expression(italic(w[Lr])), ylab=expression(Optimal~italic(w[Li])),
     cex.lab=1.3, col="blue")
points(data$wLr, data$resHmax, type="l", cex.lab=1.3, col="red", lwd=2)
abline(a=0, b=1, lwd=1)
with(data.frame(wLr=wLr[1]), {
  segments(wLr, 0.14, wLr, wLr, lty=2, lwd=1, col="blue")
  segments(wLr, 0.24, wLr, wLr, lty=2, lwd=1, col="red")
  text(wLr+0.004, 0.21, "b", cex=1.5)
})
with(data.frame(wLr=wLr[2]), {
  segments(wLr, 0.14, wLr, wLr, lty=2, lwd=1, col="blue")
  segments(wLr, 0.24, wLr, wLr, lty=2, lwd=1, col="red")
  text(wLr+0.004, 0.21, "c", cex=1.5)
})
with(data.frame(wLr=wLr[3]), {
  segments(wLr, 0.14, wLr, wLr, lty=2, lwd=1, col="blue")
  segments(wLr, 0.24, wLr, wLr, lty=2, lwd=1, col="red")
  text(wLr+0.004, 0.21, "d", cex=1.5)
})
text(0.14+0.1*0.05, 0.24-0.1*0.06, "a", cex=1.5)
legend("bottomright", title="Optimum within", expression("[0, "*italic(w[Lr])*"]", "["*italic(w[Lr])*", 1]"),
       lty=c(1, 1), lwd=c(2, 2), col=c("blue", "red"))
box()

# Given wLr, varing wLi
with(data.frame(wLr=wLr[1]),
     {averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
     curve(averBif1, 0.14, wLr,
           xlim=c(0.14, 0.24), ylim=c(0, 15), yaxt="n",
           xlab=expression(italic(w[Li])), ylab=expression(bar(italic(A[Ni]))~(mu*mol~m^-2~s^-1)),
           lty=2, cex.lab=1.3, col="blue")
     curve(averBif1, wLr, 0.24,
           lty=2, add=T, col="red")
     segments(data[data$wLr==wLr, "resLmax"], 0, data[data$wLr==wLr, "resLmax"], data[data$wLr==wLr, "resLvalue"], lwd=1, col="blue")
     segments(data[data$wLr==wLr, "resHmax"], 0, data[data$wLr==wLr, "resHmax"], data[data$wLr==wLr, "resHvalue"], lwd=1, col="red")
     points(wLr, averBif1(wLr), pch=4, cex=1.5)
     axis(2, at=c(0, 5, 10, 15), pos=0.14, lwd=2)
     text(0.14+0.1*0.05, 15*0.94, "b", cex=1.5)
     box()
     })

with(data.frame(wLr=wLr[2]),
     {averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
     curve(averBif1, 0.14, wLr,
           xlim=c(0.14, 0.24), ylim=c(0, 15), yaxt="n",
           xlab=expression(italic(w[Li])), ylab=expression(bar(italic(A[Ni]))~(mu*mol~m^-2~s^-1)),
           lty=2, cex.lab=1.3, col="blue")
     curve(averBif1, wLr, 0.24,
           lty=2, add=T, col="red")
     segments(data[data$wLr==wLr, "resLmax"], 0, data[data$wLr==wLr, "resLmax"], data[data$wLr==wLr, "resLvalue"], lwd=1, col="blue")
     segments(data[data$wLr==wLr, "resHmax"], 0, data[data$wLr==wLr, "resHmax"], data[data$wLr==wLr, "resHvalue"], lwd=1, col="red")
     points(wLr, averBif1(wLr), pch=4, cex=1.5)
     axis(2, at=c(0, 5, 10, 15), pos=0.14, lwd=2)
     text(0.14+0.1*0.05, 15*0.94, "c", cex=1.5)
     box()
     })

with(data.frame(wLr=wLr[3]),
     {averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
     curve(averBif1, 0.14, wLr,
           xlim=c(0.14, 0.24), ylim=c(0, 15), yaxt="n",
           xlab=expression(italic(w[Li])), ylab=expression(bar(italic(A[Ni]))~(mu*mol~m^-2~s^-1)),
           lty=2, cex.lab=1.3, col="blue")
     curve(averBif1, wLr, 0.24,
           lty=2, add=T, col="red")
     segments(data[data$wLr==wLr, "resLmax"], 0, data[data$wLr==wLr, "resLmax"], data[data$wLr==wLr, "resLvalue"], lwd=1, col="blue")
     segments(data[data$wLr==wLr, "resHmax"], 0, data[data$wLr==wLr, "resHmax"], data[data$wLr==wLr, "resHvalue"], lwd=1, col="red")
     points(wLr, averBif1(wLr), pch=4, cex=1.5)
     axis(2, at=c(0, 5, 10, 15), pos=0.14, lwd=2)
     text(0.14+0.1*0.05, 15*0.94, "d", cex=1.5)
     box()
     })

dev.copy2pdf(file = "Figures/Figure 8.pdf")
