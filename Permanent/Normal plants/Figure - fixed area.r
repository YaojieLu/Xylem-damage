
############################ Normal plants ############################
############################## Permanent ##############################
source("Permanent/Normal plants/Functions.r")
data <- read.csv("Permanent/Normal plants/LHS maxmum & RHS maximum (0.025, 3000).csv")
wLr <- data$wLr
wLiL <- with(data, ifelse(resLvalue>resHvalue, resLmax, NA))
wLiH <- with(data, ifelse(resLvalue>resHvalue, NA, resHmax))


# Figure
windows(8, 12)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(2, 1))
plot(wLr, wLiL,
     xlab=expression(italic(w[Lr])), ylab=expression(Optimal~italic(w[Li])),
     xlim=c(0, 1), ylim=c(0, 1), type="l", cex.lab=1.3)

points(wLr, wLiH, type="l")

abline(a=0, b=1, lwd=1)

plot(wLr, wLiL,
     xlab=expression(italic(w[Lr])), ylab=expression(Optimal~italic(w[Li])),
     xlim=c(0.14, 0.24), ylim=c(0.14, 0.24), type="l", cex.lab=1.3)

points(wLr, wLiH, type="l")

abline(a=0, b=1, lwd=1)
