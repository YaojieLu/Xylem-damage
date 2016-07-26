
############################ Normal plants ############################
############################## Permanent ##############################
options(digits=3)
source("Permanent/Normal plants/Functions.r")
data <- read.csv("Permanent/Normal plants/LHS maxmum & RHS maximum (0.025, 3000).csv")

dataL <- subset(data, wLr<=0.3 & resLvalue>resHvalue, select=c(wLr, resLmax))
dataH <- subset(data, wLr<=0.3 & resLvalue<resHvalue, select=c(wLr, resHmax))

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(4, 4, 2, 2), mfrow=c(1,1))
plot(dataL$wLr, dataL$resLmax, xlim=c(0.14, 0.3), ylim=c(0.16, 0.19), type="l",
     xlab=expression(italic(w[Lr])), ylab=expression(italic(w[Li])),
     cex.lab=1.3, col="blue", lwd=2)
points(dataH$wLr, dataH$resHmax, type="l", cex.lab=1.3, col="red", lwd=2)
abline(a=0, b=1)