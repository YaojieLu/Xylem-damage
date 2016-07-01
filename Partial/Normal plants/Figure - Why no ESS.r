
########################### Normal plants ###########################
############################## Partial ##############################
options(digits=3)
source("Partial/Normal plants/Functions.r")
datal <- read.csv("Partial/Normal plants/ESS lower wLi.csv")
datah <- read.csv("Partial/Normal plants/ESS higher wLi.csv")

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 4, 1, 1), mfrow=c(1, 1))
# LHS max vs. RHS max
plot(datal$wLr, datal$resLmax, xlim=c(0.14, 0.22), lwd=2, ylim=c(0.14, 0.22), type="l",
     xlab=expression(italic(w[Lr])), ylab=expression(Optimal~italic(w[Li])),
     cex.lab=1.3, col="blue")

points(datah$wLr, datah$resHmax, type="l", cex.lab=1.3, col="red", lwd=2)
abline(a=0, b=1, lwd=1)

legend("bottomright", title="Optimize between", expression("[0, "*italic(w[Lr])*"]", "["*italic(w[Lr])*", 1]"),
       lty=c(1, 1), lwd=c(2, 2), col=c("blue", "red"))
