
################################ Permanent ################################
############################ Vulnerable plants ############################
options(digits=20)
source("Permanent/Vulnerable plants/Functions.r")

# Results
ca <- 400
k <- 0.025
MAP <- 216

x <- seq(0.218, 0.219, by=0.0002)
resLmax <- vector(mode="numeric", length=length(x))
resHmax <- vector(mode="numeric", length=length(x))
resLvalue <- vector(mode="numeric", length=length(x))
resHvalue <- vector(mode="numeric", length=length(x))

for (i in 1:length(x)){
  resL <- optbelowf(x[i])
  resLmax[i] <- resL$maximum
  resLvalue[i] <- resL$objective
  resH <- optabovef(x[i])
  resHmax[i] <- resH$maximum
  resHvalue[i] <- resH$objective
  message(x[i])
}
res <- as.data.frame(cbind(x, resLmax, resHmax, resLvalue, resHvalue))
colnames(res) <- c("wLr", "resLmax", "resHmax", "resLvalue", "resHvalue")
#write.csv(res, "Permanent/Vulnerable plants/LHS maxmum & RHS maximum.csv", row.names = FALSE)

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(4, 4, 3, 2), mfrow=c(1,1))
plot(res$wLr, res$resLmax, xlim=c(head(res$wLr, 1), tail(res$wLr, 1)), lwd=2, ylim=c(head(res$wLr, 1), tail(res$wLr, 1)), type="l",
     xlab=expression(italic(w[Lr])), ylab=expression(italic(w[Li])),
     cex.lab=1.3, col="blue")
points(res$wLr, res$resHmax, type="l", cex.lab=1.3, col="red", lwd=2)
abline(a=0, b=1)
