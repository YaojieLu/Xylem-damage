
################################ Permanent ################################
############################ Vulnerable plants ############################
options(digits=20)
source("Permanent/Vulnerable plants/Functions.r")

# Results
ca <- 400
k <- 0.025
MAP <- 365

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", mar=c(4, 4, 2, 2), mfrow=c(1,1))
with(data.frame(wLr=c(0.2196)),
     {f <- Vectorize(function(wLi)averBif(wLi, wLr))
     curve(f, 0.2191, 0.2201)
     abline(v=wLr)})
