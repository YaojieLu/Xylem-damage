
source("Partial/gswL(w).r")

ca <- 400
k <- 0.025
MAP <- 365
pkx <- 0.5

# fL
fLf <- function(wL,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05,
                gamma=1/((MAP/365/k)/1000)*nZ){
  
  wLL <- wLLf(wL)$root
  gswLf1 <- Vectorize(function(w)gswLf(w, wL))
  Ef <- function(w){h*VPD*gswLf1(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.4)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wLL)/(1-wLL)-k*integralrEf(w)*1/(1-wLL))*1/(1-wLL)
  res1 <- integrate(fnoc, wLL, 1, rel.tol=.Machine$double.eps^0.4)
  res2 <- 1/k*exp(-k*integrate(rEf, wLL, 1, rel.tol=.Machine$double.eps^0.4)$value*1/(1-wLL))
  res <- res2/(res2+res1$value)
  message("fL(", wL, ")=", res)
  return(res)
}

# Result
fLf1 <- Vectorize(fLf)
wL <- seq(0.12, 0.145, by=0.005)
fL <- fLf1(wL)
res <- cbind(wL, fL)
write.csv(res, "Partial/fL.csv", row.names = FALSE)

# Figure
data <- read.csv("Partial/fL.csv")
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(3.5, 3.5, 1, 1.5), mfrow=c(1,1))
plot(data$wL, data$fL, xlim=c(0.12, 0.145), ylim=c(0.8, 0.9), type='l', cex.lab=1.3,
     xlab=expression(w[L]), ylab=expression(f[L]))
