
############################ Normal plants ############################
############################## Permanent ##############################
options(digits=20)
source("Permanent/Normal plants/Functions.r")

# Auxiliary function
# Optimization among wLi < wLr
f1 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(0.001, wLr), tol=.Machine$double.eps, maximum=T)
  return(res)
}
# Optimization among wLi > wLr
f2 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr))
  res <- optimize(averBif1, c(wLr, 0.999), tol=.Machine$double.eps, maximum=T)
  return(res)
}

# Results
ca <- 400
k <- 0.025
MAP <- 3000

x <- seq(0.139, 0.99, by=0.001)
resLmax <- vector(mode="numeric", length=length(x))
resHmax <- vector(mode="numeric", length=length(x))
resLvalue <- vector(mode="numeric", length=length(x))
resHvalue <- vector(mode="numeric", length=length(x))

for (i in 1:length(x)){
  resL <- f1(x[i])
  resLmax[i] <- resL$maximum
  resLvalue[i] <- resL$objective
  resH <- f2(x[i])
  resHmax[i] <- resH$maximum
  resHvalue[i] <- resH$objective
  message(x[i])
}
res <- cbind(x, resLmax, resHmax, resLvalue, resHvalue)
colnames(res) <- c("wLr", "resLmax", "resHmax", "resLvalue", "resHvalue")
write.csv(res, "Permanent/Normal plants/LHS maxmum & RHS maximum (0.025, 3000).csv", row.names = FALSE)
