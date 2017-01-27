
############################# Partial #############################
############################# Normal #############################
options(digits=20)
source("Partial/Normal plants/Functions.r")

ca <- 400

f1 <- function(wLr){
  averBif1 <- Vectorize(function(wLi)averBif(wLi, wLr, ca))
  res1 <- optimize(averBif1, c(0.13, wLr), tol=.Machine$double.eps, maximum=T)
  res <- c(wLr, res1$maximum, res1$objective)
  return(res)
}
f2 <- Vectorize(f1)

# Results
pkx <- 0.9
ca <- 400
k <- 0.025
MAP <- 3000

x <- seq(0.14, 0.22, by=0.001)
res1 <- f2(x)
res2 <- t(res1)
colnames(res2) <- c("wLr", "resLmax", "resLvalue")
write.csv(res2, "Partial/Normal plants/ESS lower wLi.csv", row.names=FALSE)
