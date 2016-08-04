
################################ Permanent ################################
############################ Vulnerable plants ############################
options(digits=20)
source("Permanent/Vulnerable plants/Functions.r")

ca <- 400
k <- 0.05
MAP <- 222

fLf1 <- Vectorize(fLf)
fLf1(0.227)
curve(fLf1, 0.21, 0.3)
