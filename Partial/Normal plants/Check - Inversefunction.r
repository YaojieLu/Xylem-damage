
############################# Partial #############################
############################# Normal #############################
source("Partial/Normal plants/Functions.r")

ca <- 400
pkx <- 0.95
wL <- 0.18

wLL <- wLLf(wL, ca)

f1 <- Vectorize(function(w)gsmaxfm(w, wL))
f2 <- Vectorize(function(w)gswLf(w, wL, ca))
inversegsmaxfm <- function(gs)optimize(function(w)(gsmaxfm(w, wL)-gs)^2, c(wLL$root, 1), tol=.Machine$double.eps)$minimum
inversegswLf <- function(gs)optimize(function(w)(gswLf(w, wL, ca)-gs)^2, c(wLL$root, 1), tol=.Machine$double.eps)$minimum
f3 <- Vectorize(inversegsmaxfm)
f4 <- Vectorize(inversegswLf)

x3 <- seq(f1(wLL$root), f1(1), by=(f1(1)-f1(wLL$root))/100)
y3 <- f3(x3)
x4 <- seq(f2(wLL$root), f2(1), by=(f2(1)-f2(wLL$root))/100)
y4 <- f4(x4)

curve(f1, wLL$root, 1)
points(y3, x3)
curve(f2, wLL$root, 1)
points(y4, x4)
