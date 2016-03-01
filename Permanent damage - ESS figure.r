
gsf1 <- function(w,
                 a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                 pe=-1.58*10^-3, b=4.38, h2=h/1000, kmax=5, c=5.71, d=10.05){
  h2/(h*VPD)*pe*kmax*(w^(-b)-w0^(-b))*exp(-(-pe*w0^(-b)/d)^c)}
gsf2 <- Vectorize(function(w)ifelse(w>w0, gsf1(w), 0))

# Figures
windows(8, 6)
par(mgp = c(2.2, 1, 0), xaxs = "i", yaxs = "i", lwd = 2, mar=c(4, 4, 2.5, 1), mfrow=c(1,1))
w0 <- 0.61887864661125191
curve(gsf2,0, 1, ylim=c(0, 0.25))
w0 <- 0.63197148074076115
curve(gsf2,0, 1, add=T)
w0 <- 0.62412701543078464
curve(gsf2,0, 1, add=T)
w0 <- 0.61820904273854638
curve(gsf2,0, 1, add=T)
w0 <- 0.63025774898667186
curve(gsf2,0, 1, add=T)
