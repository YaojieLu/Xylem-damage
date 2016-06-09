
options(digits=20)

# gswL(w)
gswLf <- function(w, wL,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                  pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05){
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  
  pxmin <- pe*wL^(-b)
  res <- (pe*w^-b-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}

# fL
fLf <- function(wL,
                a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=5.71, d=10.05,
                gamma=1/((MAP/365/k)/1000)*nZ){
  
  gswLf1 <- Vectorize(function(w)ifelse(w<wL+1e-10, 0, gswLf(w, wL)))
  
  Ef <- function(w){h*VPD*gswLf1(w)}
  rEf <- function(w){1/Ef(w)}
  integralrEf <- Vectorize(function(w){integrate(rEf, w, 1, rel.tol=.Machine$double.eps^0.25)$value})
  
  fnoc <- function(w)1/Ef(w)*exp(-gamma*(w-wL-1e-10)/(1-wL-1e-10)-k*integralrEf(w)*1/(1-wLr-1e-10))*1/(1-wLr-1e-10)
  res1 <- integrate(fnoc, wL+1e-10, 1, rel.tol=.Machine$double.eps^0.25)
  fLnoc <- 1/k*exp(-k*integrate(rEf, wL+1e-10, 1, rel.tol=.Machine$double.eps^0.25)$value*1/(1-wLr-1e-10))
  res <- fLnoc/(res1$value+fLnoc)
  return(res)
}

# Result
ca <- 400
k <- 0.05
MAP <- 1825
fLf1 <- Vectorize(fLf)
#curve(fLf1, 0.12, 0.16)
fLf1(0.129)
