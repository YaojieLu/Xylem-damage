
# Check http://www.sosmath.com/calculus/improper/testconv/testconv.html
options(digits=20)

# Scenario 2, family ESS
gswf <- Vectorize(function(w, wL=0.25,
                  a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1,
                  h=l*a*LAI/nZ*p, VPD=0.02, h2=l*LAI/nZ*p/1000,
                  pe=-1.58*10^-3, b=4.38, kxmax=5, c=9.53, d=1.28){
  pxmin <- pe*wL^(-b)
  kxmin <- kxmax*exp(-(-pxmin/d)^c)
  ps <- pe*w^(-b)
  res <- (ps-pxmin)*h2*kxmin/(h*VPD)
  return(res)
})
# Chapter 1
#gswf <- function(w)0.008680221*(w-0.2)^0.397299129
Ewf <- function(w, a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02)h*VPD*gswf(w)

f <- function(w)1/Ewf(w)
g <- function(w)1/(1e-5*(w-0.25)^0.9)

curve(g, 0.25, 0.25+1e-12)
curve(f, 0.25, 0.25+1e-12, add=T, col="red")

integrate(g, 0.25, 0.25+1e-5, rel.tol=.Machine$double.eps^0.25)
integrate(f, 0.25, 0.25+1e-5, rel.tol=.Machine$double.eps^0.25)
integrate(f, 0.25, 0.25+1e-5, rel.tol=.Machine$double.eps^0.01)
