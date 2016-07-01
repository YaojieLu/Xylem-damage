
options(digits=3)
#######################################################################################################
# Scenario 1
#######################################################################################################
# gsmax
gsmaxf1 <- function(w,
                    a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                    pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)-((ps-x)*h2*kxf(x))
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  ps <- pe*w^(-b)
  pxmin <- pxminf(w)
  res <- (ps-pxmin)*h2*kxf(pxmin)/(h*VPD)
  return(res)
}
# Af(gs, ca)
Af1 <- function(gs, ca,
                Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))
# mf(w, gs)
mf1 <- function(w, gs,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54, h3=10){
  
  # xylem conductance function
  kxf <- function(x)kxmax*exp(-(-x/d)^c)
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)(ps-x)*h2*kxf(x)
    res <- optimize(f1, c(-20,0), tol=.Machine$double.eps, maximum=T)$maximum
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxf(x)-h*VPD*gs)^2
    res <- optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum
    return(res)
  }
  
  px <- pxf(w, gs)
  kx <- kxf(px)
  PLC <- 1-kx/kxmax
  res <- h3*PLC
  return(res)
}

# B(w, gs)
Bf1 <- function(w, gs, ca)Af1(gs, ca)-mf1(w, gs)

# ESS stomatal behaviour function
ESSf1 <- function(w, ca){
  f1 <- function(gs)Bf1(w, gs, ca)
  res <- optimize(f1, c(0, gsmaxf1(w)), tol=.Machine$double.eps, maximum=T)
  return(res$maximum)
}
#######################################################################################################
# Scenario 3
#######################################################################################################
# the original PLC(px)
PLCf3 <- function(px, c=2.64, d=3.54)1-exp(-(-px/d)^c)

# modified gsmax
gsmaxfm3 <- function(w, wL,
                     a=1.6, nZ=0.5, p=43200, l=1.8e-5, LAI=1, h=l*a*LAI/nZ*p, VPD=0.02,
                     pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54){
  
  pxL <- pe*wL^(-b)
  # modified PLC
  PLCfm <- function(x)PLCf3(pxL)-(PLCf3(pxL)-PLCf3(x))*pkx
  # modified xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  
  ps <- pe*w^(-b)
  f1 <- function(x)(ps-x)*h2*kxfm(x)/(h*VPD)
  
  res <- ifelse(pxL<ps, optimize(f1, c(pxL,ps), tol=.Machine$double.eps, maximum=T)$objective, 0)
  return(res)
}

# Af(gs)
Af3 <- function(gs, Vcmax=50, cp=30, Km=703, Rd=1, LAI=1)LAI*1/2*(Vcmax+(Km+ca)*gs-Rd-((Vcmax)^2+2*Vcmax*(Km-ca+2*cp)*gs+((ca+Km)*gs+Rd)^2-2*Rd*Vcmax)^(1/2))

# modified mf(w, gs)
mf3 <- function(w, gs, wL,
                a=1.6, LAI=1, nZ=0.5, p=43200, l=1.8e-5, h=l*a*LAI/nZ*p, VPD=0.02,
                pe=-1.58*10^-3, b=4.38, h2=h/1000, kxmax=5, c=2.64, d=3.54, h3=10){
  
  pxL <- pe*wL^(-b)
  # modified PLC
  PLCfm <- function(x)PLCf3(pxL)-(PLCf3(pxL)-PLCf3(x))*pkx
  # xylem conductance function
  kxfm <- function(x)kxmax*(1-PLCfm(x))
  # minimum xylem water potential function for given w
  pxminf <- function(w){
    ps <- pe*w^(-b)
    f1 <- function(x)(ps-x)*h2*kxfm(x)
    res <- ifelse(pxL<ps, optimize(f1, c(pxL, ps), tol=.Machine$double.eps, maximum=T)$maximum, ps)
    return(res)
  }
  # xylem water potential function
  pxf <- function(w, gs){
    ps <- pe*w^(-b)
    pxmin <- pxminf(w)
    f1 <- function(x)((ps-x)*h2*kxfm(x)-h*VPD*gs)^2
    res <- ifelse(pxmin<ps, optimize(f1, c(pxmin, ps), tol=.Machine$double.eps)$minimum, ps)
    return(res)
  }
  
  px <- pxf(w, gs)
  kx <- kxfm(px)
  PLC <- 1-kx/kxmax
  res <- h3*(PLC-PLCfm(0))
  return(res)
}

# modified B(w, gs)
Bf3 <- function(w, gs, wL)Af3(gs)-mf3(w, gs, wL)

# family ESS
gswLf3 <- function(w, wL){
  Bfm1 <- function(gs)Bf3(w, gs, wL)
  gsmaxfm1 <- function(w)gsmaxfm3(w, wL)
  res <- ifelse(0<gsmaxfm1(w), optimize(Bfm1, c(0, gsmaxfm1(w)), tol=.Machine$double.eps, maximum=T)$maximum, 0)
  return(res)
}

# LHS endpoint
wLLf3 <- function(wL){
  Bfm1 <- function(w)Bf3(w, gswLf3(w, wL), wL)
  res <- uniroot(Bfm1, c(wL, 1), tol=.Machine$double.eps)
  return(res)
}
#######################################################################################################

#######################################################################################################
## wL
## Initialize
#ca <- 400
#pkx <- 0.5
#wL <- c(0.2, 0.3)
#Cols <- c("blue", "orange")
#
#ESSf11 <- Vectorize(function(w)ESSf1(w, ca))
#
## Figure
#windows(8, 6)
#par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1, 1))
#curve(ESSf11, 0.17251375297104052, 1,
#      xlim=c(0, 1), ylim=c(0, 0.5),
#      xlab=expression(italic(w)),
#      ylab=expression(italic(g[s])~(mol~m^-2~s^-1)),
#      cex.lab=1.3)
#
#with(list(wL=wL[1], Col=Cols[1]), {
#  gswLf31 <- Vectorize(function(w)gswLf3(w, wL))
#  gsmaxfm31 <- Vectorize(function(w)gsmaxfm3(w, wL))
#  wLL <- wLLf3(wL)
#  curve(gswLf31, wLL$root, 1, col=Col, add=T)
#  curve(gsmaxfm31, wL, 1, col=Col, add=T, type="p")
#})
#with(list(wL=wL[2], Col=Cols[2]), {
#  gswLf31 <- Vectorize(function(w)gswLf3(w, wL))
#  gsmaxfm31 <- Vectorize(function(w)gsmaxfm3(w, wL))
#  wLL <- wLLf3(wL)
#  curve(gswLf31, wLL$root, 1, col=Col, add=T)
#  curve(gsmaxfm31, wL, 1, col=Col, add=T, type="p")
#})
#
#legend("topleft", c("wL=0.2, gsmax, S3", "wL=0.3, gsmax, S3", "wL=0.2, ESS, S3", "wL=0.3, ESS, S3", "ESS, S1"), lty=c(NA, NA, 1, 1, 1), pch=c(1, 1, NA, NA, NA), col=c(Cols, Cols, "black"))
#
#######################################################################################################
# pkx
# Initialize
ca <- 400
pkx1 <- c(0.5, 0.99)
wL <- 0.2
Cols <- c("blue", "orange")

ESSf11 <- Vectorize(function(w)ESSf1(w, ca))

# Figure
windows(8, 6)
par(mgp=c(2.2, 1, 0), xaxs="i", yaxs="i", lwd=2, mar=c(4, 4, 1, 1), mfrow=c(1, 1))
curve(ESSf11, 0.17251375297104052, 1,
      xlim=c(0, 1), ylim=c(0, 0.5),
      xlab=expression(italic(w)),
      ylab=expression(italic(g[s])~(mol~m^-2~s^-1)),
      cex.lab=1.3)

with(list(pkx <- c(pkx1[1]), Col=Cols[1]), {
  gswLf31 <- Vectorize(function(w)gswLf3(w, wL))
  gsmaxfm31 <- Vectorize(function(w)gsmaxfm3(w, wL))
  wLL <- wLLf3(wL)
  curve(gswLf31, wLL$root, 1, col=Col, add=T)
  curve(gsmaxfm31, wL, 1, col=Col, add=T, type="p")
})
with(list(pkx <- c(pkx1[2]), Col=Cols[2]), {
  gswLf31 <- Vectorize(function(w)gswLf3(w, wL))
  gsmaxfm31 <- Vectorize(function(w)gsmaxfm3(w, wL))
  wLL <- wLLf3(wL)
  curve(gswLf31, wLL$root, 1, col=Col, add=T)
  curve(gsmaxfm31, wL, 1, col=Col, add=T, type="p")
})

legend("topleft", c("pkx=0.5, gsmax, S3", "pkx=0.99, gsmax, S3", "pkx=0.5, ESS, S3", "pkx=0.99, ESS, S3", "ESS, S1"), lty=c(NA, NA, 1, 1, 1), pch=c(1, 1, NA, NA, NA), col=c(Cols, Cols, "black"))
