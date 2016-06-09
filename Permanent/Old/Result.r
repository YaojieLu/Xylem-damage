
source("Permanent/ESS.r")

#environmental conditions
ca <- c(400)  # Atmospheric CO2 concentration (ppm)
k <- c(0.025, 0.05, 0.1) # Rainfall frequency (per day)
MAP <- seq(1, 10, by=1)*365 # MAP=MDP*365; MAP: mean annual precipitation; MDP: mean daily precipitation
env <- as.vector(expand.grid(ca, k, MAP))

# Initialize
dvs <- matrix(nrow=nrow(env), ncol=1)

# Run every parameter combination
for(i in 1:nrow(env)){
  begin <- proc.time()
  
  ca <- env[i, 1]
  k <- env[i, 2]
  MAP <- env[i, 3]
  dvs[i, 1] <- optimize(f2, c(0.125, 0.132), tol=.Machine$double.eps)$minimum
  
  end <- proc.time()
  message(sprintf("%s/%s completed in %.2f min",i, nrow(env), (end[3]-begin[3])/60))
}

# Collect results
res <- cbind(env, dvs)
colnames(res) <- c("ca", "k", "MAP", "wL") 
write.csv(res, "Permanent/wL.csv", row.names = FALSE)
