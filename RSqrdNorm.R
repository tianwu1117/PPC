#The aim of this script is to calculate the R squared between the estimated and true polygenic score
#Weights to construct PRS are the posterior mean

RSqrdNorm <- function(h2,m,n,pi0) {
  #overall  
  VarBetaOverall <- h2/m
  VarBetaHatOverall <- VarBetaOverall + 1/n
  #beta=0(h0)
  VarBetaH0 <- 0
  VarBetaHatH0 <- 1/n
  SdVarBetaHatH0 <- sqrt(VarBetaHatH0)
  #beta!=0(h1)
  VarBetaH1 <- h2/(m*(1 - pi0))
  VarBetaHatH1 <- VarBetaH1 + 1/n
  SdVarBetaHatH1 <- sqrt(VarBetaHatH1)
  
  #simulate BetaHat
  BetaHat <- seq(-5*SdVarBetaHatH1 , 5*SdVarBetaHatH1 , length=1000)
  
  #f(BetaHat)
  height_h1 <- dnorm(BetaHat, mean = 0, sd = SdVarBetaHatH1) * (1 - pi0)
  height_h0 <- dnorm(BetaHat, mean = 0, sd = SdVarBetaHatH0) * pi0
  height_total <-  height_h0 + height_h1
  sum_height_total <- sum(height_total)
  
  #weight
  weight <- (height_h1/height_total) * (VarBetaH1/VarBetaHatH1) * BetaHat #conditional/posterior expectation

  #below to calculate R squared when weight is conditional expectation
  #E(weight*beta|betahat)
  E_weight_given_betahat <- weight * VarBetaH1 / VarBetaHatH1 * BetaHat * height_h1 / height_total
  #cov(weight,beta|betahat)=sum(f*E(weight*beta|betahat))/sum(f)
  cov_beta_weight <- height_total * E_weight_given_betahat
  var_weight <- weight * weight * height_total
    
  #denominators
  sum_cov_beta_weight <- sum(cov_beta_weight)
  sum_var_weight <- sum(var_weight)
    
  cov.weight.beta <- sum_cov_beta_weight / sum_height_total
  var.weight <- sum_var_weight / sum_height_total
  
  r_squared <- (cov.weight.beta / sqrt(var.weight * VarBetaOverall))**2
  r_squared_OLSE <- h2/(h2+m/n)
  return(c(r_squared_OLSE,r_squared))
}