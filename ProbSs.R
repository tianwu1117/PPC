#The aim of this script is to return the probability of detecting at least X significant SNPs
ProbExpNo <- function(alpha,X,h2,m,n,pi0){
  if(h2/((1-pi0)*m) + 1/n < 0 | X > m*(1-pi0)){
   stop("Error")
  }else{
  sd.beta.hat.h1 <- sqrt(h2/((1-pi0)*m) + 1/n)
  T.critical <- qnorm(alpha/2,0,sqrt(1/n),lower.tail = F)
#Probability of finding at least 1 significant true SNP
  Power <- pnorm(T.critical,0,sd.beta.hat.h1,lower.tail = F) + pnorm(-T.critical,0,sd.beta.hat.h1)
  Type2 <- 1- Power
  P.X <- 1-pbinom((X-1),round(m * (1-pi0)), Power)

  return(paste0(round(P.X*100,2),"%"))
  }
}
