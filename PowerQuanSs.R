#The aim of this script is to return the desirable sample size
##prob.given is the given probability to detect X sig. SNPs

PowerQuanSs <- function(prob.given,X,alpha,h2,m,pi0){
target <- function(ss) {
   sd.beta.hat.h1 <- sqrt(h2/((1-pi0)*m) + 1/ss)
   T.critical <- qnorm(alpha/2,0,sqrt(1/ss),lower.tail = F)
   exp.power  <- pnorm(T.critical,0,sd.beta.hat.h1,lower.tail = F) + pnorm(-T.critical,0,sd.beta.hat.h1)
   return(pbinom((X-1),round(m*(1-pi0)), exp.power) - (1-prob.given))
     }
# To find the ss that makes target near 0
ss.low <- 100
ss.high <- 10^9
ss.mid <- (ss.low+ss.high)/2

 while(abs(target(ss.mid))>10^-8){
   if (target(ss.mid) < 0){
     ss.high <- ss.mid
   }else{
     ss.low <- ss.mid
   }
   ss.mid <- (ss.low+ss.high)/2
 }
 return(max(1,round(ss.mid)))
}
