#The aim of this function is to calculate the phenotypic variance
#under the extreme sample selection design

#TL is the lower threshold of the phenotype, TU is the upper threshold of the phenotype
#P.neg is the proportion of the extremely small samples in total sample
samvar <- function(TL,TU,P.neg) {

P.pos <- 1-P.neg
  
exp.sam.neg <- -dnorm(TL)/pnorm(TL)
exp.sam.pos <- dnorm(TU)/(1-pnorm(TU))

var.sam.neg <- 1 - TL * dnorm(TL)/pnorm(TL)-(dnorm(TL)/pnorm(TL))^2
var.sam.pos <- 1 + TU * dnorm(TU)/(1-pnorm(TU)) - (dnorm(TU)/(1-pnorm(TU)))^2

var.extreme <- var.sam.neg * P.neg + var.sam.pos * P.pos + 
              exp.sam.neg^2 * (1-P.neg)* P.neg + exp.sam.pos^2 * (1-P.pos) *P.pos -
              2 * exp.sam.neg * exp.sam.pos * P.neg * P.pos
return(var.extreme)
}