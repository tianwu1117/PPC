library(data.table)
#The aim of this script is to calculate the 1)expected variance explained by the significant / true significant SNPs
#2) calculate the variance of variance explained by the significant / true significant SNPs

#to calculate the expectation and variance of Power
Exp.No.Exp.Var <- function(h2,m,n,pi0,alpha){
  var_beta_h1 <- h2/((1-pi0)*m)
  sd.beta.h1 <- sqrt(var_beta_h1)

  var_beta_hat_h1 <- h2/((1-pi0)*m) + 1/n
  sd.beta.hat.h1 <- sqrt(var_beta_hat_h1)

  T.critical <- qnorm(alpha/2,0,sqrt(1/n),lower.tail = F)
  #split beta.hat into beta.hat.i and put them in a data table
  beta_hat <- seq(-15*sd.beta.hat.h1, 15*sd.beta.hat.h1, length=10000)

  Area <- pi0 * pnorm(T.critical, 0, sqrt(1/n),lower.tail = F) + (1-pi0) * pnorm(T.critical,0,sd.beta.hat.h1,lower.tail = F)
  Exp.no <- round(m * 2 * Area, 2)

  Area.true <- (1-pi0) * pnorm(T.critical,0,sd.beta.hat.h1,lower.tail = F)
  Exp.True.no <- round(m * 2 * Area.true, 2)
###########################################################################################################################
  #expectation and variance of power
###########################################################################################################################
  beta <- seq(-5*sd.beta.h1, 5*sd.beta.h1, length=5000)
  T.critical.chisqrd <- qchisq(1-alpha,1)

  ##to construct data table for BETA
  DT <- data.table(beta=beta)
  DT[, dens.beta := (1-pi0)*dnorm(beta,0,sd.beta.h1)]

  #NCP
  DT[, NCP:= n*beta^2]
  DT[, power := pchisq(T.critical.chisqrd,1,ncp=NCP, lower.tail = FALSE)]
  DT[, power.sqrd := power^2]

  #total frequency - denominator of expectation of power and power sqrd
  total.freq <- sum(DT$dens.beta)
  #nominator of expectation of power
  DT[, fi.powi := dens.beta * power]
  #Expectation of power
  exp.power <- sum(DT$fi.powi)/total.freq

  #nominator of the expectation of power squared
  DT[, fi.powi.sqrd := dens.beta * power.sqrd]

  #Expectation of power squared
  exp.pow.sqrd <- sum(DT$fi.powi.sqrd)/total.freq

  #variance of power
  var.power <- exp.pow.sqrd-exp.power^2

###########################################################################################################################
  #The expected variance explained by significant SNPs
###########################################################################################################################
  DT <- data.table(beta_hat=beta_hat)
  #the *ratio* of height is the probability that for given beta.hat, that beta comes from certain distribution (null or alternative)
  DT[, height_h1 := dnorm(beta_hat, mean = 0, sd = sd.beta.hat.h1) * (1 - pi0)]
  DT[, height_h0 := dnorm(beta_hat, mean = 0, sd = sqrt(1/n)) * pi0]
  DT[, height_total:= height_h0 + height_h1]

  #judge whether beta.hat is significant,0/1
  DT[, Sig_beta_hat := ifelse(abs(beta_hat)>T.critical,1,0)]

  #frequency of beta.hat. If beta.hat is not significant, frequency will be 0.
  DT[, Sig_height_total := height_total * Sig_beta_hat]

  #different replacement of beta
    ##the apparent one
    DT[, beta_hat_sqrd := beta_hat^2]
  #The apparent variance explained: using beta.hat^2 to replace beta^2 when beta.hat is significant
    DT[, Exp.sig.beta_hat_sqrd :=  beta_hat_sqrd *  Sig_height_total]

  #To calculate the variance of variance explained - beta.hat^4, single term of the nominator
    DT[, Exp.sig.beta_hat_poweroffour := beta_hat^4 * Sig_height_total]

    ##the conditional expectation = K*beta.hat, where K = var(beta)/var(beta.hat)*height ratio
    DT[, CE:= var_beta_h1 / var_beta_hat_h1 * height_h1 / height_total * beta_hat]

  ##To calculate the variance explained by significant SNPs, E(beta.hat^2) - (E(beta.hat))^2.
  ##To replace beta.hat^2 by E(beta^2|beta.hat,H1j) = (var(beta|beta.hat,H1j)+(E(beta|beta.hat,H1j))^2)*P(H1j|zj)
    DT[, CE2 := ( var_beta_h1 / var_beta_hat_h1 / n + (var_beta_h1 / var_beta_hat_h1 * beta_hat)^2) * (height_h1/height_total) ]
  # ##TRUE1: using E(beta|beta.hat)^2 to replace beta^2, single term of the nominator
  # DT[, Exp.CE.sqrd := CE^2 * Sig_height_total]

    ##TRUE2: using E(beta^2|beta.hat) to directly replace beta^2,single term of the nominator
    DT[, Exp.CE2 := CE2 * Sig_height_total]

    ##E((beta.sqrd|beta.hat)^2), single term of the nominator
    DT[, Exp.CE2.sqrd := CE2^2* Sig_height_total]

    #The expected variance explained of a single SNP: to calculate the expected variance explained by significant SNPs
    single.Exp.var.exp <- sum(DT[,Exp.sig.beta_hat_sqrd]) / sum(DT[,Sig_height_total])

    # single.true.var.exp <- sum(DT[,Exp.CE.sqrd]) / sum(DT[,Sig_height_total])
    #expected TRUE variance explained
    single.true2.var.exp <- sum(DT[,Exp.CE2]) / sum(DT[,Sig_height_total])

    #to calculate the variance of variance explained by significant SNPs
    single.var.var.exp.first <- sum(DT[,Exp.sig.beta_hat_poweroffour])/sum(DT[,Sig_height_total])
    single.var.var.exp.second <- sum(DT[,Exp.sig.beta_hat_sqrd])/sum(DT[,Sig_height_total])
    single.var.var.exp <- single.var.var.exp.first - single.var.var.exp.second^2

    #to calculate the variance of variance explained by TRUE significant SNPs
    single.var.true.var.exp.first <- sum(DT[,Exp.CE2.sqrd])/sum(DT[,Sig_height_total])
    single.var.true.var.exp <- single.var.true.var.exp.first - single.true2.var.exp^2

    ###########################################################################################################################
    #results
    ###########################################################################################################################

    #variance of no. of expected number of significant SNPs
    var.no <- m*pi0*alpha*(1-alpha)+m*(1-pi0)*(exp.power*(1-exp.power)-var.power)
    SE.no.sig.SNP <- sqrt(var.no)

    var.true.no <- m*(1-pi0)*(exp.power*(1-exp.power)-var.power)
    SE.no.true.sig.SNP <- sqrt(var.true.no)

    #variance of variance explained by significant SNPs, Var(Y) = E(Var(Y|X)) + Var(E(Y|X)), X is the number of significant SNPs
    #A reference:https://ocw.mit.edu/resources/res-6-012-introduction-to-probability-spring-2018/part-i-the-fundamentals/variance-of-the-sum-of-a-random-number-of-random-variables/

    var.var.exp2 <- Exp.no*single.var.var.exp + var.no*single.Exp.var.exp^2
    var.true2.var.exp <- Exp.True.no*single.var.true.var.exp + var.true.no*single.true2.var.exp^2

    #SD of variance explained by significant SNPs
    SE.var.exp <- sqrt(var.var.exp2)
    #SD of variance explained by TRUE significant SNPs
    SE.true2.var.exp <- sqrt(var.true2.var.exp)

    #Res1: expected variance explained by significant SNPs
    res.exp.var.exp <- single.Exp.var.exp*Exp.no
    #Res2: expected variance explained by TRUE significant SNPs
    res.exp.var.exp.true <- single.true2.var.exp*Exp.True.no

    return(
      c(
      formatC(Exp.no, format="f", big.mark=",", digits=2),        #[1] The expected no. of significant SNPs
      formatC(SE.no.sig.SNP, format="f", big.mark=",", digits=2),        #[2] The SE of the number of significant SNPs
      paste(round(res.exp.var.exp*100,2),"%",sep = ""),      #[3] The expected variance explained by the significant SNPs
      paste(round(SE.var.exp*100,2),"%",sep = ""),           #[4] The SE of the variance explained by the significant SNPs
      formatC(Exp.True.no, format="f", big.mark=",", digits=2),          #[5] The expected no. of TRUE significant SNPs
      formatC(SE.no.true.sig.SNP,format="f", big.mark=",", digits=2),   #[6] The SE of the no. of TRUE significant SNPs
      paste(round(res.exp.var.exp.true*100,2),"%",sep = ""), #[7] The expected variance explained by the TRUE significant SNPs
      paste(round(SE.true2.var.exp*100,2),"%",sep = "")      #[8] The SE of the variance explained by the TRUE significant SNPs
        )
           )

}
