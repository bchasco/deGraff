library(mvtnorm)
nc <- 4 #number of chemicals observed 
ni <- 100 #number of experiments
proc_sig2_c <- 0.5 #variability in the samples
multi_logit <- function(coef){
  #coef: list of multinomial coefficients
  #predictor variables
  ml_i <- c(0,coef)
  z <- exp(ml_i)/(1+sum(exp(ml_i[2:length(ml_i)])))
  return(z)
}

beta_c <- c(-0.5,1,1.5)
Sigma <- matrix(0,nc-1,nc-1)
diag(Sigma) <- proc_sig2_c
c_i <- rmvnorm(ni,beta_c,Sigma)

p_c <- matrix(NA,ni,nc) 
for(i in 1:ni){
  p_c[i,] <- multi_logit(c_i[i,])
}


#Estimate the bacteria mixture from chemical mixture
nb <- 5 #number of bacteria
b_c <- rbind(c(0.25,0.5,1.,1.5)*.01, #mean bacteria effects
             matrix(runif((nb-1)*nc),nb-1,nc)) #chemical interaction effects)

p_b <- matrix(NA,ni,nb) #matrix of bacteria ratios

multi_logit2 <- function(coef,x){
  #coef: list of multinomial coefficients
  #predictor variables
  z_tmp <- c(1,rep(0,ncol(coef)))
  for(ii in 1:ncol(coef)){ #move across bacteria
    z_tmp[ii+1] <- exp(sum(coef[,ii]*x))
  }
  z <- z_tmp/sum(z_tmp)
  return(z)
}

p_x <- cbind(rep(1,ni),p_c) #get the chemical concentrations
#Use the multinomial to take the chemical concentrations and predict biological concentrations
for(i in 1:ni){
  p_b[i,] <- multi_logit2(b_c,p_x[i,])
}

#Now estimate carbon flux from bacteria mixture as a function of temperature
b_int <- runif(nb,min=1,max=10)
b_slope <- runif(nb,min=-1,max=1)
temps <- c(0,5,10,15)
t <- rep(temps,each=ni/length(temps))
nu <- rep(0,ni)
for(i in 1:ni){
  nu[i] <- sum(b_int + (b_slope * t[i])*p_b[i,])
}
plot(nu)
