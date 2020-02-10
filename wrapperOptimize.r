library(TMB)

#Create the data you need.
source("createData.r")

data <- list(p_c = as.matrix(t(p_c)), #observed chemical proportions
             p_b = as.matrix(t(p_b)), #observed biological proportions 
             c_i = nu, #observed carbon data
             temp = t) #temperature experiment

parameters <- list(beta_c = rep(0,nrow(data$p_c)-1), #chemical ratio parameters
                   A_par = matrix(0,nrow(data$p_c),nrow(data$p_b)-1), #biological parameters
                   alpha = rep(0,nrow(data$p_b)), #carbon intercept
                   gamma = rep(0,nrow(data$p_b)), #carbon slope
                   ln_sig = 0)# log-normal likelihood sd

try(dyn.unload("carbon"))
# compile("carbon.cpp")
dyn.load("carbon")
obj <- MakeADFun(data = data,
                 parameters = parameters,
                 DLL="carbon")

out <- nlminb(obj$par,obj$fn,obj$gr)
rep <- obj$report()
# SD <- sderport(obj)

par(mfrow=c(2,2))
plot(data$c_i,
     ylab="'Carbon' (units unkn)")
lines(exp(rep$chat_i))

matplot(t(data$p_c),
        pch=1,
        ylab="Chemical proportions")
matlines(t(rep$phat_c),
         lty=1,
         lwd=2)

matplot(t(data$p_b),
        pch=1,
        ylab="Bacteria proportions")
matlines(t(rep$phat_b),
         lty=1,
         lwd=2)

