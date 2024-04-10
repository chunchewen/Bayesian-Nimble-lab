library(mcmcplots)      
library(nimble)
library(spam)
library(lme4)
#------------------------------------------------------------------------------#
rm(list = ls())
#------------------------------------------------------------------------------#
set.seed(04142022)
n<-100                           # total subjects
nis<-rep(12,n)                   # repeated measure per subject (balanced)
id<-rep(1:n,nis)                 # id
N<-length(id)                    # total observations

# Fixed Effects
t<-rep(0,N)
for (i in 1:n) t[id==i]<-0:(nis[i]-1)# time variable: t=0,1,...,11
X<-cbind(1,t)                   # design matrix 
truebeta<-beta<-c(.25,-.5)      # true beta 
p<-ncol(X)                      # number of parameters


# Random Effects
truesigmab<-sigmab<-matrix(c(.50,.10,
                             .10,.15),2,2,byrow = T) # Allow cor. b/t random intercept/slope
trueB<-rmvnorm(n,sigma=truesigmab)

trueb1<-b1<-trueB[,1]           # true random intercept 
trueb2<-b2<-trueB[,2]           # true random slope     

B1<-rep(b1,nis)
B2<-rep(b2,nis)

# Precision Error
truesigmae<-sigmae<-4

# Outcome
mu<-X%*%beta+B1+B2*t
y<-rnorm(N,mu,sqrt(sigmae))              
                         
# Frequentist approach
m <- lmer(y~t+(1+t|id))
summary(m)

#--------#
# Nimble #
#--------#
code <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dnorm(mu[i],tau=taue)
    mu[i] <- beta0+beta1*t[i]+b1[id[i]]+b2[id[i]]*t[i]
  }
  # Random effect
  # update b=(b1i,b2i)^T ~ MVNORM(0,covb)
  for (j in 1:n) {
    B[j,1:2] ~ dmnorm(B0[1:2],cov=covb[1:2,1:2]) # N(0,sigmab) 
    b1[j]<-B[j,1]
    b2[j]<-B[j,2]}
  # Prior 
  covb[1:2,1:2] ~ dinvwish(C0[1:2,1:2],k)        # IW(diag(2),3)    
  beta0  ~ dnorm(0, sd = 100)
  beta1  ~ dnorm(0, sd = 100)
  taue ~ dgamma(0.001,0.001)
  vare <- 1/taue         
})
constants <- list(N=N, n = n,id=id, B0=rep(0,2), C0=diag(2),k=3)
data <- list(y = y, t = t)
inits <- list(beta0=0,beta1=0,covb=diag(2),B=matrix(0,n,2),taue=1)

# nimbleModel: check code specification is correct
model <- nimbleModel(code, constants = constants, data = data, inits = inits)


# Configure and Compile Model
modelConf <- configureMCMC(model, print = TRUE)
modelConf$addMonitors("vare") # add monitor error variance
Cmodel <- compileNimble(model)


# Compile MCMC Object
modelMCMC <- buildMCMC(modelConf)
modelMCMC <- compileNimble(modelMCMC, project = model)


# Run MCMC
niter <- 20000     # Number of iterations
burn  <- niter/2   # Number of burn-ins
thin  <- 20        # Thinning interval   
set.seed(1)
samples <- runMCMC(modelMCMC, niter = niter,nburnin = burn, nchains = 1,thin = thin,summary = TRUE)

# Extract samples
samples$summary  # estimate table


# Check MCMC plot
out<-samples$samples[,c("beta0","beta1",
                        "covb[1, 1]","covb[2, 2]","vare","taue")]
mcmcplot(out,filename = "output")



