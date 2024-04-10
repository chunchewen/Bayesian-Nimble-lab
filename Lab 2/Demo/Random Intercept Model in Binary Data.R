library(mcmcplots)      
library(nimble)
library(lme4)

#------------------------------------------------------------------------------#
rm(list = ls())
#------------------------------------------------------------------------------#
set.seed(04142022)
n<-500                           # total subjects
nis<-rep(12,n)                   # repeated measure per subject (balanced)
id<-rep(1:n,nis)                 # id
N<-length(id)                    # total observations

# Fixed Effects
t<-rep(0,N)
for (i in 1:n) t[id==i]<-0:(nis[i]-1)# time variable
X<-cbind(1,t)                    # design matrix 
p<-ncol(X)                       # number of parameters
truebeta<-beta<-c(.25,-.5)       # true beta 

# Random Effects
truesigmab<-sigmab<-0.5
trueb<-b<-rnorm(n,sd=sqrt(truesigmab))
B<-rep(b,nis)

# Binary outcome
eta<-X%*%beta+B
mu<-exp(eta)/(1+exp(eta))
y<-rbinom(N,1,mu)                  

# Frequentest approach
m <- glmer(y~t+(1|id),family = "binomial")
summary(m)

#--------#
# Nimble #
#--------#
code <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dbin(prob=p[i],size=1)
    logit(p[i]) <- beta0+beta1*t[i]+b[id[i]] # logit fn
  }
  # Random effect
  for(j in 1:n){b[j] ~ dnorm(0,tau=taub)} 
  # Prior
  taub  ~ dgamma(0.001,0.001)
  beta0  ~ dnorm(0, sd = 100)
  beta1  ~ dnorm(0, sd = 100)
  varb <- 1/taub 
})
constants <- list(N=N, n = n,id=id)
data <- list(y = y, t = t)
inits <- list(beta0=0,beta1=0,b=rep(0,n),taub=1)

# nimbleModel: check code specification is correct
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Configure the MCMC
modelConf <- configureMCMC(model, print = TRUE)
modelConf$addMonitors("varb") # add monitor random effect variance

# Compile the model
Cmodel <- compileNimble(model)

# Build MCMC object
modelMCMC <- buildMCMC(modelConf)
# Compile MCMC object
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
out<-samples$samples[,c("beta0","beta1","taub","varb")]
mcmcplot(out,filename = "output")





