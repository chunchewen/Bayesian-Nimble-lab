library(mcmcplots)      
library(nimble)
library(lme4)
#------------------------------------------------------------------------------#
rm(list = ls())
#------------------------------------------------------------------------------#
# Generate Data
set.seed(04142022)
n<-100                           # total subjects
nis<-rep(12,n)                   # repeated measure per subject (balanced)
id<-rep(1:n,nis)                 # id = (1,1,..,1,2,2,...,2,...,n,n,...,n)
N<-length(id)                    # total observations


# Fixed Effects
t<-rep(0,N)
for (i in 1:n) t[id==i]<-0:(nis[i]-1)# time variable
X<-cbind(1,t)                        # design matrix 
truebeta<-beta<-c(.25,-.5)           # true beta 
p<-ncol(X)                           # number of parameters

# Random effect variance
truesigmab<-sigmab<-0.5
trueB<-b<-rnorm(n,sd=sqrt(truesigmab))
B<-rep(b,nis)


# Error variance
truesigmae<-sigmae<-4

# Outcome Model
mu<-X%*%beta+B
y<-rnorm(N,mu,sqrt(sigmae))  

# Frequentist Approach
m <- lmer(y~t+(1|id))
summary(m)


#--------#
# Nimble #
#--------#
code <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dnorm(mu[i],tau=taue)
    mu[i] <- beta0+beta1*t[i]+b[id[i]]
  }
  # Random Effect
  for(j in 1:n){b[j] ~ dnorm(0,tau=taub)}  # Random effects
  # Prior 
  taub  ~ dgamma(0.0001,0.0001)            # Random effect precision 
  taue  ~ dgamma(0.0001,0.0001)            # Error precision
  beta0  ~ dnorm(0, sd = 100)              
  beta1  ~ dnorm(0, sd = 100)
  varb <- 1/taub                           # Random effect variance 
  vare <- 1/taue                           # Error variance
})

constants <- list(N=N, n = n,id=id)        # Note: ID variable
data <- list(y = y, t = t)
inits <- list(beta0=0,beta1=0,b=rep(0,n),taub=1,taue=1)

# nimbleModel: check code specification is correct
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Configure and Compile Model
modelConf <- configureMCMC(model, print = TRUE)
modelConf$addMonitors("varb","vare") # add monitor random effect/error variance
Cmodel <- compileNimble(model)


# Compile MCMC Object
modelMCMC <- buildMCMC(modelConf)
modelMCMC <- compileNimble(modelMCMC, project = model)


# Run MCMC
niter <- 10000     # Number of iterations
burn  <- niter/2   # Number of burn-ins
thin  <- 20        # Thinning interval   
set.seed(1)
samples <- runMCMC(modelMCMC, niter = niter,nburnin = burn, nchains = 1,thin = thin,summary = TRUE)


# Extract Sample
samples$summary  # estimate table


# MCMC Diagnostics
# Check MCMC plot
out<-samples$samples[,c("beta0","beta1","taub","varb","taue","vare")]
mcmcplot(out,filename = "output")


