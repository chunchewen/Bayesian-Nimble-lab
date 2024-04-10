# Spatial Effect 
# South Carolina County Adjacency Matrix

#Normal Regression
library(nimble, warn.conflicts = FALSE)
library(mcmcplots)
library(PICBayes) # for SC adjacency matrix
library(spam)
library(tmap)
library(tigris)
library(RColorBrewer)
#------------------------------------------------------------------------------#
rm(list=ls())
#------------------------------------------------------------------------------#
set.seed(032210)
n<-46 				                   # Number of county
nis<-sample(1:250,n,T)           # Number of patient in each county
id<-rep(1:n,nis)                 # County ID
N<-length(id)                    # Number of total patients  

# Covariate
x1<-rnorm(N)
x2<-rbinom(N,1,0.5)
X<-cbind(1,x1,x2)
p<-ncol(X)

# Fixed Effect Parms
beta<-c(-0.25,0.65,-0.50)

# Spatial random effect
data(C)
A<-C                                            # Adjancency matrix
D<-apply(A,1,sum)	                              # No. neighbors
Q<-as.spam(diag(D))-as.spam(A) + diag(.0001,n)

varb<-truevarb<-2
taub<-truetaub<-1/varb                          # precision in ICAR
covb<-varb*solve(Q)      
b<-c(mvtnorm::rmvnorm(1,sigma=covb))                 
trueb<-b-mean(b)                               # Center for identifiability 

# Error variance
sigmae<-4

# Response
mu<-X%*%beta+rep(b-mean(b),nis)
y<-rnorm(N,mu,sqrt(sigmae)) 

#--------#
# Nimble #
#--------#
# Adjacency matrix for spatial effect
tmp<-as.carAdjacency(A)
adj<-tmp$adj              # Adjacency vectors
weights<-tmp$weights      # Weights
num<-tmp$num              # Number of neighbors for each county
L<-length(adj)

#--------#
# Nimble #
#--------#
code <- nimbleCode({
  for(i in 1:N) {
    y[i] ~ dnorm(mu[i], tau = taue)
    mu[i]<-inprod(beta[1:p],X[i,1:p])+b[id[i]]
  }
  for (k in 1:p){beta[k]~dnorm(0,sd=100)}  
  b[1:n] ~ dcar_normal(adj=adj[1:L],weights=weights[1:L],
                       num=num[1:n],tau=taub, 
                       zero_mean=1)           # ICAR prior
  # zero_mean=1: zero-mean constraint
  taub  ~ dgamma(0.0001,0.0001)
  varb <- 1/taub 
  taue  ~ dgamma(0.0001,0.0001)
  vare <- 1/taue
})
constants <- list(N=N, n = n, p=p, L=L, id=id,adj=adj,weights=weights,num=num)
data <- list(y = y, X=X)
inits <- list(beta=rep(0,p),taub=1,taue=1, b=rnorm(n))
model <- nimbleModel(code, constants = constants, data = data, inits = inits)

# Configure and Compile Model
modelConf <- configureMCMC(model, print = TRUE)
modelConf$addMonitors("varb","vare","b")
Cmodel <- compileNimble(model)


# Compile MCMC Object
modelMCMC <- buildMCMC(modelConf)
modelMCMC <- compileNimble(modelMCMC, project = model)


# Run MCMC
niter <- 10000     # Number of iterations
burn  <- niter/2   # Number of burn-ins
thin  <- 5         # Thinning interval   
set.seed(1)
samples <- runMCMC(modelMCMC, niter = niter,nburnin = burn, nchains = 1,thin = thin,summary = TRUE)


# Extract Sample
samples$summary[47:53,]  # estimate table

# MCMC Diagnostics
out<-samples$samples[,c("beta[1]","beta[2]","beta[3]","taue","vare","varb","taub")]
mcmcplot(out,filename = "output")

# Observed and Expected Map
dat<-data.frame(countyid=1:n,
                ymean=tapply(y,id,mean))
mb<-apply(samples$samples[,1:46],2,mean)
dat$mb<-mb

scfips<-c("45001","45003","45005","45007","45009",
          "45011","45013","45015","45017","45019","45021",
          "45023","45025","45027","45029","45031","45033",
          "45035","45037","45039","45041","45043","45045",
          "45047","45049","45051","45053","45055","45057",
          "45059","45061","45063","45065","45067","45069",
          "45071","45073","45075","45077","45079","45081",
          "45083","45085","45087","45089","45091")

dat$fips<-scfips
dat$GEOID<-as.character(dat$fips)
sc_counties_map <- counties(state = c('South Carolina'))
sc_data <- left_join(sc_counties_map, dat, by = 'GEOID')
pal <- brewer.pal(5,"BuGn")  # specify the palette colors

tm_shape(sc_data)+
  tm_fill(c("ymean"), midpoint = c(NA), title = c(expression(paste(""))), palette = pal, style = "quantile")+
  tm_facets(free.scales = TRUE, nrow = 1)+
  tm_layout(title = "Quintile Map: Observed Mean Outcome",
            title.snap.to.legend = TRUE,
            title.size = 0.8,
            title.position = c("right", "bottom"),
            legend.outside = FALSE,
            legend.position = c(0.75, 0.88),
            legend.text.size = 0.5,
            main.title.fontface = "bold",
            main.title.position = "center")+
  tm_borders(alpha = 0.3, lwd = 1)

tm_shape(sc_data)+
  tm_fill(c("mb"), midpoint = c(NA), title = c(expression(paste(""))), palette = pal, style = "quantile")+
  tm_facets(free.scales = TRUE, nrow = 1)+
  tm_layout(title = "Quintile Map: Posterior Random Effects",
            title.snap.to.legend = TRUE,
            title.size = 0.8,
            title.position = c("right", "bottom"),
            legend.outside = FALSE,
            legend.position = c(0.75, 0.88),
            legend.text.size = 0.5,
            main.title.fontface = "bold",
            main.title.position = "center")+
  tm_borders(alpha = 0.3, lwd = 1)





