library(mcmcplots)      
library(nimble)
library(ggplot2)

# Generate Data
set.seed(04142022)
n<-500                           # total subjects

# Fixed Effects
t<-sample(seq(10,60,by=0.5),n,replace=T) # time variable

# Set CP occurs at day 30
kappa<-30                       # cp
xg<-(t-kappa)*(t>=kappa)         # spline function

X<-cbind(1,t,xg)                     # design matrix 
truebeta<-beta<-c(.45,-.8,1.0)       # true beta 
p<-ncol(X)                           # number of parameters

# Error variance
truesigmae<-sigmae<-2

# Outcome Model
mu<-X%*%beta
y<-rnorm(n,mu,sqrt(sigmae)) 

# Scatter plot
dat<-data.frame(y=y,t=t)

ggplot(dat,aes(x=t,y=y))+
  geom_point()+ 
  geom_smooth(method='lm',col="blue")+ # linear regression
  geom_smooth(col="red")               # LOESS smooth regression


#--------#
# Nimble #
#--------#
code <- nimbleCode({
  for(i in 1:n) {
    y[i] ~ dnorm(mu[i],tau=taue)
    mu[i] <- beta0+beta1*t[i]+beta2*(t[i]-kappa)*(t[i]>=kappa)
  }
  # Prior 
  taue  ~ dgamma(0.0001,0.0001)            # Error precision
  beta0  ~ dnorm(0, sd = 100)              
  beta1  ~ dnorm(0, sd = 100)
  beta2  ~ dnorm(0, sd = 100)
  vare <- 1/taue                           # Error variance
  kappa ~ T(dnorm(35,sd = 100),10,60)      # kappa ~ T(mu=35,sigma2=100,L0=10,U0=60)
})

constants <- list(n = n)        
data <- list(y = y, t = t)
inits <- list(beta0=0,beta1=0,beta2=0,kappa=35,taue=1) # Init for kappa

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
niter <- 30000     # Number of iterations
burn  <- niter/2   # Number of burn-ins
thin  <- 20        # Thinning interval   
set.seed(1)
samples <- runMCMC(modelMCMC, niter = niter,nburnin = burn, nchains = 1,thin = thin,summary = TRUE)

# Extract Sample
samples$summary  # estimate table

# MCMC Diagnostics
# Check MCMC plot
out<-samples$samples[,c("beta0","beta1","beta2","taue","vare","kappa")]
mcmcplot(out,filename = "output")


# Fitted Lines
mcmc.out<-samples$samples
Beta0<-mcmc.out[,"beta0"]
Beta1<-mcmc.out[,"beta1"]
Beta2<-mcmc.out[,"beta2"]
Kappa<-mcmc.out[,"kappa"]

lastit<-(niter-burn)/thin
ts<-seq(10,60,by=1)              # Create time variable
YPRED<-matrix(0,lastit,length(ts))

# Calcualte the fitted line in each MCMC
for (j in 1:lastit){
  beta0<-Beta0[j]
  beta1<-Beta1[j]
  beta2<-Beta2[j]
  kappa<-Kappa[j]
  
  ypred<-beta0+beta1*ts+beta2*(ts-kappa)*(ts>=kappa)
  
  YPRED[j,]<-ypred
}

# Overall fitted line (95% CrI)
ypred<-colMeans(YPRED)
upper<-apply(YPRED,2,quantile,.975)
lower<-apply(YPRED,2,quantile,.025)

# Selected MCMC draws
nsub<-100                   # Choose 100 MCMC draws 
ysub<-sample(1:lastit,nsub) # Randomly sample 100 MCMC draws

YPRED_SUB<-YPRED[ysub,]

# Dataset for overall fitted line
dplot<-data.frame(t=ts,
                  mu=ypred,
                  lower=lower,
                  upper=upper)

# Dataset for selected MCMC draws
dplot_sub<-data.frame(t=rep(ts,nsub),
                      mu=c(t(YPRED_SUB)),
                      id=rep(1:nsub,each=length(ts)))

ggplot(dplot_sub,aes(x=t,y=mu))+
  geom_line(aes(group=id,col="Selected 100 MCMC Draws"))+
  geom_point(data=dat,aes(x=t,y=y))+
  geom_line(data=dplot,aes(x=t,y=mu,col="Posterior Fitted Line"),linewidth=1)+
  geom_vline(aes(xintercept = mean(Kappa),col="Estiamte CP"))+
  scale_x_continuous(breaks = seq(10,60,by=5))+
  xlab("days")+ylab("Y")+
  scale_color_manual(breaks = c("Selected 100 MCMC Draws","Posterior Fitted Line","Estiamte CP"), 
                     values = c("gray50","blue","red"))+
  guides(color = guide_legend(title="Estiamted Lines",
                              override.aes = list(
                                linetype = c("solid","solid","solid"),
                                linewidth = c(1,1,1),
                                shape=c(NA,NA,NA)),
                              reverse = F), fill="none",shape="none",linetype="none",group="none")+
  theme(legend.key = element_rect(fill = "white"),
        legend.position = c(0.75,0.85),legend.text=element_text(size=15),
        legend.title=element_text(size=15),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.margin=unit(c(0.5, 1, 1, 0.5), "lines"))









