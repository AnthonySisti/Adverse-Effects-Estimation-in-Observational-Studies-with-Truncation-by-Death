library(runjags)
library(rjags)
library(R2jags)
library(mgcv)
library(boot)
library(coda)
library(splines2)
library(ResourceSelection)
library(readxl)

#Read in Data, create variable designating who was assigned treatment and control
Synth.X <- read_xls("your_path_to/Synth_Analy_X.xls")
#Synth.X has 2016 observations and 15 variables, 5 binary
#This is a reduced (columnwise) synthetic data set compared to the real one used in the analysi
#to allow for reasonable computation times
#Reduce columns further as needed when testing if you want reduced computation time,
#or modify the code to use bayes glm as was done in the simulation.
treatment <- c(rep(0, nrow(Synth.X)/2),rep(0, nrow(Synth.X)/2))
X.Sim <- as.matrix(cbind(1,Synth.X))
Nobs<-nrow(X.Sim)



# Define Functions that will compute composite ordinal outcome based on adverse event and death
assign.group = function(heart, death){
  group = vector(length=length(heart))
  for(i in 1:length(heart)){
    if(heart[i]==0 & death[i]==0){
      group[i] = 1}
    if(heart[i]==1 & death[i]==0){
      group[i]=2}
    if(heart[i]==0 & death[i]==1){
      group[i]=3}
    if(heart[i]==1 & death[i]==1){
      group[i]=4} }
  return(group)
}

logit <- function(x) qlogis(x)
inv.logit <- function(x) plogis(x)



### Create Outcome Data for Data Analysis Example ###


# Set model parameters and coefficient values
int.treat.hf <- 1
int.control.hf <- 1
int.treat.death <- 0.5
int.control.death <- 0.75
gamma.treat.true <- 1.5
gamma.control.true <-0.5
# Here we are sampling the model coefficients from a normal distribution
set.seed(1)
beta.treat.true <- c(int.treat.hf , rnorm(ncol(Synth.X),0,0.25))
beta.control.true <- c(int.control.hf, rnorm(ncol(Synth.X),0,0.25)) 
beta2.treat.true <- c(int.treat.death , rnorm(ncol(Synth.X),0,0.25))
beta2.control.true = c(int.control.death , rnorm(ncol(Synth.X),0,0.25))
set.seed(NULL)


# Calculate the true expected probabilites of adverse events under each intervention
true.probs.treat.hf = inv.logit(X.Sim%*%beta.treat.true)
true.probs.control.hf = inv.logit(X.Sim%*%beta.control.true)

# Calculate the true expected probabilities of each level of composite ordinal outcome under treatment
true.probs.treat.yesdeath.yesHF = inv.logit(X.Sim%*%beta2.treat.true+rep(1,nrow(X.Sim))*gamma.treat.true)
true.probs.treat.yesdeath.noHF = inv.logit(X.Sim%*%beta2.treat.true+rep(0,nrow(X.Sim))*gamma.treat.true)
true.probs.treat.nodeath.yesHF  = 1- inv.logit(X.Sim%*%beta2.treat.true+rep(1,nrow(X.Sim))*gamma.treat.true)
true.probs.treat.nodeath.noHF = 1- inv.logit(X.Sim%*%beta2.treat.true+rep(0,nrow(X.Sim))*gamma.treat.true)

true.probs.treat.G1 <- true.probs.treat.nodeath.noHF*(1-true.probs.treat.hf)
true.probs.treat.G2 <- true.probs.treat.nodeath.yesHF*(true.probs.treat.hf)
true.probs.treat.G3 <- true.probs.treat.yesdeath.noHF*(1-true.probs.treat.hf)
true.probs.treat.G4 <- true.probs.treat.yesdeath.yesHF*(true.probs.treat.hf)

expected.treat.G1 <- mean(true.probs.treat.G1)
expected.treat.G2 <- mean(true.probs.treat.G2)
expected.treat.G3 <- mean(true.probs.treat.G3)
expected.treat.G4 <- mean(true.probs.treat.G4)


# Calculate the true expected probabilities of each level of composite ordinal outcome under control
true.probs.control.yesdeath.yesHF = inv.logit(X.Sim%*%beta2.control.true+rep(1,nrow(X.Sim))*gamma.control.true)
true.probs.control.yesdeath.noHF = inv.logit(X.Sim%*%beta2.control.true+rep(0,nrow(X.Sim))*gamma.control.true)
true.probs.control.nodeath.yesHF  = 1- inv.logit(X.Sim%*%beta2.control.true+rep(1,nrow(X.Sim))*gamma.control.true)
true.probs.control.nodeath.noHF = 1- inv.logit(X.Sim%*%beta2.control.true+rep(0,nrow(X.Sim))*gamma.control.true)

true.probs.control.G1 <- true.probs.control.nodeath.noHF*(1-true.probs.control.hf)
true.probs.control.G2 <- true.probs.control.nodeath.yesHF*(true.probs.control.hf)
true.probs.control.G3 <- true.probs.control.yesdeath.noHF*(1-true.probs.control.hf)
true.probs.control.G4 <- true.probs.control.yesdeath.yesHF*(true.probs.control.hf)

expected.control.G1 <- mean(true.probs.control.G1)
expected.control.G2 <- mean(true.probs.control.G2)
expected.control.G3 <- mean(true.probs.control.G3)
expected.control.G4 <- mean(true.probs.control.G4)


#Calculate values of traditional estimands
trueExpITTheart <- ((expected.treat.G2+expected.treat.G4) - (expected.control.G2+expected.control.G4))
trueExpITTdeath <- ((expected.treat.G3+expected.treat.G4)) - (expected.control.G3+expected.control.G4)
trueExpCompBin <- ((expected.treat.G2+expected.treat.G3+expected.treat.G4)- (expected.control.G2+expected.control.G3+expected.control.G4))
trueExpSACE <- ((expected.treat.G2)/(expected.treat.G1+expected.treat.G2)) - ((expected.control.G2)/(expected.control.G1+expected.control.G2))


#Calculate values of composite ordinal estimands
expected.diff.G1 <- expected.treat.G1-expected.control.G1
expected.diff.G2 <- expected.treat.G2-expected.control.G2
expected.diff.G3 <- expected.treat.G3-expected.control.G3
expected.diff.G4 <- expected.treat.G4-expected.control.G4

true_Pr_G1_gr_G0 <- true.probs.control.G1*(true.probs.treat.G2+true.probs.treat.G3+true.probs.treat.G4)+
  true.probs.control.G2*(true.probs.treat.G3+true.probs.treat.G4)+
  true.probs.control.G3*(true.probs.treat.G4)

true_Pr_G0_gr_G1 <- true.probs.treat.G1*(true.probs.control.G2+true.probs.control.G3+true.probs.control.G4)+
  true.probs.treat.G2*(true.probs.control.G3+true.probs.control.G4)+
  true.probs.treat.G3*(true.probs.control.G4)


trueExpDiff_k10k01 <- colMeans(true_Pr_G1_gr_G0)-colMeans(true_Pr_G0_gr_G1)
trueExpRR_k10k01 <- colMeans(true_Pr_G1_gr_G0)/colMeans(true_Pr_G0_gr_G1)


# View true expected values of each estimand for this study
trueExpITTheart
trueExpITTdeath
trueExpCompBin
trueExpSACE

expected.diff.G1
expected.diff.G2
expected.diff.G3
expected.diff.G4

trueExpDiff_k10k01
trueExpRR_k10k01


# Define the propensity score model

true.prop.score  <- inv.logit(as.matrix(X.Sim)%*%as.matrix(rnorm(ncol(X.Sim),0,0.25)))


# sample  treatment and control adverse events and death based on true probabilities
samp.treat.hf <- rbinom(length(true.probs.treat.hf),1,true.probs.treat.hf)
samp.control.hf <- rbinom(length(true.probs.control.hf),1,true.probs.control.hf)

true.probs.treat.death <- inv.logit(X.Sim%*%beta2.treat.true + samp.treat.hf*gamma.treat.true)
true.probs.control.death <- inv.logit(X.Sim%*%beta2.control.true + samp.control.hf*gamma.control.true)

samp.treat.death <- rbinom(length(true.probs.treat.death),1,true.probs.treat.death)
samp.control.death <- rbinom(length(true.probs.control.death),1,true.probs.control.death)


# define observed and unobserved adverse event and death outcomes based on treatment assignment
Wobs <- rbinom(Nobs,1,prob =  true.prop.score) 
TreatAssign <- Wobs==1

obs.treat.hf <- samp.treat.hf[TreatAssign]
obs.control.hf <- samp.control.hf[!TreatAssign]

mis.treat.hf <- samp.treat.hf[!TreatAssign]
mis.control.hf <- samp.control.hf[TreatAssign]

obs.treat.death <- samp.treat.death[TreatAssign]
obs.control.death <- samp.control.death[!TreatAssign]

mis.treat.death <- samp.treat.death[!TreatAssign]
mis.control.death <- samp.control.death[TreatAssign]

samp.treat.group <- assign.group(samp.treat.hf,samp.treat.death)
samp.control.group <- assign.group(samp.control.hf,samp.control.death)

obs.treat.group <- assign.group(obs.treat.hf,obs.treat.death)
obs.control.group <- assign.group(obs.control.hf,obs.control.death)




### Conducting the Analysis ###

# Estimate propensity score based on generated treatment assignment
X.P <- as.data.frame(X.Sim[,-1])
X.P$W <- Wobs
prop_fit <- glm(W~., data=X.P, family = binomial())
est.prop <- predict(prop_fit, newdata=X.P, type="response")


# Construct Spline on the estimated propensity score
# The algorithm below can help ensure there are 3 treatment and
# control individuals in each subclass

check.hf <- samp.treat.hf
check.hf[!TreatAssign]<-samp.control.hf[!TreatAssign]
check.death <- samp.treat.death
check.death[!TreatAssign]<-samp.control.death[!TreatAssign]
subclasses <-7
flag_knots <- TRUE
while (flag_knots==TRUE & subclasses>=4) {
  subclasses <-subclasses-1
  use_knots<- subclasses-1
  knots = quantile(logit(est.prop),(1:(subclasses-1)/subclasses))
  
  check_knots <- c(min(logit(est.prop)),knots,max(logit(est.prop)))
  check_sum_hf_treat <- rep(NA,subclasses)
  check_sum_hf_control <- rep(NA,subclasses)
  check_sum_death_treat <- rep(NA,subclasses)
  check_sum_death_control <- rep(NA,subclasses)
  for (q in 1:subclasses) {
    # q=1
    check_sum_hf_control[q] <- sum(check.hf[!TreatAssign & logit(est.prop) >= check_knots[q] & logit(est.prop) <= check_knots[q+1]])
    check_sum_hf_treat[q] <- sum(check.hf[TreatAssign & logit(est.prop) >= check_knots[q] & logit(est.prop) <= check_knots[q+1]])
    
    check_sum_death_control[q]<- sum(check.death[!TreatAssign & logit(est.prop) >= check_knots[q] & logit(est.prop) <= check_knots[q+1]])
    check_sum_death_treat[q] <- sum(check.death[TreatAssign & logit(est.prop) >= check_knots[q] & logit(est.prop) <= check_knots[q+1]])
  }
  
  flag_knots <- any(c(check_sum_hf_control,
                      check_sum_hf_treat,
                      check_sum_death_treat,
                      check_sum_death_control)<3)
  
}



b.prop = naturalSpline(logit(est.prop),
                    knots= knots, intercept=FALSE)

b.treat <- b.prop[Wobs==1,]
b.control <- b.prop[Wobs==0,]


#Add spline to data set containing covariates. Remove intercept and one covariate

X.treat = cbind(b.treat, X.Sim[Wobs==1,-c(1,2)])
X.control = cbind(b.control,X.Sim[Wobs==0,-c(1,2)])



# Create Data lists in order to fit Bayesian model in JAGS

datalist.treat=list(X=as.matrix(X.treat), death=obs.treat.death, s=ncol(b.treat),
                    heart = obs.treat.hf, n=dim(X.treat)[1], m=dim(X.treat)[2])
datalist.control=list(X=as.matrix(X.control), death=obs.control.death,  s=ncol(b.control),
                    heart = obs.control.hf, n=dim(X.control)[1], m=dim(X.control)[2])

n=dim(X.treat)[1] 
m=dim(X.treat)[2]


#Define Model for use with JAGS, we can use a non-informative prior or shrinkage prior


model.string.normal="
model{

eta.heart = X%*%beta
for(i in 1:n){
prob.heart[i]=ilogit(eta.heart[i])
heart[i]~dbern(prob.heart[i])
}

eta.death = X%*%beta2 + gamma*heart
for(i in 1:n){
prob.death[i]=ilogit(eta.death[i])
death[i]~dbern(prob.death[i])
}

gamma ~ dnorm(0,1/100)

for(i in 1:s){beta[i] ~ dnorm(0, 1/100)
              beta2[i] ~ dnorm(0,1/100)}
              
for(i in (s+1):m){beta[i] ~ dnorm(0, 1/100)
              beta2[i] ~ dnorm(0,1/100)}
}

"





model.string.shrink="
model{

eta.heart = X%*%beta
for(i in 1:n){
prob.heart[i]=ilogit(eta.heart[i])
heart[i]~dbern(prob.heart[i])
}

eta.death = X%*%beta2 + gamma*heart
for(i in 1:n){
prob.death[i]=ilogit(eta.death[i])
death[i]~dbern(prob.death[i])
}


r.a.c~dt(0, 1, 1)T(0,)
r.d.c~dt(0, 1, 1)T(0,)

gamma ~ dnorm(0,1/100)

for(i in 1:s){beta[i] ~ dnorm(0, 1/100)
              beta2[i] ~ dnorm(0,1/100)}

for(i in (s+1):m){beta[i] ~ dnorm(0, r.a.c/100)
              beta2[i] ~ dnorm(0,r.d.c/100)}
}

"


# Run JAGS to fit the Bayesian model on the control observations and the treatment observations
# The MCMC samples are highly correlated due to spline, thin needs to be large

out.control = run.jags(model.string.normal, data=datalist.control, 
                       monitor=c("beta","beta2", "gamma"),
                       burnin=500,adapt=100, sample=100, thin=50, n.chains=3)

out.treat = run.jags(model.string.normal, data=datalist.treat, 
                     monitor=c("beta","beta2","gamma"),
                     burnin=500,adapt=100, sample=100, thin=50, n.chains=3)

control.summary <- as.data.frame(summary(out.control))
treat.summary <- as.data.frame(summary(out.treat))


#Check Convergence, adjust thin as needed
hist(control.summary$psrf)
control.summary[order(control.summary$psrf,decreasing = T),]
hist(treat.summary$psrf)
treat.summary[order(treat.summary$psrf,decreasing = T),]

#Trace Plots
res.control = as.matrix(out.control$mcmc)
res.treat = as.matrix(out.treat$mcmc)
plot(1:nrow(res.control), res.control[,colnames(res.control) == "beta2[10]"], type = "l")
plot(1:nrow(res.treat), res.treat[,colnames(res.treat) == "beta[10]"], type = "l")


#Super-population causal inference, using the empirical distribution
#of covariates to approximate the population distribution
beta.treat = t(res.treat[,1:m])
beta.control = t(res.control[,1:m])

beta2.treat = t(res.treat[,(m+1):(2*m)])
beta2.control = t(res.control[, (m+1):(2*m)])

gamma.treat = res.treat[, 2*m+1]
gamma.control = res.control[, 2*m+1]

X.Est <- rbind(X.treat,X.control)

samp.probs.treat.hf <-  inv.logit(X.Est%*%beta.treat)

samp.probs.treat.yesdeath.yesHF = inv.logit(X.Est%*%beta2.treat+matrix(1,nrow = nrow(X.Est), ncol =length(gamma.treat))%*%diag(gamma.treat))
samp.probs.treat.yesdeath.noHF = inv.logit(X.Est%*%beta2.treat)
samp.probs.treat.nodeath.yesHF  = 1- samp.probs.treat.yesdeath.yesHF
samp.probs.treat.nodeath.noHF = 1- samp.probs.treat.yesdeath.noHF 

samp.probs.treat.G1 <- samp.probs.treat.nodeath.noHF*(1-samp.probs.treat.hf)
samp.probs.treat.G2 <- samp.probs.treat.nodeath.yesHF*(samp.probs.treat.hf)
samp.probs.treat.G3 <- samp.probs.treat.yesdeath.noHF*(1-samp.probs.treat.hf)
samp.probs.treat.G4 <- samp.probs.treat.yesdeath.yesHF*(samp.probs.treat.hf)

samp.probs.control.hf<-  inv.logit(X.Est%*%beta.control)

samp.probs.control.yesdeath.yesHF = inv.logit(X.Est%*%beta2.control+matrix(1,nrow = nrow(X.Est), ncol =length(gamma.control))%*%diag(gamma.control))
samp.probs.control.yesdeath.noHF = inv.logit(X.Est%*%beta2.control)
samp.probs.control.nodeath.yesHF  = 1- samp.probs.control.yesdeath.yesHF
samp.probs.control.nodeath.noHF = 1- samp.probs.control.yesdeath.noHF 

samp.probs.control.G1 <- samp.probs.control.nodeath.noHF*(1-samp.probs.control.hf)
samp.probs.control.G2 <- samp.probs.control.nodeath.yesHF*(samp.probs.control.hf)
samp.probs.control.G3 <- samp.probs.control.yesdeath.noHF*(1-samp.probs.control.hf)
samp.probs.control.G4 <- samp.probs.control.yesdeath.yesHF*(samp.probs.control.hf)


#Estimate P(G(1)>G(0)) and P(G(0)>G(1)) and 
Pr_G1_gr_G0 <- samp.probs.control.G1*(samp.probs.treat.G2+samp.probs.treat.G3+samp.probs.treat.G4)+
  samp.probs.control.G2*(samp.probs.treat.G3+samp.probs.treat.G4)+
  samp.probs.control.G3*(samp.probs.treat.G4)

Pr_G0_gr_G1 <- samp.probs.treat.G1*(samp.probs.control.G2+samp.probs.control.G3+samp.probs.control.G4)+
  samp.probs.treat.G2*(samp.probs.control.G3+samp.probs.control.G4)+
  samp.probs.treat.G3*(samp.probs.control.G4)

k10_k01_diff <- colMeans(Pr_G1_gr_G0 -Pr_G0_gr_G1)

k10_k01_RR <- colMeans(Pr_G1_gr_G0)/colMeans(Pr_G0_gr_G1)


#Distribution of adverse event ITT effect
treat.hf.samp <-colMeans(samp.probs.treat.G2)+colMeans(samp.probs.treat.G4)
mean(treat.hf.samp)
quantile(treat.hf.samp,c(0.025,0.975))

control.hf.samp <-colMeans(samp.probs.control.G2)+colMeans(samp.probs.control.G4)
mean(control.hf.samp)
quantile(control.hf.samp,c(0.025,0.975))

mean(treat.hf.samp-control.hf.samp)
quantile(treat.hf.samp-control.hf.samp, c(0.025, 0.975))


#Distribution of death ITT effect
treat.death.samp <-colMeans(samp.probs.treat.G3)+colMeans(samp.probs.treat.G4)
mean(treat.death.samp)
quantile(treat.death.samp,c(0.025,0.975))

control.death.samp <-colMeans(samp.probs.control.G3)+colMeans(samp.probs.control.G4)
mean(control.death.samp)
quantile(control.death.samp,c(0.025,0.975))

mean(treat.death.samp-control.death.samp)
quantile(treat.death.samp-control.death.samp, c(0.025, 0.975))


#Distribution of complosite binary ITT effect
treat.compBin.samp <-colMeans(samp.probs.treat.G2)+colMeans(samp.probs.treat.G4)+colMeans(samp.probs.control.G3)
mean(treat.compBin.samp)
quantile(treat.compBin.samp,c(0.025,0.975))

control.compBin.samp <-colMeans(samp.probs.control.G2)+colMeans(samp.probs.control.G4)+colMeans(samp.probs.control.G3)
mean(control.compBin.samp)
quantile(control.compBin.samp,c(0.025,0.975))

mean(treat.compBin.samp-control.compBin.samp)
quantile(treat.compBin.samp-control.compBin.samp, c(0.025, 0.975))


#Distribution of SACE estimand
treat.SACE.samp <-colMeans(samp.probs.treat.G2)/colMeans(samp.probs.treat.G2+samp.probs.treat.G1)
mean(treat.SACE.samp)
quantile(treat.SACE.samp,c(0.025,0.975))

control.SACE.samp <-colMeans(samp.probs.control.G2)/colMeans(samp.probs.control.G2+samp.probs.control.G1)
mean(control.SACE.samp)
quantile(control.SACE.samp,c(0.025,0.975))

mean(treat.SACE.samp-control.SACE.samp)
quantile(treat.SACE.samp-control.SACE.samp, c(0.025, 0.975))


#Distribution of level 3 of composite ordinal outcome difference estimand
G3.diff.samp <- colMeans(samp.probs.treat.G3-samp.probs.control.G3)
mean(G3.diff.samp)
quantile(G3.diff.samp, c(0.025, 0.975))


#Distribution of \kappa_10 - \kappa_01 estimand
mean(k10_k01_diff)
quantile(k10_k01_diff, c(0.025,0.975))





#Finite Sample Estimands Multiply Imputing unobserved poteintial outcomes

A_1 <- matrix(NA, nrow = nrow(X.Est), ncol = ncol(beta.treat))
A_0 <- matrix(NA, nrow = nrow(X.Est), ncol = ncol(beta.treat))
D_1 <- matrix(NA, nrow = nrow(X.Est), ncol = ncol(beta.treat))
D_0 <- matrix(NA, nrow = nrow(X.Est), ncol = ncol(beta.treat))
G_1 <- matrix(NA, nrow = nrow(X.Est), ncol = ncol(beta.treat))
G_0 <- matrix(NA, nrow = nrow(X.Est), ncol = ncol(beta.treat))

n_treat <- sum(Wobs==1)
n_control <- sum(Wobs==0)
samp.logit.treat.death <- X.Est%*%beta2.treat
samp.logit.control.death <- X.Est%*%beta2.control

for (j in 1:ncol(A_1)) {
  A_1[Wobs==1,j] <- obs.treat.hf
  A_1[Wobs==0,j] <- rbinom(n_control, 1, samp.probs.treat.hf[Wobs==0,j])
  
  A_0[Wobs==0,j] <- obs.control.hf
  A_0[Wobs==1,j] <- rbinom(n_treat, 1, samp.probs.control.hf[Wobs==1,j])
  
  
  D_1[Wobs==1,j] <- obs.treat.death
  finite.treat.cf.probs.death <- inv.logit(samp.logit.treat.death[Wobs==0,j] + gamma.treat[j]*A_1[Wobs==0,j])
  D_1[Wobs==0,j] <- rbinom(n_control, 1, finite.treat.cf.probs.death)
  
  D_0[Wobs==0,j] <- obs.control.death
  finite.control.cf.probs.death <- inv.logit(samp.logit.control.death[Wobs==1,j] + gamma.control[j]*A_0[Wobs==1,j])
  D_0[Wobs==1,j] <- rbinom(n_treat, 1, finite.control.cf.probs.death)
  
  G_1[,j] <- assign.group(A_1[,j],D_1[,j])
  G_0[,j] <- assign.group(A_0[,j],D_0[,j])
}


# Adverse event ITT finite sample estimand

fs.ITT.hf <-  colMeans(A_1-A_0)
mean(fs.ITT.hf)
quantile(fs.ITT.hf,c(0.025,0.975))


# Death ITT finite sample estimand

fs.ITT.death <-  colMeans(D_1-D_0)
mean(fs.ITT.death)
quantile(fs.ITT.death,c(0.025,0.975))


# kappa_10 - kappa_01  finite sample estimand

fs.k10_k01_diff <- colMeans(G_1>G_0)-colMeans(G_0>G_1)
mean(fs.k10_k01_diff)
quantile(fs.k10_k01_diff,c(0.025,0.975))





### Multiple Imputation for Population Estimand Example #####

#Difference in heart failure

M <- ncol(A_1)
N <- nrow(A_1)

sp.ITT.hf <-  mean(colMeans(A_1-A_0))
U_bar <- mean(apply(A_1-A_0,2,var)/N)
B <- 1/(M-1) * sum((colMeans(A_1-A_0)-sp.ITT.hf)^2)
TotV<- U_bar + (1+1/M)*B
t_df <- (M-1)*(1+U_bar/((1+1/M)*B))

#estimate
sp.ITT.hf 
# 95% interval
c(sp.ITT.hf - sqrt(TotV)*qt(0.975,t_df), sp.ITT.hf + sqrt(TotV)*qt(0.975,t_df))

#Compare to Bayesian credible interval for analytically derived super-population estimand
quantile(treat.hf.samp-control.hf.samp, c(0.025, 0.975))

#Compare to finite sample
quantile(fs.ITT.hf,c(0.025,0.975))



#Difference in death

M <- ncol(D_1)
N <- nrow(D_1)

sp.ITT.death <-  mean(colMeans(D_1-D_0))
U_bar <- mean(apply(D_1-D_0,2,var)/N)
B <- 1/(M-1) * sum((colMeans(D_1-D_0)-sp.ITT.death)^2)
TotV<- U_bar + (1+1/M)*B
t_df <- (M-1)*(1+U_bar/((1+1/M)*B))

#estimate
sp.ITT.death 
# 95% interval
c(sp.ITT.death - sqrt(TotV)*qt(0.975,t_df),sp.ITT.death + sqrt(TotV)*qt(0.975,t_df))

#Compare to Bayesian credible interval for analytically derived super-population estimand
quantile(treat.death.samp-control.death.samp, c(0.025, 0.975))

#Compare to finite sample
quantile(fs.ITT.death,c(0.025,0.975))




##### Sensitivity Analysis for super population estimands ####

#First we define a function that includes an unobserved covariate in the estimation procedure
#The estimand of interest is standardized the relative treatment effect (kappa_10 - kappa01)/SE(kappa_10 - kappa01)

Sensitivity <- function(delta, # vector containing 2 parameters 
                        mu_Z, # standardized bias between treatment and control
                        X.t =X.treat, X.c=X.control,  # data matricies
                        res.t=res.treat, res.c=res.control # output of jags models containing posterior samples of model parameters
                        ){
  delta_a <- delta[1]
  delta_d <- delta[2]
  mu_Z_1 <- 0
  mu_Z_0 <- mu_Z
  
  X.control <- X.c
  X.treat <- X.t
  res.treat <- res.t
  res.control <- res.c
  
  n_0 <- nrow(X.control)
  n_1 <- nrow(X.treat)  
  
  Z_0 <- rnorm(n_0,mu_Z_0)
  Z_1 <- rnorm(n_1,mu_Z_1)
  
  X.treat.Est <- cbind(X.treat,Z_1)
  X.control.Est <- cbind(X.control,Z_0)
  X.Est <- rbind(X.treat.Est,X.control.Est)
  
  beta.treat = rbind(t(res.treat[,1:m]),delta_a)
  beta.control = rbind(t(res.control[,1:m]),delta_a)
  
  beta2.treat = rbind(t(res.treat[,(m+1):(2*m)]),delta_d)
  beta2.control = rbind(t(res.control[, (m+1):(2*m)]),delta_d)
  
  gamma.treat = res.treat[, 2*m+1]
  gamma.control = res.control[, 2*m+1]
  
  samp.probs.treat.hf <-  inv.logit(X.Est%*%beta.treat)
  XB2_treat <-X.Est%*%beta2.treat
  
  samp.probs.treat.yesdeath.yesHF = inv.logit(XB2_treat+matrix(1,nrow = nrow(X.Est), ncol =length(gamma.treat))%*%diag(gamma.treat))
  samp.probs.treat.yesdeath.noHF = inv.logit(XB2_treat)
  samp.probs.treat.nodeath.yesHF  = 1- samp.probs.treat.yesdeath.yesHF
  samp.probs.treat.nodeath.noHF = 1- samp.probs.treat.yesdeath.noHF 
  
  samp.probs.treat.G1 <- samp.probs.treat.nodeath.noHF*(1-samp.probs.treat.hf)
  samp.probs.treat.G2 <- samp.probs.treat.nodeath.yesHF*(samp.probs.treat.hf)
  samp.probs.treat.G3 <- samp.probs.treat.yesdeath.noHF*(1-samp.probs.treat.hf)
  samp.probs.treat.G4 <- samp.probs.treat.yesdeath.yesHF*(samp.probs.treat.hf)
  
  samp.probs.control.hf<-  inv.logit(X.Est%*%beta.control)
  XB2_control <-X.Est%*%beta2.control
  
  samp.probs.control.yesdeath.yesHF = inv.logit(XB2_control+matrix(1,nrow = nrow(X.Est), ncol =length(gamma.control))%*%diag(gamma.control))
  samp.probs.control.yesdeath.noHF = inv.logit(XB2_control)
  samp.probs.control.nodeath.yesHF  = 1- samp.probs.control.yesdeath.yesHF
  samp.probs.control.nodeath.noHF = 1- samp.probs.control.yesdeath.noHF 
  
  samp.probs.control.G1 <- samp.probs.control.nodeath.noHF*(1-samp.probs.control.hf)
  samp.probs.control.G2 <- samp.probs.control.nodeath.yesHF*(samp.probs.control.hf)
  samp.probs.control.G3 <- samp.probs.control.yesdeath.noHF*(1-samp.probs.control.hf)
  samp.probs.control.G4 <- samp.probs.control.yesdeath.yesHF*(samp.probs.control.hf)

  Pr_G1_gr_G0 <- samp.probs.control.G1*(samp.probs.treat.G2+samp.probs.treat.G3+samp.probs.treat.G4)+
    samp.probs.control.G2*(samp.probs.treat.G3+samp.probs.treat.G4)+
    samp.probs.control.G3*(samp.probs.treat.G4)
  
  Pr_G0_gr_G1 <- samp.probs.treat.G1*(samp.probs.control.G2+samp.probs.control.G3+samp.probs.control.G4)+
    samp.probs.treat.G2*(samp.probs.control.G3+samp.probs.control.G4)+
    samp.probs.treat.G3*(samp.probs.control.G4)
  
  k10_k01_diff <- colMeans(Pr_G1_gr_G0 -Pr_G0_gr_G1)
  
  return(mean(k10_k01_diff)/sd(k10_k01_diff))
}


# We define the standaridize bias and grid of delta values for the conditional log odds ratio
# for adverse event and death that we want to investigate

B <- -1
delta <- cbind(rep(seq(-1,1,length.out =25) ,25),rep(seq(-1,1,length.out =25) ,each=25))


# We apply the sensitivity function over the grid of delta values at our desired standardized bias
# This can take some time
k10_k01_diff_Sens <- apply(delta,1, Sensitivity, mu_Z =B)


# We now plot the results

sens_data_plot <- data.frame(delta_a = delta[,1],
                             delta_d = delta[,2],
                             estimand = k10_k01_diff_Sens)


library(ggplot2)
library(latex2exp)


# Estimated Stadardized Relative Treatment Effect:
est.val <- mean(k10_k01_diff)/sd(k10_k01_diff)

# Plot the values of the esitmand for different combinations of delta values
# We may be interested if the effect becomes or drops below statsitical significance or how the effect changes with differnt delta in general

ggplot(data =sens_data_plot, aes(x=delta_a, y=delta_d))+
  geom_raster(aes( fill= estimand), interpolate=TRUE)+
  xlab(TeX("$\\delta_a$"))+labs(fill = "RTE/SE") +
  ylab(TeX("$\\delta_d$"))+
  scale_fill_gradient2(low="blue", mid="white", high="red", # This defines the heat map, the range of colors that will display the change in estimate: red means higher, blue means lower
                       midpoint=est.val, limits=c(1,10) # We can set midpoint to be 0, or to the estimated value, adjust limits to consider range of interest
                       ) +
  theme_bw()+ggtitle(TeX(r"(Sensitivity of Standardized RTE: $\mu^z_1$ = 1)"))+
  geom_point(aes(x=0,y=0),colour="black")


#We can also plot a histogram of the standardized RTE over all combinations
# of delta values we specified in the funciton

ggplot(data =sens_data_plot, aes(x=estimand))+
  geom_histogram(color="black", fill= "steelblue", bins=30)+
  labs(title = TeX(r"(Range of Standardized RTE when $\mu^z_1$ = 1 for all $\delta_a$ and $\delta_d$ pairs)"),
       x=TeX(r"(RTE/SE )"))+theme_bw()

