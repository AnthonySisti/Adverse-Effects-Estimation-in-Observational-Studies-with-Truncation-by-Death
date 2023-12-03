library(runjags)
library(rjags)
library(R2jags)
library(mgcv)
library(coda)
library(splines2)
library(ResourceSelection)
library(CBPS)
library(MASS)
library(readxl)


#Read in Data, create variable designating who was assigned treatment and control
Synth.X <- read.xla("your_path_to/Synth_CaseStudy_X.xls")
treatment <- c(rep(0, nrow(Synth.X)/2),rep(0, nrow(Synth.X)/2))
X.Sim <- as.matrix(cbind(1,Synth.X))
Nobs<-nrow(X.Sim)


# Set model parameters and coefficient values
int.treat.hf <- 1
int.control.hf <- 1
int.treat.death <- 0.5
int.control.death <- 0.75
gamma.treat.true <- 1.5
gamma.control.true <-0.5
beta.treat.true <- c(int.treat.hf , 0.21,  0.32,  0.21, -2.33, -1.44) 
beta.control.true <- c(int.control.hf ,  0.21,  0.07,  0.21, -2.62, -2.00) 
beta2.treat.true <- c(int.treat.death, -0.03,  0.18, -0.13, -1.08, -0.19)
beta2.control.true = c(int.control.death, 0.20, 0.22,  0.08, -1.08, -0.55)


# Define link function for generating the propensity score and outcomes.
# When c=1, the burr.link function corresponds to the logistic link function
logit <- function(x) qlogis(x)
inv.logit <- function(x) plogis(x)
burr.link <- function(x,c=1){
  return(1-(1+exp(x))^(-1*c))
}



# Calculate the true expected probabilites of adverse events under each intervention
true.probs.treat.hf = burr.link(X.Sim%*%beta.treat.true)
true.probs.control.hf = burr.link(X.Sim%*%beta.control.true)

# Calculate the true expected probabilities of each level of composite ordinal outcome under treatment
true.probs.treat.yesdeath.yesHF = burr.link(X.Sim%*%beta2.treat.true+rep(1,nrow(X.Sim))*gamma.treat.true)
true.probs.treat.yesdeath.noHF = burr.link(X.Sim%*%beta2.treat.true+rep(0,nrow(X.Sim))*gamma.treat.true)
true.probs.treat.nodeath.yesHF  = 1- burr.link(X.Sim%*%beta2.treat.true+rep(1,nrow(X.Sim))*gamma.treat.true)
true.probs.treat.nodeath.noHF = 1- burr.link(X.Sim%*%beta2.treat.true+rep(0,nrow(X.Sim))*gamma.treat.true)

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



# View true values of each estimand for this case study
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


# Define Functions that will compute composite ordinal outcome and propensity score in simulation 
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


Prop_Scores <- function(W=Wobs, X=X.Sim[,-1]){
  X.P <- as.data.frame(X)
  X.P$W <- W
  prop_fit <- glm(W~., data=X.P, family = binomial())
  p.scores <- predict(prop_fit, newdata=X.P, type="response")
  return(p.scores)
}

DR_Est <- function(Yobs,X=ordered.X,W= ordered.treat,prop=ordered.pscore){

  out_est_t <- predict(suppressWarnings(glm(Yobs~.,data=X,family = binomial(), subset = W==1)), 
                       newdata = X, type="response")
  out_est_c <- predict(suppressWarnings(glm(Yobs~.,data=X,family = binomial(), subset = W==0)),
                       newdata = X,type="response")
  DR <- mean((W*Yobs-(W-prop)*out_est_t)/prop)-  mean(((1-W)*Yobs+(W-prop)*out_est_c)/(1-prop))
  
  I_i <- ((W*Yobs-(W-prop)*out_est_t)/prop - ((1-W)*Yobs+(W-prop)*out_est_c)/(1-prop)) - DR
  
  DR_SE <- sqrt(length(I_i)^(-2) * sum(I_i^2))
  
  DR_L <- DR - 1.96*DR_SE
  DR_U <- DR + 1.96*DR_SE
  return(c(DR_L, DR,DR_U))}



#Set Number of Simulations
Nsim <- 300

#Define matricies to collect simulation information
pp.group.control <- matrix(NA, nrow=Nobs, ncol = 1500)
pp.group.treat <- matrix(NA, nrow=Nobs, ncol = 1500)
checkCover <- matrix(NA, nrow = 10, ncol = Nsim)
checkInterval <- matrix(NA, nrow = 10, ncol =Nsim)
checkMSE <- matrix(NA, nrow = 10, ncol = Nsim)
checkCover_freq <- matrix(NA, nrow = 7, ncol = Nsim)
checkInterval_freq <- matrix(NA, nrow = 7, ncol =Nsim)
checkMSE_freq <- matrix(NA, nrow = 7, ncol = Nsim)

#Begin simulation
for (i in 1:Nsim) { 
start.time <- Sys.time()
# Generate true propensity score for simulation "i"
# Notice, coefficients are sampled uniquely for each "i" and can lead 
# to different results between simulations if a seed is not set prior to the simulation 
true.prop.score  <- burr.link(as.matrix(X.Sim)%*%as.matrix(rnorm(ncol(X.Sim),0,1)))


# sample  treatment  and control adverse events based on true probabilities
samp.treat.hf <- rbinom(length(true.probs.treat.hf),1,true.probs.treat.hf)
samp.control.hf <- rbinom(length(true.probs.control.hf),1,true.probs.control.hf)


# create observed and missing adverse event and death outcomes based on treatment assignment
Wobs <- rbinom(Nobs,1,prob =  true.prop.score) 
TreatAssign <- Wobs==1

obs.treat.hf <- samp.treat.hf[TreatAssign]
obs.control.hf <- samp.control.hf[!TreatAssign]

mis.treat.hf <- samp.treat.hf[!TreatAssign]
mis.control.hf <- samp.control.hf[TreatAssign]

true.probs.treat.death <- burr.link(X.Sim%*%beta2.treat.true + samp.treat.hf*gamma.treat.true)
true.probs.control.death <- burr.link(X.Sim%*%beta2.control.true + samp.control.hf*gamma.control.true)

samp.treat.death <- rbinom(length(true.probs.treat.death),1,true.probs.treat.death)
samp.control.death <- rbinom(length(true.probs.control.death),1,true.probs.control.death)

obs.treat.death <- samp.treat.death[TreatAssign]
obs.control.death <- samp.control.death[!TreatAssign]

mis.treat.death <- samp.treat.death[!TreatAssign]
mis.control.death <- samp.control.death[TreatAssign]

samp.treat.group <- assign.group(samp.treat.hf,samp.treat.death)
samp.control.group <- assign.group(samp.control.hf,samp.control.death)

obs.treat.group <- assign.group(obs.treat.hf,obs.treat.death)
obs.control.group <- assign.group(obs.control.hf,obs.control.death)


# Estimate propensity score based on genreated treatment assignment
est.prop <- Prop_Scores(Wobs,X.Sim[,-1])


# Separate data set into treatment and control subjects
X.treatSim <- X.Sim[TreatAssign,-1]
X.controlSim <- X.Sim[!TreatAssign,-1]


# Create matrix to define the spline on logit of the propensity score and add to data sets
# Algorithm below starts at 6 subclasses removes subclasses until three remain
# to ensure that at least three members of treatment and control are in each subclass
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



b.whole = naturalSpline(logit(est.prop),
                        knots= knots, intercept=F)

X.treatSim <- cbind(b.whole[TreatAssign,], X.treatSim[,-1])
X.controlSim <- cbind(b.whole[!TreatAssign,], X.controlSim[,-1])



# Obtain MLE and Fisher information matrix for sampling from asymptotic posterior distribution
# of the proposed Bayesian model's parameters

X.control.check.hf <- as.data.frame(cbind(obs.control.hf,X.controlSim))
  control.check.hf <- summary(glm(obs.control.hf~.-1, data = X.control.check.hf, family = binomial()))
beta_init.control <- control.check.hf$coefficients[,1]
beta_cov.control <-   control.check.hf$cov.unscaled

X.control.check.death <- as.data.frame(cbind(obs.control.death,X.controlSim, obs.control.hf))
control.check.death <- summary(glm(obs.control.death~.-1, data = X.control.check.death, family = binomial()))
beta2_init.control <- control.check.death$coefficients[,1]
beta2_cov.control <- control.check.death$cov.unscaled

X.treat.check.hf <- as.data.frame(cbind(obs.treat.hf,X.treatSim))
treat.check.hf <- summary(glm(obs.treat.hf~.-1, data = X.treat.check.hf, family = binomial()))
beta_init.treat <- treat.check.hf$coefficients[,1]
beta_cov.treat <-   treat.check.hf$cov.unscaled

X.treat.check.death <- as.data.frame(cbind(obs.treat.death,X.treatSim, obs.treat.hf))
treat.check.death <- summary(glm(obs.treat.death~.-1, data = X.treat.check.death, family = binomial()))
beta2_init.treat <- treat.check.death$coefficients[,1]
beta2_cov.treat <- treat.check.death$cov.unscaled


# Sample model parameters from approximating MVN distributions
beta.control.est<- as.matrix(mvrnorm(500, beta_init.control,beta_cov.control))

beta2.control.gamma <- mvrnorm(500, beta2_init.control,beta2_cov.control)
gamma.control.est <- beta2.control.gamma[,ncol(beta2.control.gamma)]
beta2.control.est <- as.matrix(beta2.control.gamma[,-ncol(beta2.control.gamma)])

beta.treat.est<- as.matrix(mvrnorm(500, beta_init.treat,beta_cov.treat))

beta2.treat.gamma <- mvrnorm(500, beta2_init.treat,beta2_cov.treat)
gamma.treat.est <- beta2.treat.gamma[,ncol(beta2.treat.gamma)]
beta2.treat.est <- as.matrix(beta2.treat.gamma[,-ncol(beta2.treat.gamma)])


# Combine treatment and control data sets for estimation of super-population outcomes
X.controlEst <- X.controlSim
X.treatEst <-  X.treatSim
X.Est <- rbind(X.treatEst, X.controlEst)

#Calculate posterior estimates of probability of outcomes
samp.probs.treat.hf <-  inv.logit(tcrossprod(X.Est,beta.treat.est))
samp.probs.control.hf<-  inv.logit(tcrossprod(X.Est,beta.control.est))

XBeta2.treat <- tcrossprod(X.Est,beta2.treat.est)
samp.probs.treat.yesdeath.yesHF = inv.logit(XBeta2.treat+t(replicate(Nobs, gamma.treat.est)))
samp.probs.treat.yesdeath.noHF = inv.logit(XBeta2.treat)
samp.probs.treat.nodeath.yesHF  = 1- samp.probs.treat.yesdeath.yesHF
samp.probs.treat.nodeath.noHF = 1- samp.probs.treat.yesdeath.noHF

samp.probs.treat.G1 <- samp.probs.treat.nodeath.noHF*(1-samp.probs.treat.hf)
samp.probs.treat.G2 <- samp.probs.treat.nodeath.yesHF*(samp.probs.treat.hf)
samp.probs.treat.G3 <- samp.probs.treat.yesdeath.noHF*(1-samp.probs.treat.hf)
samp.probs.treat.G4 <- samp.probs.treat.yesdeath.yesHF*(samp.probs.treat.hf)

samp.treat.G1 <- colSums(samp.probs.treat.G1)
samp.treat.G2 <- colSums(samp.probs.treat.G2)
samp.treat.G3 <- colSums(samp.probs.treat.G3)
samp.treat.G4 <- colSums(samp.probs.treat.G4)


XBeta2.control <- tcrossprod(X.Est,beta2.control.est)
samp.probs.control.yesdeath.yesHF = inv.logit(XBeta2.control+t(replicate(Nobs, gamma.control.est)))
samp.probs.control.yesdeath.noHF = inv.logit(XBeta2.control)
samp.probs.control.nodeath.yesHF  = 1- samp.probs.control.yesdeath.yesHF
samp.probs.control.nodeath.noHF = 1- samp.probs.control.yesdeath.noHF

samp.probs.control.G1 <- samp.probs.control.nodeath.noHF*(1-samp.probs.control.hf)
samp.probs.control.G2 <- samp.probs.control.nodeath.yesHF*(samp.probs.control.hf)
samp.probs.control.G3 <- samp.probs.control.yesdeath.noHF*(1-samp.probs.control.hf)
samp.probs.control.G4 <- samp.probs.control.yesdeath.yesHF*(samp.probs.control.hf)

samp.control.G1 <- colSums(samp.probs.control.G1)
samp.control.G2 <- colSums(samp.probs.control.G2)
samp.control.G3 <- colSums(samp.probs.control.G3)
samp.control.G4 <- colSums(samp.probs.control.G4)


samp.diff.G1 <- samp.treat.G1/Nobs-samp.control.G1/Nobs
samp.diff.G2 <- samp.treat.G2/Nobs-samp.control.G2/Nobs
samp.diff.G3 <- samp.treat.G3/Nobs-samp.control.G3/Nobs
samp.diff.G4 <- samp.treat.G4/Nobs-samp.control.G4/Nobs


sampExpITTheart <- (samp.treat.G2+samp.treat.G4)/Nobs - (samp.control.G2+samp.control.G4)/Nobs
sampExpITTdeath <- (samp.treat.G3+samp.treat.G4)/Nobs - (samp.control.G3+samp.control.G4)/Nobs
sampExpCompBin <- (samp.treat.G2+samp.treat.G3+samp.treat.G4)/Nobs - (samp.control.G2+samp.control.G3+samp.control.G4)/Nobs
sampExpSACE <- ((samp.treat.G2)/(samp.treat.G1+samp.treat.G2)) - ((samp.control.G2)/(samp.control.G1+samp.control.G2))



Pr_G1_gr_G0 <- samp.probs.control.G1*(samp.probs.treat.G2+samp.probs.treat.G3+samp.probs.treat.G4)+
  samp.probs.control.G2*(samp.probs.treat.G3+samp.probs.treat.G4)+
  samp.probs.control.G3*(samp.probs.treat.G4)

Pr_G0_gr_G1 <- samp.probs.treat.G1*(samp.probs.control.G2+samp.probs.control.G3+samp.probs.control.G4)+
  samp.probs.treat.G2*(samp.probs.control.G3+samp.probs.control.G4)+
  samp.probs.treat.G3*(samp.probs.control.G4)


sampExpDiff_k10k01 <- colMeans(Pr_G1_gr_G0)-colMeans(Pr_G0_gr_G1)
sampExpRR_k10k01 <- colMeans(Pr_G1_gr_G0)/colMeans(Pr_G0_gr_G1)


# Check Coverage of true estimand values by the Bayesian 95% credible interval
checkCover[1,i] <- as.numeric(trueExpITTheart >= quantile(sampExpITTheart, c(0.025,0.975))[1] &
                                trueExpITTheart <= quantile(sampExpITTheart, c(0.025,0.975))[2] )

checkCover[2,i] <- as.numeric(trueExpITTdeath >= quantile(sampExpITTdeath, c(0.025,0.975))[1] &
                                trueExpITTdeath <= quantile(sampExpITTdeath, c(0.025,0.975))[2] )

checkCover[3,i] <- as.numeric(trueExpCompBin >= quantile(sampExpCompBin , c(0.025,0.975))[1] &
                                trueExpCompBin <= quantile(sampExpCompBin, c(0.025,0.975))[2] )

checkCover[4,i] <- as.numeric(trueExpSACE >= quantile(sampExpSACE , c(0.025,0.975))[1] &
                                trueExpSACE <= quantile(sampExpSACE , c(0.025,0.975))[2] )

checkCover[5,i] <- as.numeric(expected.diff.G1 >= quantile(samp.diff.G1 , c(0.025,0.975))[1] &
                                expected.diff.G1 <= quantile(samp.diff.G1 , c(0.025,0.975))[2] )

checkCover[6,i] <- as.numeric(expected.diff.G2 >= quantile(samp.diff.G2 , c(0.025,0.975))[1] &
                                expected.diff.G2 <= quantile(samp.diff.G2 , c(0.025,0.975))[2] )

checkCover[7,i] <- as.numeric(expected.diff.G3 >= quantile(samp.diff.G3 , c(0.025,0.975))[1] &
                                expected.diff.G3 <= quantile(samp.diff.G3 , c(0.025,0.975))[2] )

checkCover[8,i] <- as.numeric(expected.diff.G4 >= quantile(samp.diff.G4 , c(0.025,0.975))[1] &
                                expected.diff.G4 <= quantile(samp.diff.G4 , c(0.025,0.975))[2] )

checkCover[9,i] <- as.numeric(trueExpDiff_k10k01 >= quantile(sampExpDiff_k10k01 , c(0.025,0.975))[1] &
                                trueExpDiff_k10k01 <= quantile(sampExpDiff_k10k01 , c(0.025,0.975))[2] )

checkCover[10,i] <- as.numeric(trueExpRR_k10k01 >= quantile(sampExpRR_k10k01 , c(0.025,0.975))[1] &
                                 trueExpRR_k10k01<= quantile(sampExpRR_k10k01 , c(0.025,0.975))[2] )


# Calculate interval length of Bayesian 95% credible interval
checkInterval[1, i] <- quantile(sampExpITTheart, c(0.025,0.975))[2] - quantile(sampExpITTheart, c(0.025,0.975))[1]

checkInterval[2, i] <- quantile(sampExpITTdeath, c(0.025,0.975))[2] - quantile(sampExpITTdeath, c(0.025,0.975))[1]

checkInterval[3, i] <- quantile(sampExpCompBin, c(0.025,0.975))[2] - quantile(sampExpCompBin, c(0.025,0.975))[1]

checkInterval[4, i] <- quantile(sampExpSACE , c(0.025,0.975))[2] - quantile(sampExpSACE , c(0.025,0.975))[1]

checkInterval[5, i] <- quantile(samp.diff.G1  , c(0.025,0.975))[2] - quantile(samp.diff.G1  , c(0.025,0.975))[1]

checkInterval[6, i] <- quantile(samp.diff.G2  , c(0.025,0.975))[2] - quantile(samp.diff.G2  , c(0.025,0.975))[1]

checkInterval[7, i] <- quantile(samp.diff.G3  , c(0.025,0.975))[2] - quantile(samp.diff.G3  , c(0.025,0.975))[1]

checkInterval[8, i] <- quantile(samp.diff.G4  , c(0.025,0.975))[2] - quantile(samp.diff.G4  , c(0.025,0.975))[1]

checkInterval[9, i] <- quantile(sampExpDiff_k10k01 , c(0.025,0.975))[2] - quantile(sampExpDiff_k10k01 , c(0.025,0.975))[1]

checkInterval[10, i] <- quantile(sampExpRR_k10k01  , c(0.025,0.975))[2] - quantile(sampExpRR_k10k01  , c(0.025,0.975))[1]


# Calculate bias of estimate from Bayesian model, later used to compute MSE as well
checkMSE[1, i] <- mean(sampExpITTheart) - trueExpITTheart

checkMSE[2, i] <- mean(sampExpITTdeath) - trueExpITTdeath

checkMSE[3, i] <- mean(sampExpCompBin) - trueExpCompBin

checkMSE[4, i] <- mean(sampExpSACE) - trueExpSACE

checkMSE[5, i] <- mean(samp.diff.G1) - expected.diff.G1

checkMSE[6, i] <- mean(samp.diff.G2) - expected.diff.G2

checkMSE[7, i] <- mean(samp.diff.G3) - expected.diff.G3

checkMSE[8, i] <- mean(samp.diff.G4) - expected.diff.G4

checkMSE[9, i] <- mean(sampExpDiff_k10k01) - trueExpDiff_k10k01

checkMSE[10, i] <- mean(sampExpRR_k10k01) - trueExpRR_k10k01


# Compute DR estimates of traditional and composite ordinal estimands
obs.hf <- c(obs.treat.hf,obs.control.hf)
obs.death <- c(obs.treat.death,obs.control.death)
obs.CompBin  <- as.numeric(obs.death+obs.hf>0)
obs.DiffG1 <- as.numeric(obs.death+obs.hf==0)
obs.DiffG2 <- as.numeric(obs.death-obs.hf==-1)
obs.DiffG3 <- as.numeric(obs.death-obs.hf==1)
obs.DiffG4 <- as.numeric(obs.death+obs.hf==2)
ordered.treat <- c(rep(1,sum(Wobs)),rep(0,sum(1-Wobs)))
ordered.pscore <- c(est.prop[TreatAssign],est.prop[!TreatAssign])
ordered.X <- as.data.frame(X.Sim[c(which(TreatAssign),which(!TreatAssign)),-1])

freqITTheart <- DR_Est(obs.hf)
freqITTdeath <- DR_Est(obs.death)
freqCompBin <- DR_Est(obs.CompBin)
freqDiffG1 <- DR_Est(obs.DiffG1)
freqDiffG2 <- DR_Est(obs.DiffG2)
freqDiffG3 <- DR_Est(obs.DiffG3)
freqDiffG4 <- DR_Est(obs.DiffG4)

# Check coverage of true parameter by  DR confidence interval
checkCover_freq[1,i] <- as.numeric(trueExpITTheart > freqITTheart[1] & trueExpITTheart < freqITTheart[3])
checkCover_freq[2,i] <- as.numeric(trueExpITTdeath > freqITTdeath[1] & trueExpITTdeath < freqITTdeath[3])
checkCover_freq[3,i] <- as.numeric(trueExpCompBin > freqCompBin[1] & trueExpCompBin < freqCompBin[3])
checkCover_freq[4,i] <- as.numeric(expected.diff.G1 > freqDiffG1[1] & expected.diff.G1 < freqDiffG1[3])
checkCover_freq[5,i] <- as.numeric(expected.diff.G2 > freqDiffG2[1] & expected.diff.G2 < freqDiffG2[3])
checkCover_freq[6,i] <- as.numeric(expected.diff.G3 > freqDiffG3[1] & expected.diff.G3 < freqDiffG3[3])
checkCover_freq[7,i] <- as.numeric(expected.diff.G4 > freqDiffG4[1] & expected.diff.G4 < freqDiffG4[3])

# Calculate interval length of DR confidence interval
checkInterval_freq[1,i] <- freqITTheart[3]-freqITTheart[1]
checkInterval_freq[2,i] <- freqITTdeath[3]-freqITTdeath[1]
checkInterval_freq[3,i] <- freqCompBin[3]-freqCompBin[1]
checkInterval_freq[4,i] <- freqDiffG1[3]-freqDiffG1[1]
checkInterval_freq[5,i] <- freqDiffG2[3]-freqDiffG2[1]
checkInterval_freq[6,i] <- freqDiffG3[3]-freqDiffG3[1]
checkInterval_freq[7,i] <- freqDiffG4[3]-freqDiffG4[1]

# Calculate bias of DR estimate later to be used to calculate MSE
checkMSE_freq[1,i] <- freqITTheart[2] - trueExpITTheart
checkMSE_freq[2,i] <- freqITTdeath[2] - trueExpITTdeath
checkMSE_freq[3,i] <- freqCompBin[2]- trueExpCompBin
checkMSE_freq[4,i] <- freqDiffG1[2] - expected.diff.G1
checkMSE_freq[5,i] <- freqDiffG2[2] - expected.diff.G2
checkMSE_freq[6,i] <- freqDiffG3[2] - expected.diff.G3
checkMSE_freq[7,i] <- freqDiffG4[2] - expected.diff.G4

end.time<-Sys.time()
#Print metrics to track simulation results and progress
print(i)
if(i>1){
  print("Bayesian Method Coverage:")
  print(rowMeans(checkCover[,1:i]))
  print("DR Coverage:")
  print(rowMeans(checkCover_freq[,1:i]))}
}

# Coverage, Interval Width, Bias, RMSE for each estimand under Bayesian Model
rowMeans(checkCover)
rowMeans(checkInterval)
rowMeans(checkMSE)
sqrt(rowMeans((checkMSE)^2))


# Coverage, Interval Width, Bias, RMSE for each estimand under IPW
rowMeans(checkCover_freq)
rowMeans(checkInterval_freq)
rowMeans(checkMSE_freq)
sqrt(rowMeans((checkMSE_freq)^2))



