
#----------------------------------- CIM ################################

library(tidyverse)
library(parallel)

 

########################## Project 1 ####################################

##------------------------------- Question 1 ---------------------##
library("DAAG") 
data(nassCDS) 
names(nassCDS) 

table(nassCDS$dead)
summary(nassCDS$ageOFocc)
nassCDS$Dead<- ifelse(nassCDS$dead=="dead", 1,0)

### 1. Estimate the model using a classical GLM approach. 
# Fit the GLM
model <- glm(Dead ~ ageOFocc, family = binomial, data = nassCDS)
# Summary of the model
summary(model)

beta_0= coef(model)[1]
beta_1= coef(model)[2]

###  2.  Estimate 𝑋50. Use non parametric bootstrap to estimate the 
###   distribution of 𝑋50 and to construct a 95% C.I. for the 𝑋50

#pi(ageOFocc)= exp(beta_0 = beta_1*ageOFocc)/ (1+ exp(beta_0 = beta_1*ageOFocc))
X_50= - beta_0/beta_1

# non parametric bootstrap

#-------------------------------------------------------------------------
# Setup the number of cores 
num_cores <- detectCores() - 2 

# Define the parallel bootstrap function
bootstrap_parallel <- function(i, n, nassCDS, index) {
  index.b <- sample(index, n, replace = TRUE)
  nassCDS.b <- nassCDS[index.b, ]
  fit.glm.b <- glm(nassCDS.b$Dead ~ nassCDS.b$ageOFocc, family = binomial)
  beta_0.b <- summary(fit.glm.b)$coeff[1, 1]
  beta_1.b <- summary(fit.glm.b)$coeff[2, 1]
  X_50.b <- -beta_0.b / beta_1.b
  return(X_50.b)
}

# Setup parallel computation
B <- 100  
n <- length(nassCDS$Dead)  
index <- 1:n  

# Use parLapply for parallel execution
X_50.b<- parLapply(cl = makeCluster(num_cores), X = 1:B, 
                             fun = bootstrap_parallel, n = n, nassCDS = nassCDS, index = index)

# Stop the cluster once done
stopCluster(cl)
X_50.b<-unlist(X_50.b)
mean(X_50.b)
hist(X_50.b, nclass=50)
#-----------------------------------------------------------------------

set.seed(123)
n<-length(nassCDS$Dead)
B<-100
beta_0.b<-beta_1.b<- X_50.b<-c(1:B)
index<-c(1:n)

for (i in 1:B)
{
  index.b<-sample(index,n,replace=TRUE)
  nassCDS.b<-nassCDS[index.b,]
  fit.glm.b<-glm(nassCDS.b$Dead~nassCDS.b$ageOFocc,family = binomial)
  beta_0.b<-summary(fit.glm.b)$coeff[1,1]
  beta_1.b<-summary(fit.glm.b)$coeff[2,1]
  X_50.b[i]<- -beta_0.b/beta_1.b
  
}

hist(X_50.b,nclass=50)
quantile(X_50.b,probs=c(0.025,0.975))
mean(X_50.b)

### 3. Estimate the OR (for a unit increased in age)
OR= exp(beta_1)

# non parametric bootstrap to construct a 95% C.I. for the OR : percentile 
set.seed(123)
n<-length(nassCDS$Dead)
B<-100
beta_1.b<- OR.b<-c(1:B)
index<-c(1:n)

for (i in 1:B)
{
  index.b<-sample(index,n,replace=TRUE)
  nassCDS.b<-nassCDS[index.b,]
  fit.glm.b<-glm(nassCDS.b$Dead~nassCDS.b$ageOFocc,family = binomial)
  beta_1.b<-summary(fit.glm.b)$coeff[2,1]
  OR.b[i]<- exp(beta_1.b)
  
}

mean(OR.b)
hist(OR.b,nclass=50)
quantile(OR.b,probs=c(0.025,0.975))

# non parametric bootstrap to construct a 95% C.I. for the OR :  bootstrap t interval methods 
set.seed(123)
n<-length(nassCDS$Dead)
B<-100
beta_1.b<- t.OR<- se.OR<- c(1:B)
index<-c(1:n)
OR= exp(beta_1)

for(b in 1:B)
{
  index.b<-sample(index,n,replace=TRUE)
  nassCDS.b<-nassCDS[index.b,]
  fit.glm.b<-glm(nassCDS.b$Dead~nassCDS.b$ageOFocc,family = binomial)
  beta_1.b<-summary(fit.glm.b)$coeff[2,1]
  OR.b[i]<- exp(beta_1.b)
  se.OR[i]<-sqrt(var(OR.b[i])/n)
  t.OR[i]<-(mean(OR.b[i])-OR)/se.OR[i]
}
quantile(na.omit(t.OR)  ,probs=c(0.05,0.95))

hist(t.OR,probability=T,nclass=100,ylim=c(0,0.5),xlim=c(-10,10))
xx<-seq(from=-3,to=3,length=1000)
dx2<-dt(xx,6)
lines(xx,dx2,col=2)
lines(rep(qt(0.95,6),2),c(0,0.4),col=2)
calpha<-quantile(t.OR,probs=c(0.95))
lines(rep(calpha,2),c(0,0.4),col=4)
qt(0.95,6)
quantile(na.omit(t.OR),probs=c(0.95))


### 4. Use parametric bootstrap  to test the null hypothesis 𝐻0: 𝑂𝑅 = 1

prob.i<-exp(beta_0)/(1+ exp(beta_0))
prob.i

B<-100
beta0.b<-beta1.b<- OR.b<-c(1:B)
Dead.b<- c(1:n)
prob<- rep(prob.i, n)

for (i in 1:B)
{
  for(j in 1:n)
  {
    Dead.b[j]<-rbinom(1,1,prob[j])
  }
  fit.glm.b<-glm(Dead.b~nassCDS$ageOFocc,family = binomial)
  beta0.b[i]<-summary(fit.glm.b)$coeff[1,1]
  beta1.b[i]<-summary(fit.glm.b)$coeff[2,1]
  OR.b[i]<- exp(summary(fit.glm.b)$coeff[2,1])
}

hist(OR.b,nclass=50)
quantile(OR.b,probs=c(0.025,0.975))

Pvalue<-(1+sum(OR.b >1)/(B+1))
Pvalue


### 5. parametric bootstrap to calculate the standard error for 𝜋̂33 and construct a 90% C.I. for 𝜋33
B<-100
beta0.b<-beta1.b<-pi33.b<- se_pi33.b <- c(1:B)
Dead.b<-c(1:n)
for (i in 1:B)
{
  for(j in 1:n)
  {
    Dead.b[j]<-rbinom(1,1,exp(beta_0+beta_1*nassCDS$ageOFocc[j])/(1+ exp(beta_0+beta_1*nassCDS$ageOFocc[j])))
  }
  fit.glm.b<-glm(Dead.b~nassCDS$ageOFocc,family = binomial)
  beta0.b[i]<-summary(fit.glm.b)$coeff[1,1]
  beta1.b[i]<-summary(fit.glm.b)$coeff[2,1]
  pi33.b[i]<- exp(beta0.b[i]+beta1.b[i]*33)/(1+ exp(beta0.b[i]+beta1.b[i]*33))
  se_pi33.b[i]<-sqrt(var(pi33.b)/n) 
}

hist(pi33.b,nclass=50)
mean(pi33.b)
quantile(pi33.b,probs=c(0.05,0.950))

hist(se_pi33.b,nclass=50)
mean(se_pi33.b)
quantile(se_pi33.b,probs=c(0.05,0.950))




##------------------------------- Question 2 ---------------------##

### 1. Estimate the model using the R package glmRob
library(robust)
model.rob<- glmRob(Dead ~ ageOFocc, family = binomial(),
                   data = nassCDS, method = "cubif")
# Summary of the model
summary(model.rob)

betarob_0= coef(model.rob)[1]
betarob_1= coef(model.rob)[2]


### 2. Use non parametric bootstrap to estimate the SE for the intercept and slope
set.seed(123)
x<- nassCDS$ageOFocc
y<- nassCDS$Dead
n<-length(nassCDS$Dead)
B<-30
beta0<-beta1<-se.beta0<-se.beta1<-c(1:B)
index<-c(1:n)

for(i in 1:B)
{
  index.b<-sample(x,n,replace=TRUE)
  ageOFocc.b<-x[index.b]
  Dead.b<-y[index.b]
  fit.rob.b<- glmRob(Dead.b ~ ageOFocc.b, family = binomial(), method = "cubif")
  beta0[i]<-summary(fit.rob.b)$coeff[1,1]
  beta1[i]<-summary(fit.rob.b)$coeff[2,1]

  #se.beta0[i]<- sqrt(var(beta0)/n) 
  #se.beta1[i]<- sqrt(var(beta1)/n) 
}

hist(beta0,nclass=50)
mean(beta0)
se.beta0<-  sqrt(var(beta0)/length(beta0))


hist(beta1,nclass=50)
mean(beta1)
se.beta1<- sqrt(var(beta1)/length(beta1))

### 3.Jackknife and the bootstrap procedures to estimate the bias and MSE for the intercept and slope 

#------------------------------------- GLM 1
## Bias

#bootstrap
beta0r.obs<-betarob_0
beta1r.obs<- betarob_1
n<-length(nassCDS$Dead)
B<-30
index<-c(1:n)
x<- nassCDS$ageOFocc
y<- nassCDS$Dead
n<-length(nassCDS$Dead)
beta0r<-beta1r<-mse_beta0.b <- mse_beta1.b<-  c(1:B)

for(i in 1:B)
{
  index.b<-sample(x,n,replace=TRUE)
  ageOFocc.b<-x[index.b]
  Dead.b<-y[index.b]
  # Fit the model
  fit.rob.b<- glm(Dead.b ~ ageOFocc.b, family = binomial)
  #Extract coefficient
  beta0r[i]<-summary(fit.rob.b)$coeff[1,1]
  beta1r[i]<-summary(fit.rob.b)$coeff[2,1]
  #MSE for Parameters
  mse_beta0.b[i] <- (beta0r[i] - beta_0)^2
  mse_beta1.b[i] <- (beta1r[i] - beta_1)^2
}
   # bias
hist(beta0r,nclass=50)
beta0r.b <- mean(beta0r)
bias_beta0.r <- beta0r.b - beta0r.obs
bias_beta0.r

hist(beta1r,nclass=50)
beta1r.b <- mean(beta1r)
bias_beta1.r <- beta1r.b - beta1r.obs
bias_beta1

   # MSE
hist(mse_beta0.b,nclass=50)
mse_beta0.boot <- mean( mse_beta0.b)

hist(mse_beta1.b,nclass=50)
mse_beta1.boot <- mean(mse_beta1.b)


#Jackknife
beta0r.obs<-betarob_0
beta1r.obs<- betarob_1
n<-length(nassCDS$Dead)
index<-c(1:n)
x<- nassCDS$ageOFocc
y<- nassCDS$Dead
n<-length(nassCDS$Dead)
beta0rj<-beta1rj<- mse_beta0.j <- mse_beta1.j <-c(1:n)

for(i in 1:n)
{

  ageOFocc.j<-x[- c(i)]
  Dead.j<-y[- c(i)]
  # Fit the model
  fit.rob.j<- glm(Dead.j ~ ageOFocc.j, family = binomial)
  #Extract coefficient
  beta0rj[i]<-summary(fit.rob.j)$coeff[1,1]
  beta1rj[i]<-summary(fit.rob.j)$coeff[2,1]
  #MSE for Parameters
  mse_beta0.j[i] <- (beta0rj[i] - beta_0)^2
  mse_beta1.b[i] <- (beta1rj[i] - beta_1)^2
}

  # Bias
hist(beta0rj,nclass=50)
beta0r.j <- mean(beta0rj)
bias_beta0.rj <- beta0r.j - beta0r.obs
bias_beta0.rj

hist(beta1rj,nclass=50)
beta1r.j <- mean(beta1rj)
bias_beta1.rj <- beta1r.j - beta1r.obs
bias_beta1.rj

  # MSE
# MSE
hist(mse_beta0.j,nclass=50)
mse_beta0.jack <- mean( mse_beta0.j)

hist(mse_beta1.j,nclass=50)
mse_beta1.jack <- mean(mse_beta1.j)
#------------------------------------------


#------------------------------------- GLM 2: robust GLM
## Bias

#bootstrap
beta0r.obs<-betarob_0
beta1r.obs<- betarob_1
n<-length(nassCDS$Dead)
B<-100
index<-c(1:n)
x<- nassCDS$ageOFocc
y<- nassCDS$Dead
n<-length(nassCDS$Dead)
beta0r<-beta1r<- mse_beta0.b <- mse_beta1.b<- c(1:B)

for(i in 1:B)
{
  index.b<-sample(x,n,replace=TRUE)
  ageOFocc.b<-x[index.b]
  Dead.b<-y[index.b]
  # Fit the model
  fit.rob.b<- glmRob(Dead.b ~ ageOFocc.b, family = binomial(), method = "cubif")
  #Extract coefficient
  beta0r[i]<-summary(fit.rob.b)$coeff[1,1]
  beta1r[i]<-summary(fit.rob.b)$coeff[2,1]
  #MSE for Parameters
  mse_beta0.b[i] <- (beta0r[i] - beta_0)^2
  mse_beta1.b[i] <- (beta1r[i] - beta_1)^2
}
# bias
hist(beta0r,nclass=50)
beta0r.b <- mean(beta0r)
bias_beta0.r <- beta0r.b - beta0r.obs
bias_beta0.r

hist(beta1r,nclass=50)
beta1r.b <- mean(beta1r)
bias_beta1.r <- beta1r.b - beta1r.obs
bias_beta1

# MSE
hist(mse_beta0.b,nclass=50)
mse_beta0.boot <- mean( mse_beta0.b)

hist(mse_beta1.b,nclass=50)
mse_beta1.boot <- mean(mse_beta1.b)


#Jackknife
beta0r.obs<-betarob_0
beta1r.obs<- betarob_1
n<-length(nassCDS$Dead)
index<-c(1:n)
x<- nassCDS$ageOFocc
y<- nassCDS$Dead
n<-length(nassCDS$Dead)
beta0rj<-beta1rj<- mse_beta0.j <- mse_beta1.j <- c(1:n)

for(i in 1:n)
{
  cat(i)
  ageOFocc.j<-x[- c(i)]
  Dead.j<-y[- c(i)]
  # Fit the model
  fit.rob.j<- glmRob(Dead.j ~ ageOFocc.j, family = binomial(), method = "cubif")
  #Extract coefficient
  beta0rj[i]<-summary(fit.rob.j)$coeff[1,1]
  beta1rj[i]<-summary(fit.rob.j)$coeff[2,1]
  #MSE for Parameters
  mse_beta0.j[i] <- (beta0rj[i] - beta_0)^2
  mse_beta1.b[i] <- (beta1rj[i] - beta_1)^2
}

# Bias
hist(beta0rj,nclass=50)
beta0r.j <- mean(beta0rj)
bias_beta0.rj <- beta0r.j - beta0r.obs
bias_beta0.rj

hist(beta1rj,nclass=50)
beta1r.j <- mean(beta1rj)
bias_beta1.rj <- beta1r.j - beta1r.obs
bias_beta1.rj

# MSE
# MSE
hist(mse_beta0.j,nclass=50)
mse_beta0.jack <- mean( mse_beta0.j)

hist(mse_beta1.j,nclass=50)
mse_beta1.jack <- mean(mse_beta1.j)
#------------------------------------------

################################ Question 3

### 1. Define the observation unit (𝑋𝑖,𝑌𝑖) for the question
XY<- table(nassCDS$airbag, nassCDS$dead)

### 2. Odd ratio
library(epitools)

# Create the 2x2 table
XY <- matrix(c(11058, 13825, 669, 511), 
                nrow = 2, 
                byrow = FALSE,
                dimnames = list(airbag = c("None", "Airbag"), 
                                dead = c("Alive", "Dead")))

# Calculate the odds ratio and 95% CI
OR <- oddsratio(as.table(XY), method = "wald")

# Interpretation
   # if the OR is greater than 1, airbags are associated with increased odds of survival.
   # If the OR is less than 1, airbags are associated with decreased odds of survival.

### 3. Use parametric bootstrap to construct a construct a 95%confidence interval for the  OR.
  #Simulate B samples from the multinomial distribution using the observed cell proportions. For each bootstrap sample:
  #Reconstruct the contingency table.
  #Compute the odds ratio

# Convert counts to proportions
total <- sum(XY)
proportions <- XY / total

# Function to calculate the odds ratio
calc_or <- function(table) {
  alive_none <- table[1, 1]
  dead_none <- table[1, 2]
  alive_airbag <- table[2, 1]
  dead_airbag <- table[2, 2]
  (alive_airbag / dead_airbag) / (alive_none / dead_none)
}

# Parametric bootstrap
set.seed(123)
B <- 1000
ORs <- numeric(B)

for (i in 1:B) {
  # Generate bootstrap sample
  bootstrap_sample <- rmultinom(1, size = total, prob = as.vector(proportions))
  bootstrap_table <- matrix(bootstrap_sample, nrow = 2, byrow = TRUE)
  
  # Calculate odds ratio for the bootstrap sample
  mean(ORs)
  ORs[i] <- calc_or(bootstrap_table)
}

# 95% Confidence Interval
hist(ORs,nclass=50)
OR<- mean(ORs)
ci <- quantile(ORs, probs = c(0.025, 0.975))

### 4.Use permutations test to test the hypothesis that airbags in the car DO NOT influence 
    # the accident outcome (testing OR=1)

# Observed chi-square test statistic
observed_test <- chisq.test(XY, correct = FALSE)
observed_stat <- observed_test$statistic

# Permutation test
set.seed(123) 
B <- 10000
permuted_stats <- numeric(B)

# Marginal totals
row_totals <- rowSums(XY)
col_totals <- colSums(XY)
total <- sum(XY)

for (i in 1:B) {
  # Generate a permuted table under the null hypothesis
  permuted_table <- matrix(rmultinom(1, total, prob = outer(row_totals, col_totals, "*") / total^2),
                           nrow = 2)
  
  # Compute chi-square test statistic for permuted table
  permuted_test <- chisq.test(permuted_table, correct = FALSE)
  permuted_stats[i] <- permuted_test$statistic
}

# Calculate p-value
hist(permuted_stats,nclass=50,probability=TRUE)
p_value<- (1+sum(permuted_stats>observed_stat))/(B+1)


### Comparison of  chi-square test statistic 

# Plot the distributions
hist(permuted_stats, breaks = 30, probability = TRUE, 
     col = "skyblue", main = "Comparison of Chi-Square Distributions",
     xlab = "Chi-Square Test Statistic")
curve(dchisq(x, df = 1), col = "red", lwd = 2, add = TRUE)

legend("topright", legend = c("Permutation Test (Empirical)", "Chi-Square (Theoretical)"),
       col = c("skyblue", "red"), lwd = 2)



############################## Question 4 : Dead and sex

### 1. Estimate the proportion of male (𝜋𝑀) and female (𝜋𝐹) that died in the accidents

table(nassCDS$dead,nassCDS$sex )

# Data
f_alive <- 11784
f_dead <- 464
m_alive <- 13253
m_dead <- 716

# Total counts
f_total <- f_alive + f_dead
m_total <- m_alive + m_dead

# Proportions
pi_F <- f_dead / f_total
pi_M <- m_dead / m_total

### 2.Test the hypothesis that the proportion of male and female that died in an accident are equal
# Data
f_dead <- 464
f_total <- 11784 + 464
m_dead <- 716
m_total <- 13253 + 716

# Proportions
pi_F <- f_dead / f_total
pi_M <- m_dead / m_total

# Pooled proportion
pooled_pi <- (f_dead + m_dead) / (f_total + m_total)

# Test statistic
z_stat <- (pi_M - pi_F) / sqrt(pooled_pi * (1 - pooled_pi) * (1 / f_total + 1 / m_total))

# p-value for two-sided test
p_value <- 2 * (1 - pnorm(abs(z_stat)))


# Decision
if (p_value < 0.05) {
  cat("Reject the null hypothesis: The proportions are significantly different.\n")
} else {
  cat("Fail to reject the null hypothesis: No significant difference in proportions.\n")
}


## NB: using R function prop.test
# Perform the two-proportion z-test using prop.test
Test_prop <- prop.test(c(f_dead, m_dead), c(f_total, m_total), alternative = "two.sided")

### 3. Use parametric bootstrap to test the hypothesis that the proportion of male and 
    #female that died in an accidents are equal against a two sided alternative. 

set.seed(123) 
# Bootstrap settings
B <- 10000  
boot_stat <- numeric(B)

# Bootstrap process
for (i in 1:B) {
  # Simulate deaths based on the pooled proportion
  f_boot <- rbinom(1, f_total, pooled_pi)
  m_boot <- rbinom(1, m_total, pooled_pi)
  
  # Compute the difference in proportions for the bootstrap sample
  boot_stat[i] <- (f_boot / f_total) - (m_boot / m_total)
}

# Observed difference in proportions
obs_diff <- pi_F - pi_M

# Compute p-value by comparing observed difference to bootstrap distribution
p_value <-(1+sum(abs(boot_stat) > abs(obs_diff)))/(B+1)

# Conclusion based on p-value
if (p_value < 0.05) {
  cat("Reject the null hypothesis: The proportions are significantly different.\n")
} else {
  cat("Fail to reject the null hypothesis: No significant difference in proportions.\n")
}


### 4. Use non-parametric bootstrap construct a 95% confident interval for 𝜋𝑀 − 𝜋�
set.seed(123)

# Observed difference in proportions
obs_diff <- pi_M - pi_F

# Bootstrap settings
B <- 10000  
boot_diff <- numeric(B)

# Non-parametric bootstrap process
for (i in 1:B) {
  # Resample males and females with replacement
  f_boot <- sample(c(rep(1, f_dead), rep(0, f_total - f_dead)), f_total, replace = TRUE)
  m_boot <- sample(c(rep(1, m_dead), rep(0, m_total - m_dead)), m_total, replace = TRUE)
  
  # Compute the difference in proportions for the bootstrap sample
  boot_diff[i] <- mean(f_boot) - mean(m_boot)
}

# Calculate the 95% confidence interval from the bootstrap distribution
hist(boot_diff,nclass=100)
quantile(boot_diff,probs=c(0.025,0.975))
ci_lower <- quantile(boot_diff, 0.025)
ci_upper <- quantile(boot_diff, 0.975)



#########################################################################
library(tidyverse)
library(magrittr)

#-------------------

#---------- Project 2
#===============================================================================
#---------------- Question 1
library(glm2)
data(crabs)
names(crabs)
plot(crabs$Width, crabs$Satellites)
head(crabs) 



# 1.
crab_mod_1 <- glm(Satellites ~ Width, family = poisson(), data = crabs)
crab_mod_2 <- glm(Satellites ~ Width + GoodSpine, family = poisson(), data = crabs)
# likelihood ratio test
lrt_obs <- -2*(logLik(crab_mod_1)[1] - logLik(crab_mod_2)[1])
# P=value
1 - pchisq(lrt_obs, df = 1)
# or
lr_test <- anova(crab_mod_1, crab_mod_2, test = "Chisq")
print(lr_test)



# 2. 
lambda <- mean(crabs$Satellites)
n <- nrow(crabs)
B <- 1000

# Parametric
lrt_para <- c()
set.seed(2023)
for (i in 1:B) {
  Sat_boot <- rpois(n, lambda)
  crab_boot_1 <- glm(Sat_boot ~ crabs$Width, family = poisson())
  crab_boot_2 <- glm(Sat_boot ~ crabs$Width + crabs$GoodSpine, family = poisson())
  lrt_para[i] <- -2*(logLik(crab_boot_1)[1] - logLik(crab_boot_2)[1])
}

# Non-parametric
lrt_nonpara <- c()
set.seed(2023)
for (i in 1:B) {
  Sat_boot <- sample(crabs$Satellites, n, replace = T)
  crab_boot_1 <- glm(Sat_boot ~ crabs$Width, family = poisson())
  crab_boot_2 <- glm(Sat_boot ~ crabs$Width + crabs$GoodSpine, family = poisson())
  lrt_nonpara[i] <- -2*(logLik(crab_boot_1)[1] - logLik(crab_boot_2)[1])
}
par(mfrow = c(1, 2))
hist(lrt_para, main = "Parametric")
hist(lrt_nonpara, main = "Non-parametric")
par(mfrow = c(1, 1))

# P-values
(1 + sum(lrt_para > lrt_obs))/(1 + B) 
(1 + sum(lrt_nonpara > lrt_obs))/(1 + B)

# Plot histograms of the bootstrap distributions
hist(lrt_para, breaks = 50, probability = TRUE, col = "skyblue", main = "Parametric vs Non-Parametric Bootstrap",
     xlab = "Likelihood Ratio Statistic (D)")
hist(lrt_nonpara, breaks = 50, probability = TRUE, col = rgb(1, 0, 0, 0.5), add = TRUE)

# Overlay the theoretical chi-squared distribution
curve(dchisq(x, df = 1), col = "blue", lwd = 2, add = TRUE)

legend("topright", legend = c("Parametric", "Non-Parametric", "Chi-squared"),
       col = c("skyblue", "red", "blue"), lwd = 2)


# 3
# permutaions
lrt_permu <- c()
set.seed(2023)
for (i in 1:B) {
  Sat_boot <- sample(crabs$Satellites, n, replace = F)
  crab_boot_1 <- glm(Sat_boot ~ crabs$Width, family = poisson())
  crab_boot_2 <- glm(Sat_boot ~ crabs$Width + crabs$GoodSpine, family = poisson())
  lrt_permu[i] <- -2*(logLik(crab_boot_1)[1] - logLik(crab_boot_2)[1])
}

hist(lrt_permu)
(1 + sum(lrt_permu > lrt_obs))/(1 + B) # p-value


#------------------- Question 2

extra<-c(0.7,-1.6,-0.2,-1.2,-0.1,3.4,3.7,0.8,0.0,2.0, 
         1.9,0.8,1.1,0.1,-0.1,4.4,5.5,1.6,4.6,4.3) 
group<-c(rep(1,10),rep(2,10)) 
ID<-c(1:20) 
sleep <- data.frame(extra, group, ID)

# 1
t.test(sleep$extra, sleep$group)

# 2
     # -------
t_obs <- t.test(sleep$extra, sleep$group)$statistic
B <- 1000
t_boot <- c()
set.seed(2023)
for (i in 1:B) {
  extra_boot <- sample(sleep$extra, length(sleep$extra), replace = T)
  t_boot[i] <- t.test(extra_boot, sleep$group)$statistic
}
hist(t_boot)

# 2-sided p-value
(1 + sum(abs(t_boot) > abs(t_obs)))/(1 + B)

    # ------ 
t.obs <- t.test(sleep$extra, sleep$group)$statistic
B<-1000
t.boot<-c(1:B)
nn<- nrow(sleep[sleep$group == 1, ])
mm<- nrow(sleep[sleep$group == 2, ])
nm<- nrow(sleep)
z<- sleep[,-3]
for(b in 1:B)
{
  z.b<-sample(sleep$extra, nm, replace = T)
  g1.b<-z.b[1:nn]
  g2.b<-z.b[(nn+1):nm]
  t.boot[b]<-t.test(g1.b,g2.b,,alternative="greater",var.equal=TRUE)$statistic
}
hist(t.boot,nclass=50,probability=T)
lines(c(t.obs,t.obs),c(0,1),lwd=3,col=2)
Pmc<-(1+sum(t.boot > t.obs))/(B+1)
Pmc


# 3
psych::describeBy(sleep$extra, group = sleep$group)

B <- 1000
t_boot_m <- c()
mean_gr1 <- mean(sleep$extra[sleep$group == 1])
sd_gr1 <- sd(sleep$extra[sleep$group == 1])
mean_gr2 <- mean(sleep$extra[sleep$group == 2])
sd_gr2 <- sd(sleep$extra[sleep$group == 2])

set.seed(2023)
for (i in 1:B) {
  extra_boot1 <- rnorm(10, mean_gr1, sd_gr1^2)
  extra_boot2 <- rnorm(10, mean_gr2, sd_gr2^2)
  med1 <- median(extra_boot1)
  med2 <- median(extra_boot2)
  SM <- sum(abs(extra_boot1 - med1)) + sum(abs(extra_boot2 - med2))
  t_boot_m[i] <- (med1 - med2)/SM
}
med_gr1 <- median(sleep$extra[sleep$group == 1])
med_gr2 <- median(sleep$extra[sleep$group == 2])
SM_obs <- sum(abs(sleep$extra[sleep$group == 1] - med_gr1)) + 
  sum(abs(sleep$extra[sleep$group == 2] - med_gr2))
tm_obs <- (med_gr1 - med_gr2)/SM_obs

hist(t_boot_m)
# 2-sided p-value
(1 + sum(abs(t_boot_m) > abs(tm_obs)))/(1 + B)

par(mfrow = c(1, 2))
hist(t_boot, main = "Based on t-test statistic")
hist(t_boot_m, main = "Based on robust t-test statistic")
par(mfrow = c(1, 1))

# 3 : Another alternative
  # Function to calculate the test statistic t_M
calculate_tM <- function(group1, group2) {
  M1 <- median(group1)
  M2 <- median(group2)
  SM <- sum(abs(group1 - M1)) + sum(abs(group2 - M2))
  tM <- (M1 - M2) / SM
  return(tM)
}

  # Split the data into two groups
group1 <- sleep$extra[sleep$group == 1]
group2 <- sleep$extra[sleep$group == 2]

  # Calculate observed t_M
observed_tM <- calculate_tM(group1, group2)
print(observed_tM)

# Set seed for reproducibility
set.seed(123)

# Number of bootstrap samples
B <- 1000

# Combine data under null hypothesis (no difference in medians)
combined_data <- c(group1, group2)
n1 <- length(group1)
n2 <- length(group2)

# Generate bootstrap samples and calculate t_M
bootstrap_tM <- numeric(B)
for (b in 1:B) {
  #sample_data <- sample(combined_data, size = n1 + n2, replace = TRUE)
  sample_data <- rnorm(length(combined_data), mean(combined_data), sd(combined_data)^2)
  boot_group1 <- sample_data[1:n1]
  boot_group2 <- sample_data[(n1 + 1):(n1 + n2)]
  bootstrap_tM[b] <- calculate_tM(boot_group1, boot_group2)
}

# Plot the bootstrap distribution
hist(bootstrap_tM, main = "Bootstrap Distribution of t_M",
     xlab = "t_M", col = "lightblue", breaks = 20)
abline(v = observed_tM, col = "red", lwd = 2, lty = 2)

# Calculate p-value
p_value <- (1 + sum(abs(bootstrap_tM) >= abs(observed_tM)))/(1 + B)
print(p_value)


# 4

# Plot both distributions
par(mfrow = c(1, 2))
hist(t_boot, freq = FALSE, col = rgb(0, 0, 1, 0.5), breaks = 20,
     main = "Comparison of Bootstrap Distributions", xlab = "Test Statistic (Q.2.2)")
abline(v = t_obs, col = "red", lwd = 2, lty = 2)

hist(bootstrap_tM, freq = FALSE, col = rgb(1, 0, 0, 0.5), breaks = 20, add = FALSE,xlab = "Test Statistic (Q.2.3)", main = "")

abline(v = observed_tM, col = "blue", lwd = 2, lty = 2)
par(mfrow = c(1, 1))

# Numerical Comparison
   # --- Compare Mean and Variance of Distributions
# Bootstrap mean and SD
tM_bootstrap_mean <- mean(bootstrap_tM)
tM_bootstrap_sd <- sd(bootstrap_tM)

# Theoretical t-distribution mean and SD
t_bootstrap_mean <- mean(t.boot)
t_bootstrap_sd <- sd(t.boot)

# Output comparison
cat("Bootstrap t_M Mean:", tM_bootstrap_mean, "\n")
cat("Bootstrap t_M SD:", tM_bootstrap_sd, "\n")
cat("Theoretical t-distribution Mean:", t_bootstrap_mean, "\n")
cat("Theoretical t-distribution SD:", t_bootstrap_sd, "\n")

  #-----Compare Tail Behavior
# tM Bootstrap tail proportion
tM_bootstrap_tail <- mean(abs(bootstrap_tM) > 2)

# t Bootstrap tail proportion
t_bootstrap_tail <- mean(abs(t.boot) > 2)

cat("Bootstrap Tail Proportion:", tM_bootstrap_tail, "\n")
cat("Theoretical Tail Proportion:", t_bootstrap_tail, "\n")





#---------------------  Question 3

ID <- c(1:10)
x1 <- c(0.8, -1.23, 1.25, -0.28, -0.03, 0.61, 1.43, -0.54, -0.35, -1.60)
x2 <- c(0.64, -1.69, 1.47, -0.14, -0.18, 0.43, 1.61, -0.31, -0.38, -1.82)
df <- data.frame(ID, x1, x2)

# 1
(ration_obs <- mean(df$x1)/mean(df$x2))

# Bootstrap
B <- 1000
ratio_hat_boot <- c()
n <- length(df$x1)
index <- c(1:n)
set.seed(2023) 
for (i in 1:B) {
  index.b <- sample(index, n, replace = TRUE)
  x1_boot <- df$x1[index.b]
  x2_boot <- df$x2[index.b]
  ratio_hat_boot[i] <- mean(x1_boot)/mean(x2_boot)
}

# Jackknife
ratio_hat_jn <- c()
set.seed(2023)
for (i in 1:n) {
  x1_jn <- df$x1[-i]
  x2_jn <- df$x2[-i]
  ratio_hat_jn[i] <- mean(x1_jn)/mean(x2_jn)
}

# SE
sd(ratio_hat_boot) # Bootstrap
(sum((ratio_hat_jn - mean(ratio_hat_jn))^2 * (n - 1)/n))^(0.5) # Jackknife


sqrt(var(ratio_hat_jn)*(n - 1))

##

# Ratio statistic with small constant to stabilize division
ratio_statistic <- function(x1, x2, epsilon = 1e-6) {
  mean(x1) / (mean(x2) + epsilon)
}

# Non-parametric bootstrap
bootstrap_ratio <- function(x1, x2, B, epsilon = 1e-6) {
  n <- length(x1)
  bootstrap_estimates <- numeric(B)
  
  for (b in 1:B) {
    indices <- sample(1:n, size = n, replace = TRUE)
    x1_sample <- x1[indices]
    x2_sample <- x2[indices]
    bootstrap_estimates[b] <- ratio_statistic(x1_sample, x2_sample, epsilon)
  }
  
  return(bootstrap_estimates)
}



# Bootstrap SE for different B values
bootstrap_B_values <- c(10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 7500, 10000)
bootstrap_SEs <- sapply(bootstrap_B_values, function(B) {
  bootstrap_estimates <- bootstrap_ratio(x1, x2, B)
  sd(bootstrap_estimates) # Bootstrap standard error
})




# 3
quantile(ratio_hat_boot, p = c(0.025, 0.975))


# 4
# Observed ratio of means
theta_obs <- mean(x1) / mean(x2)

# Adjust data under H0: theta = 1
x2_adjusted <- x2 * (mean(x1) / mean(x2))

# Bootstrap procedure
set.seed(2025)
B <- 10000
bootstrap_ratios <- numeric(B)

for (i in 1:B) {
  # Resample x1 and x2_adjusted
  x1_boot <- sample(x1, length(x1), replace = TRUE)
  x2_boot <- sample(x2_adjusted, length(x2_adjusted), replace = TRUE)
  
  # Compute the ratio for bootstrap sample
  bootstrap_ratios[i] <- mean(x1_boot) / mean(x2_boot)
}

# One-sided p-value
p_value <- (1 + sum(bootstrap_ratios >= theta_obs))/(1 + B)
print(p_value)

t.test(ratio_hat_boot ~ 1, alternative = "greater")?
t.test(ratio_hat_boot ~ 1, alternative = "less")?


#----------------------- Question 4
  
# 1:  95% C.I for mu_d using the classical method
n <-length(df$x2)
x <- df$x1 - df$x2
mu_d <- mean(x)
se <- sd(x)/sqrt(10)
c(mu_d - qt(0.975, df = 9)*se, mu_d + qt(0.975, df = 9)*se)


#- Another way
# Mean difference
mu_d_hat <- mean(x1) - mean(x2)

# Standard error
s1 <- sd(x1)
s2 <- sd(x2)
n <- length(x1)
SE_mu_d <- sqrt((s1^2 / n) + (s2^2 / n))

# Critical value (t*)
t_critical <- qt(0.975, df = n - 1)

# Confidence interval
CI_lower <- mu_d_hat - t_critical * SE_mu_d
CI_upper <- mu_d_hat + t_critical * SE_mu_d

#-


# 2: 
# Non-parametric bootstrap
B <- 1000
mu_hat_boot <- c()
n <- length(df$x1)
index <- c(1:n)
set.seed(2023) 
for (i in 1:B) {
  index.b <- sample(index, n, replace = TRUE)
  x1_boot <- df$x1[index.b]
  x2_boot <- df$x2[index.b]
  mu_hat_boot[i] <- mean(x1_boot) - mean(x2_boot)
}
# 95%CI (percentile)
quantile(mu_hat_boot, prob = c(0.025, 0.975))


# BCa CI
bcboot <- function(x, nboot, theta, alpha = c(0.025, 0.975)) 
{ 
  n <- length(x)
  thetahat <- theta(x)
  bootsam <- matrix(sample(x, size = n*nboot, replace = TRUE), nrow = nboot)
  thetastar <- apply(bootsam, 1, theta)
  z0 <- qnorm(sum(thetastar < thetahat)/nboot)
  u <- rep(0, n)
  for(i in 1:n){
    u[i] <- theta(x[-i])
  }
  uu <- mean(u)-u
  acc <- sum(uu^3)/(6*(sum(uu^2))^1.5)
  zalpha <- qnorm(alpha)
  tt <- pnorm(z0 + (z0 + zalpha)/(1-acc*(z0 + zalpha)))
  confpoints <- quantile(x = thetastar, probs = tt, type = 1)
  return(confpoints)
}

x <- df$x1 - df$x2
set.seed(2023)
bcboot(x, 1000, mean, alpha=c(0.025, 0.975))

 
# Bootstrap t-Method
# Calculate bootstrap standard errors
bootstrap_se <- numeric(B)
for (i in 1:B) {
  indices <- sample(1:n, n, replace = TRUE)
  x1_boot <- x1[indices]
  x2_boot <- x2[indices]
  bootstrap_se[i] <- sqrt(var(x1_boot) / n + var(x2_boot) / n)
}

# Compute t-statistics
t_statistics <- (bootstrap_mu_d - mu_d_hat) / bootstrap_se
t_crit <- quantile(t_statistics, c(0.025, 0.975))
t_CI <- mu_d_hat - t_crit * sd(bootstrap_mu_d)



# 3
# Non-parametric bootstrap
B <- 1000
t_boot <- c()
n <- length(df$x1)
x <- c(df$x1, df$x2)
set.seed(2023) 
for (i in 1:B) {
  x_boot <- sample(x, 2*n, replace = T)
  x1 <- x_boot[1:n]
  x2 <- x_boot[(n+1):(2*n)]
  t_boot[i] <- t.test(x1, x2)$statistic
}
hist(t_boot)
t_obs <- t.test(df$x1, df$x2)$statistic
(1 + sum(abs(t_boot) > t_obs))/(1 + B)



#---------- Project 3
#===============================================================================
# Question 1
#------ 1
head(chickwts)
table(chickwts$feed)
lm_mod <- lm(weight ~ feed, data = chickwts)
anova(lm_mod)
Ftest_obs <- anova(lm_mod)$`F value`[1]
# The critical F-value at a significance level of α=0.05, with 5 and 65 degrees of freedom, is approximately 2.36
Ftest_obs=15.3648
Fcritical=2.36
Ftest_obs>Fcritical # => RH0

#------ 2
  # non-parametric: not correct
Ftest_boot <- c()
B <- 1000
n <- nrow(chickwts)
set.seed(2023)
for (i in 1:B) {
  weight_boot <- sample(chickwts$weight, n, replace = T)
  lm_boot <- lm(weight_boot ~ chickwts$feed)
  Ftest_boot[i] <- anova(lm_boot)$`F value`[1]
}

png("Figures/Fig1_Fnull.png", units="in", width = 9, height = 5, res = 300)
hist(Ftest_boot, xlim = c(0, 20), 
     main = "Distribution of bootstrap F statistics under the null",
     xlab = "F statistics")
abline(v = Ftest_obs, col = "red")
dev.off()

# P-value
(1 + sum(Ftest_boot > Ftest_obs))/(B + 1)

 # semi-parametric
fit.lm <- lm(chickwts$weight ~ 1) 
ei <- fit.lm$resid
n <- length(chickwts$weight)
B <- 1000
F_values <- c(1:B)
for(i in 1:B) {
  e.boot <- sample(ei, size = n, replace = T)
  chickwts.boot <- chickwts
  chickwts.boot$weight <- fitted(fit.lm) + e.boot
  boot_anova <- aov(weight ~ feed, data = chickwts.boot)
  F_values[i] <- summary(boot_anova)[[1]]["feed", "F value"]
}

hist(F_values,nclass=50)

png("Figures/Fig1_Fnull.png", units="in", width = 9, height = 5, res = 300)
hist(F_values, xlim = c(0, 20), 
     main = "Distribution of bootstrap F statistics under the null",
     xlab = "F statistics")
abline(v = observed_F, col = "red")
dev.off()
  
# Observed F-statistic
observed_F <- anova(lm_mod)$`F value`[1]

  # P-value
(1 + sum(F_values > observed_F))/(B + 1)



#------ 3 
Ftest_per <- c()
B <- 1000
n <- nrow(chickwts)
# Index
nindex <- table(chickwts$feed)
for (i in 1:B) {
  # Permute the feed labels
  permuted_data <- chickwts
  permuted_data$feed <- sample(permuted_data$feed)
  
  # Perform ANOVA on permuted data
  permuted_anova <- aov(weight ~ feed, data = permuted_data)
  Ftest_per[i] <- summary(permuted_anova)[[1]]["feed", "F value"]
}

hist(Ftest_per,nclass=50)

png("Figures/Fig1b_Fpernmu.png", units="in", width = 9, height = 5, res = 300)
hist(Ftest_per, xlim = c(0, 20), 
     main = "Distribution of permutation F statistics under the null",
     xlab = "F statistics")
abline(v = Ftest_obs, col = "red")
dev.off()

# P-value
p_value_permutation <-(1 + sum(Ftest_per > Ftest_obs))/(B + 1)
p_value_permutation


#------ 4: Semi-parametric
chickwts2 <- chickwts %>% filter(feed %in% c("sunflower", "soybean")) %>%
  mutate(feed01 = ifelse(feed == "sunflower", 1, 0))
lm_mod2 <- lm(weight ~ feed01, data = chickwts2)
lm_res <- lm_mod2$residuals
B <- 1000
n <- nrow(chickwts2)
theta_boot <- c()
set.seed(2023)
for (i in 1:B) {
  e_boot <- sample(lm_res, n, replace = T)
  wt_boot <- coef(lm_mod2)[1] + coef(lm_mod2)[2]*chickwts2$feed01 + e_boot
  lm_boot <- lm(wt_boot ~ chickwts2$feed01)
  theta_boot[i] <- coef(lm_boot)[2]
}

hist(theta_boot)

png("Figures/F1_4Semi.png", units="in", width = 9, height = 5, res = 300)
hist(theta_boot, 
     main = "Distribution of non-parametric bootstrap for theta",
     xlab = "Theta")
abline(v = quantile(theta_boot, probs = 0.05), col = "red", lty = 2)
abline(v = quantile(theta_boot, probs = 0.95), col = "red", lty = 2)
dev.off()

# 90%CI
quantile(theta_boot, probs = c(0.05, 0.95))


# Parametric
 # Subset data for Sunflower and Soybean
sunflower <- chickwts$weight[chickwts$feed == "sunflower"]
soybean <- chickwts$weight[chickwts$feed == "soybean"]

 # Estimate the observed difference in means
theta_hat <- mean(sunflower) - mean(soybean)
theta_hat

set.seed(123) 
B <- 1000
n <- nrow(chickwts2)
bootstrap_thetas <- c()

 # Parameters for Sunflower and Soybean groups
n_sunflower <- length(sunflower)
n_soybean <- length(soybean)
mean_sunflower <- mean(sunflower)
mean_soybean <- mean(soybean)
sd_sunflower <- sd(sunflower)
sd_soybean <- sd(soybean)

 # Generate bootstrap samples
for (i in 1:B) {
  # Resample weights under normal distribution
  boot_sunflower <- rnorm(n_sunflower, mean = mean_sunflower, sd = sd_sunflower)
  boot_soybean <- rnorm(n_soybean, mean = mean_soybean, sd = sd_soybean)
  
  # Calculate bootstrap estimate of theta
  bootstrap_thetas[i] <- mean(boot_sunflower) - mean(boot_soybean)
}

mean(bootstrap_thetas)
CI_lower <- quantile(bootstrap_thetas, 0.05)
CI_upper <- quantile(bootstrap_thetas, 0.95)
CI <- c(CI_lower, CI_upper)
CI



# Question 2
library(Ecdat)
data("Computers")
names(Computers)
options(scipen = 999)
#---------- 1
ols_mod <- lm(price ~ hd, data = Computers)
summary(ols_mod)

#---------- 2

price_pred <- predict(ols_mod)
hist(price_pred, nclass=50,
     main = "The distribution of prediction",
     xlab = "Predictions")

# Error
pred_error <- Computers$price - price_pred
png("Figures/F2_2error.png", units="in", width = 9, height = 5, res = 300)
hist(pred_error, nclass=50,
     main = "The distribution of errors",
     xlab = "Errors")
dev.off()


# Mean squared error
(mse1 <- mean(pred_error^2))
# Root Mean squared error
(rmse1 <- sqrt(mse1))



#---------- 3: 10 fold cross validation
  #Perform 10-Fold Cross-Validation
library(caret)
# Define the cross-validation control
set.seed(123)  # For reproducibility
cv_control <- trainControl(method = "cv", number = 10)
# Train the model using cross-validation
cv_model <- train(price ~ hd, data = Computers, method = "lm", trControl = cv_control)

# Print the results
print(cv_model)

#-------- Manually
# Manually perform 10-fold cross-validation
set.seed(123)
folds <- createFolds(Computers$price, k = 10, list = TRUE)

cv_errors <- sapply(folds, function(test_indices) {
  # Split data into training and test sets
  train_data <- Computers[-test_indices, ]
  test_data <- Computers[test_indices, ]
  
  # Fit the model on training data
  model <- lm(price ~ hd, data = train_data)
  
  # Predict on test data
  predictions <- predict(model, newdata = test_data)
  
  # Compute Mean Squared Error
  mean((predictions - test_data$price)^2)
})

# Compute average RMSE across folds
cv_rmse <- sqrt(mean(cv_errors))
cat("10-Fold Cross-Validated RMSE:", cv_rmse, "\n")



#---------- 4 : leave one out cross validation
n <- nrow(Computers)
price_pred_cv <- c() # Q3
beta1_cv <- c() # Q4
for (i in 1:n) {
  price_cv <- Computers$price[-i]
  hd_cv <- Computers$hd[-i]
  ols_cv <- lm(price_cv ~ hd_cv)
  price_pred_cv[i] <- coef(ols_cv)[1] + coef(ols_cv)[2]*Computers$hd[i]
 # beta1_cv[i] <- coef(ols_cv)[2]
}

pred_error_cv <- Computers$price - price_pred_cv

png("Figures/F2_2bLOOerror.png", units="in", width = 9, height = 5, res = 300)
hist(pred_error_cv, 
     main = "The distribution of LOO-CV errors",
     xlab = "LOO-CV errors")
dev.off()

(mse_cv <- mean(pred_error_cv^2)); mse1
(rmse_cv <- sqrt(mse_cv)); rmse1

#
hist(beta1_cv)
png("Figures/F2_3beta1.png", units="in", width = 9, height = 5, res = 300)
plot(1:n, beta1_cv, type = "l", xlab = "Iteration", ylab = "Beta 1")
dev.off()
#-
mean(price_pred)

png("Figures/F2_4pred.png", units="in", width = 9, height = 5, res = 300)
par(mfrow = c(1, 2))
hist(price_pred, 
     main = "Predicted price",
     xlab = "Predicted price")
hist(price_pred_cv, 
     main = "LOO-CV predicted price",
     xlab = "LOO-CV predicted price")
par(mfrow = c(1, 1))
dev.off()

quantile(price_pred, probs = c(0.025, 0.975))
quantile(price_pred_cv, probs = c(0.025, 0.975))


#---------- other method
# Initialize vector to store slope estimates
n <- nrow(Computers)
beta1_estimates <- numeric(n)

# Loop through each observation
for (i in 1:n) {
  # Leave one observation out
  train_data <- Computers[-i, ]
  
  # Fit the linear regression model
  model <- lm(price ~ hd, data = train_data)
  
  # Store the slope estimate
  beta1_estimates[i] <- coef(model)["hd"]
}

# Visualize the change in slope estimates
plot(
  1:n, beta1_estimates, type = "l", col = "blue", lwd = 2,
  xlab = "Excluded Observation Index",
  ylab = "Slope Estimate (β̂1)",
  main = "Change in Slope Estimate (β̂1) with LOOCV")

# Add a horizontal line for the full model slope estimate
full_model <- lm(price ~ hd, data = Computers)
abline(h = coef(full_model)["hd"], col = "red", lty = 2, lwd = 2)

legend(
  "topright", legend = c("LOOCV Slope Estimates", "Full Model Slope Estimate"),
  col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 2)
)


#----------------- 5
set.seed(123)  # For reproducibility
B <- 1000      # Number of bootstrap iterations
n <- nrow(Computers)  # Number of observations

# Specify the range of 'hd' values to predict
hd_values <- seq(min(Computers$hd), max(Computers$hd), length.out = 100)

# Matrix to store predictions for each bootstrap sample
bootstrap_predictions <- matrix(NA, nrow = B, ncol = length(hd_values))

# Perform bootstrap resampling
for (b in 1:B) {
  # Resample data with replacement
  boot_sample <- Computers[sample(1:n, replace = TRUE), ]
  
  # Fit the linear model
  model <- lm(price ~ hd, data = boot_sample)
  
  # Predict the values for the specified range of 'hd'
  bootstrap_predictions[b, ] <- predict(model, newdata = data.frame(hd = hd_values))
}

# Calculate the 95% confidence intervals
lower_bound <- apply(bootstrap_predictions, 2, quantile, probs = 0.025)
upper_bound <- apply(bootstrap_predictions, 2, quantile, probs = 0.975)
predicted_mean <- apply(bootstrap_predictions, 2, mean)

# Plot the regression line and confidence intervals
plot(Computers$hd, Computers$price, col = "gray", pch = 20,
     xlab = "Hard Drive Size (MB)", ylab = "Price (USD)",
     main = "Regression Line with 95% Bootstrap C.I.")
lines(hd_values, predicted_mean, col = "blue", lwd = 2, lty = 1)
lines(hd_values, lower_bound, col = "red", lwd = 2, lty = 2)
lines(hd_values, upper_bound, col = "red", lwd = 2, lty = 2)

legend(
  "topright",
  legend = c("Mean Prediction", "95% C.I."),
  col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 2)
)




# Question 3

#------------- 1
set.seed(123)  # For reproducibility
B <- 1000      # Number of bootstrap iterations
n <- nrow(Computers)

# Vectors to store bootstrap SEs
se_beta0 <- numeric(B)
se_beta1 <- numeric(B)

# Perform bootstrap
for (b in 1:B) {
  # Resample the data with replacement
  boot_sample <- Computers[sample(1:n, replace = TRUE), ]
  
  # Fit the linear regression model
  model <- lm(price ~ hd, data = boot_sample)
  
  # Extract standard errors of the coefficients
  se_beta0[b] <- summary(model)$coefficients[1, 2]
  se_beta1[b] <- summary(model)$coefficients[2, 2]
}

# Calculate 95% confidence intervals
ci_beta0 <- quantile(se_beta0, probs = c(0.025, 0.975))
ci_beta1 <- quantile(se_beta1, probs = c(0.025, 0.975))

# Print results
cat("95% C.I. for SE(β̂0):", ci_beta0, "\n")
cat("95% C.I. for SE(β̂1):", ci_beta1, "\n")



#------------------ 2
set.seed(123)  # For reproducibility
B <- 1000      # Number of bootstrap iterations
n <- nrow(Computers)

# Split data into full dataset and reduced dataset (exclude hd > 2000)
full_data <- Computers
reduced_data <- Computers[Computers$hd <= 2000, ]

# Function to perform bootstrap and calculate SEs
bootstrap_se <- function(data, B) {
  se_beta0 <- numeric(B)
  se_beta1 <- numeric(B)
  
  for (b in 1:B) {
    boot_sample <- data[sample(1:nrow(data), replace = TRUE), ]
    model <- lm(price ~ hd, data = boot_sample)
    se_beta0[b] <- summary(model)$coefficients[1, 2]
    se_beta1[b] <- summary(model)$coefficients[2, 2]
  }
  
  return(list(se_beta0 = se_beta0, se_beta1 = se_beta1))
}

# Perform bootstrap for full and reduced datasets
bootstrap_full <- bootstrap_se(full_data, B)
bootstrap_reduced <- bootstrap_se(reduced_data, B)

# Calculate 95% confidence intervals
ci_full_beta0 <- quantile(bootstrap_full$se_beta0, probs = c(0.025, 0.975))
ci_full_beta1 <- quantile(bootstrap_full$se_beta1, probs = c(0.025, 0.975))
ci_reduced_beta0 <- quantile(bootstrap_reduced$se_beta0, probs = c(0.025, 0.975))
ci_reduced_beta1 <- quantile(bootstrap_reduced$se_beta1, probs = c(0.025, 0.975))

# Print results
cat("95% C.I. for SE(β̂0) - Full Data:", ci_full_beta0, "\n")
cat("95% C.I. for SE(β̂1) - Full Data:", ci_full_beta1, "\n")
cat("95% C.I. for SE(β̂0) - Reduced Data:", ci_reduced_beta0, "\n")
cat("95% C.I. for SE(β̂1) - Reduced Data:", ci_reduced_beta1, "\n")

# Plot bootstrap distributions for SE(β̂0)
par(mfrow = c(1, 2))
hist(bootstrap_full$se_beta0, col = "blue", main = "SE(β̂0) - Full Data",
     xlab = "SE(β̂0)", breaks = 20)
hist(bootstrap_reduced$se_beta0, col = "red", main = "SE(β̂0) - Reduced Data",
     xlab = "SE(β̂0)", breaks = 20)

# Plot bootstrap distributions for SE(β̂1)
hist(bootstrap_full$se_beta1, col = "blue", main = "SE(β̂1) - Full Data",
     xlab = "SE(β̂1)", breaks = 20)
hist(bootstrap_reduced$se_beta1, col = "red", main = "SE(β̂1) - Reduced Data",
     xlab = "SE(β̂1)", breaks = 20)




# Question 4: 

x<-c(0.68446806,-0.02596037,-0.90015774,0.72892605,-0.45612255, 0.19311847, -0.13297109, -0.99845382,  0.37278006, -0.20371894, -0.15468803,  0.19298230
     , -0.42755534, -0.04704525,  0.15273726, 0.03655799, 0.01315016, -0.59121428, 
     4.50955771, 2.87272653) 
length(x) 

#---------- 1
mu_x <- mean(x)
med_x <- median(x)

#---------- 2
B <- 1000
mu_boot <- c()
med_boot <- c()
n <- length(x)
for (i in 1:B) {
  x_boot <- sample(x, n, replace =  T)
  mu_boot[i] <- mean(x_boot)
  med_boot[i] <- median(x_boot)
}
hist(mu_boot)
hist(med_boot)

png("Figures/F4_2Dis.png", units="in", width = 9, height = 5, res = 300)
par(mfrow = c(1, 2))
hist(mu_boot, 
     main = "Distribution of bootstrap mean",
     xlab = "Bootstrap mean")
hist(med_boot, 
     main = "Distribution of bootstrap median",
     xlab = "Bootstrap median")
par(mfrow = c(1, 1))
dev.off()


#---------- 3
# Semi-parametric bootstrap
mod_mean <- lm(x ~ 1)
mod_res <- residuals(mod_mean)
n <- length(x)
B <- 1000
mu_boot_semi <- c()
med_boot_semi <- c()
set.seed(2023)

for (i in 1:B) {
  res_boot <- sample(mod_res, n, replace = T)
  x_boot <- coef(mod_mean) + res_boot
  mu_boot_semi[i] <- mean(x_boot)
  med_boot_semi[i] <- median(x_boot)
}

hist(mu_boot_semi)
hist(med_boot_semi)
quantile(mu_boot_semi, probs = c(0.025, 0.975))
quantile(med_boot_semi, probs = c(0.025, 0.975))


(se_mu <- sd(mu_boot_semi))
(se_med <- sd(med_boot_semi))

#---------- 4
# Bootstrap
B <- 1000
mu_boot <- c()
med_boot <- c()
n <- length(x)
set.seed(2023)
for (i in 1:B) {
  x_boot <- sample(x, n, replace =  T)
  mu_boot[i] <- mean(x_boot)
  med_boot[i] <- median(x_boot)
}

# MSE for mean
mu_error <- mu_x - mu_boot
(mu_sme <- mean(mu_error^2))
# MSE for median
med_error <- med_x - med_boot
(med_sme <- mean(med_error^2))

# Jackknife
B <- 1000
mu_jn <- c()
med_jn <- c()
n <- length(x)

for (i in 1:n) {
  x_jn <- x[-i]
  mu_jn[i] <- mean(x_jn)
  med_jn[i] <- median(x_jn)
}

# MSE for mean (adjust for inflation factor)
mu_er_jn <- mu_x - mu_jn
(mu_sme_jn <- mean(mu_er_jn^2)*(n-1))
# MSE for median (adjust for inflation factor)
med_er_jn <- med_x - med_jn
(med_sme_jn <- mean(med_er_jn^2)*(n-1))


   #------------------ other options
n <- length(x)
# Jackknife estimates for mean and median
jackknife_means <- numeric(n)
jackknife_medians <- numeric(n)

for (i in 1:n) {
  x_loo <- x[-i]  # Leave-one-out sample
  jackknife_means[i] <- mean(x_loo)
  jackknife_medians[i] <- median(x_loo)
}

# MSE calculations
mean_mse <- (n - 1) / n * sum((jackknife_means - mean(jackknife_means))^2)
median_mse <- (n - 1) / n * sum((jackknife_medians - mean(jackknife_medians))^2)

# Results
cat("MSE for Mean:", mean_mse, "\n")
cat("MSE for Median:", median_mse, "\n")

# Decide which parameter to prefer
if (mean_mse < median_mse) {
  cat("The Mean is preferred based on lower MSE.\n")
} else {
  cat("The Median is preferred based on lower MSE.\n")
}



#---------- 5
set.seed(123)  
n <- length(x)
B <- 1000 

# Initialize vectors
bootstrap_medians <- numeric(B)  
medians_below_zero <- numeric(0)  

for (b in 1:B) {
  resample <- sample(x, size = n, replace = TRUE)  
  bootstrap_median <- median(resample)            
  bootstrap_medians[b] <- bootstrap_median        
  if (bootstrap_median < 0) {
    medians_below_zero <- c(medians_below_zero, bootstrap_median)  
  }
}

# Calculate the proportion of medians < 0
pi_hat <- mean(bootstrap_medians < 0)

# Distribution of pi
hist(medians_below_zero, 
     main = "Distribution of pi_M",
     xlab = "Bootstrap pi_M")
par(mfrow = c(1, 1))
dev.off()

# 95% Confidence Interval
ci <- quantile(medians_below_zero, probs = c(0.025, 0.975))

# Results
cat("Estimated π(M < 0):", pi_hat, "\n")
cat("95% Confidence Interval for π(M < 0):", ci, "\n")















