library(faraway)
summary(spector)
attach(spector)
pairs(~grade+psi+tuce+gpa)
model1 = glm(grade ~ psi+tuce+gpa, family=binomial, spector)
summary(model1)
model2 = step(model1,trace=0)
summary(model2)
#summary(grade)
#summary(psi)
# install.packages("statmod")
library(statmod)
library(boot)
glm.diag.plots(model2)

exp(coef(model2)[2]*1)
exp(confint(model2)[2,]*1)
table(grade,psi)
#install.packages("generalhoslem")
#install.packages("reshape")
#install.packages("Mass")
#library(reshape)
#library(MASS)
#library(generalhoslem)
#logitgof(grade,fitted(model2))
summary(orings)
# simulation
# number of smaples
nsample = 800
# initialize pearson chisquare and deviance
pchisq = numeric(length=1000)
dev = numeric(length=1000)
# Generate x from normal distribution
x <- rnorm(nsample, 0, 1)
# Generate N = 1000 data sets
for(i in 1:1000){
  logitp = 0.35+x
  # y|x ~ Bernoulli(p)
  y = rbinom(nsample,1,(1+exp(-logitp))^(-1))
  model1 = glm(y~x,family=binomial)   
  pchisq[i] = sum(residuals(model1,type = "pearson")^2)
  dev[i] = summary(model1)$deviance
}
# plot result
par(mfrow=c(1,1))
hist(pchisq,prob=T,nclass=25,ylim=c(0,.07),las=1)
lines(seq(min(pchisq),max(pchisq),length=500),dchisq(seq(min(pchisq),max(pchisq),length=500),nsample-2))
hist(dev,prob=T,nclass=25,ylim=c(0,.1),las=1)
lines(seq(min(dev),max(dev),length=500),dchisq(seq(min(dev),max(dev),length=500),nsample-2))
