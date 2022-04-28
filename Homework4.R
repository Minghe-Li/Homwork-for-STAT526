library(faraway)
attach(esoph)
esoph
model = glm(cbind(ncases,ncontrols)~.^2,family=binomial,esoph)
AIC(model)
model2 = glm(cbind(ncases,ncontrols)~.,family=binomial,esoph)
AIC(model2)
model3 = glm(cbind(ncases,ncontrols)~agegp+alcgp+tobgp+agegp:alcgp,family=binomial,esoph)
AIC(model3)
model4 = glm(cbind(ncases,ncontrols)~agegp+alcgp+tobgp+tobgp:alcgp,family=binomial,esoph)
AIC(model4)
model5 = glm(cbind(ncases,ncontrols)~agegp+alcgp+tobgp+agegp:tobgp,family=binomial,esoph)
AIC(model5)

esoph_unclass <- data.frame("agegp" = unclass(agegp), "alcgp" = unclass(alcgp), "tobgp" = unclass(tobgp), "ncases" = ncases, "ncontrols" = ncontrols)
modelunclass = glm(cbind(ncases,ncontrols)~.,family=binomial,esoph_unclass)
summary(modelunclass)
par(mfrow = c(2,2))
plot(modelunclass)
summary(modelunclass)$cov.unscaled
summary(model)
# aictab(model)
modelbest = step(model,trace=FALSE)
summary(modelbest)

pyrimidines
#pyrimidines1 = pyrimidines[-10,]
attach(pyrimidines)
plot(p1.polar,activity)
plot(p1.size,activity)
plot(p1.flex,activity)
pyrimidines1 = pyrimidines[-10,]
attach(pyrimidines1)
modelpyr = glm(activity~.,family=gaussian,pyrimidines1)
summary(modelpyr)
(1-exp((modelpyr$dev-modelpyr$null)/37))/(1-exp(-modelpyr$null/37))

gaussianpredict = predict(modelpyr)
modelpyr2 = glm(activity~.,family=quasibinomial,pyrimidines1)
summary(modelpyr2)
quasibinpredict = predict(modelpyr2)
plot(gaussianpredict,quasibinpredict)
modelpyr3 = glm(logit(activity)~.,family=gaussian,pyrimidines1)
summary(modelpyr3)
logitpredict = predict(modelpyr3)
plot(logitpredict,quasibinpredict)
#install.packages("nlme")
library(nlme)
library(mgcv)
# modelbeta <- gam(activity~., family=betar(), data = pyrimidines1)
modelbeta <- gam(activity~p1.polar+p1.size+p1.flex+p1.h.doner+p1.h.acceptor+p1.pi.doner+p1.pi.acceptor+p1.polarisable+p1.sigma+p2.polar+p2.size+p2.flex+p2.h.doner+p2.h.acceptor+p2.pi.doner+p2.pi.acceptor+p2.polarisable+p2.sigma+p3.polar+p3.size+p3.flex+p3.h.doner+p3.h.acceptor+p3.pi.doner+p3.polarisable+p3.sigma, family=betar(link="logit"),data=pyrimidines1)
summary(modelbeta)

# simulation
# number of smaples
nsample = 800
mtrails = 45
# initialize pearson chisquare and deviance
pchisq = numeric(length=1000)
dev = numeric(length=1000)
# Generate x from normal distribution
x <- rnorm(nsample, 0, 1)
# Generate N = 1000 data sets
for(i in 1:1000){
  logitp = 0.35+x
  # y|x ~ Bernoulli(p)
  y = rbinom(nsample,mtrails,(1+exp(-logitp))^(-1))
  model1 = glm(cbind(mtrails-y,y)~x,family=binomial)   
  pchisq[i] = sum(residuals(model1,type = "pearson")^2)
  dev[i] = summary(model1)$deviance
}
# plot result
df1=nsample-2
par(mfrow=c(1,2))
hist(pchisq/df1,prob=T,nclass=25,las=1)
lines(seq(.5,1.5,length=500),df1*dchisq(seq(0.5*df1,1.5*df1,length=500),nsample-2),main=str1)
hist(dev/df1,prob=T,nclass=25,las=1)
lines(seq(.5,1.5,length=500),df1*dchisq(seq(0.5*df1,1.5*df1,length=500),nsample-2),main=str2)


nsample = 800
mtrails = 800
pchisq = numeric(length=1000)
dev = numeric(length=1000)
x <- rnorm(nsample, 0, 1)
for(i in 1:1000){
  logitp = 0.35+x
  y = rbinom(nsample, mtrails,(1+exp(-logitp))^(-1))
  model1 = glm(cbind(mtrails-y,y)~x,family=binomial)   
  pchisq[i] = sum(residuals(model1,type = "pearson")^2,df.residual(model1),lower=FALSE)
  dev[i] = summary(model1)$deviance
}
df1=nsample-2
plot.new
par(mfrow=c(1,1))
str1=paste('chisquare',toString(mtrails),toString(nsample),'.jpeg')
jpeg(str1)
hist(pchisq/df1,prob=T,nclass=25,las=1)
lines(seq(.5,1.5,length=500),df1*dchisq(seq(0.5*df1,1.5*df1,length=500),nsample-2),main=str1)
dev.off()

plot.new()
str2=paste('hist',toString(mtrails),toString(nsample),'.jpeg')
jpeg(str2)
hist(dev/df1,prob=T,nclass=25,las=1)
lines(seq(.5,1.5,length=500),df1*dchisq(seq(0.5*df1,1.5*df1,length=500),nsample-2),main=str2)
dev.off()

