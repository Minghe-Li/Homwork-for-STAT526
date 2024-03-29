library(faraway)
attach(salmonella)
salmonella
plot(salmonella)
Means <- vector()
Vars <- vector()
for (i in 0:5) {
  Means[i]=mean(salmonella[(i*3+1):(i*3+3),1])
  Vars[i]=var(salmonella[(i*3+1):(i*3+3),1])
}
Ratio = Vars/Means
plot(Vars,Means)
mean(salmonella[(16):(18),1])
var(salmonella[(16):(18),1])
Means = c(21.667, 18.333, 25, 42.667, 37.333, 29.667)
Vars = c(49.333, 6.333, 73, 274.333, 16.333, 126.333)
doselevel = factor(dose)
model<-glm(colonies~doselevel,family = poisson,salmonella)
summary(model)
pchisq(33.496, df=12, lower.tail=FALSE)
model2<-glm(colonies~log(dose+1),family = poisson,salmonella)
summary(model2)
pchisq(59.629, df=16, lower.tail=FALSE)
plot(model)
(dp <- sum(residuals(model,type="pearson")^2)/model$df.res)
sumary(model,dispersion=dp)
Predicted = fitted(model)
plot(colonies,Predicted)
plot(colonies,dose,col="red")
par(new=TRUE)
plot(colonies,Predicted,col="blue")
model3<-glm(colonies~dose,family = poisson,salmonella)
summary(model3)
summary(model3)$cov.unscaled

wavesolder
attach(wavesolder)
plot(y1,prebake,col="red")
par(new=TRUE)
plot(y2,prebake,col="blue")
par(new=TRUE)
plot(y3,prebake,col="yellow")
mean(wavesolder[,3])
var(wavesolder[,3])
Means = c(25.9375, 38.1875, 32.625)
Vars = c(578.3292, 2279.496, 1447.45)
plot(log(Means),log(Vars))
wavesolder11 = wavesolder[,-2]
wavesolder1 = wavesolder11[,-2]
wavesolder1
model4<-glm(y1~.,family = poisson,wavesolder1)
summary(model4)
pchisq(23.198, df=8, lower.tail=FALSE)
plot(model4)
wavesolder1_new = wavesolder1[-5,]
model5<-glm(y1~.,family = poisson,wavesolder1_new)
summary(model5)
pchisq(41.455, df=61, lower.tail=FALSE)
(dp <- sum(residuals(model5,type="pearson")^2)/model5$df.res)
sumary(model5,dispersion=dp)
