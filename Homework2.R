# Q1
library(faraway)
summary(attitude)
pairs(~attitude$rating+attitude$complaints+attitude$privileges+attitude$learning+attitude$raises+attitude$critical+attitude$advance,pch=20)
attach(attitude)
ratingc = (rating - mean(rating))/sd(rating)
complaintsc = (complaints - mean(complaints))/sd(complaints)
privilegesc = (privileges - mean(privileges))/sd(privileges)
learningc = (learning - mean(learning))/sd(learning)
raisesc = (raises - mean(raises))/sd(raises)
criticalc = (critical - mean(critical))/sd(critical)
advancec = (advance - mean(advance))/sd(advance)
modelmax = lm(rating ~ (complaintsc+privilegesc+learningc+raisesc+criticalc+advancec)^2, attitude)
modelbest = step(modelmax,trace=FALSE)
summary(modelbest)
drop1(modelbest,test="F")
# Final model : rating ~ complaints + learning + advance
boxplot(attitude)
cor(attitude)