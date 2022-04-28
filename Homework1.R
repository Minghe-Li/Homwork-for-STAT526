# Q1 
# Knowing n = 88, p-1 = 6 and R^2 = 0.48.
# df1 = 6 and  df2 = 81
# F = 12.462 So, p = 6.7*E-10

# Q2
# Full model Y[i] = Beta[0]+sum(Beta[j]*X[ij],1,p-1)+e[i]
# Degree of freedom = n-1, error df = n-p; 
# Reduced model Beta[1]=Beta[2]=...Beta[q]=0
# Y[i] = Beta[0]+sum(Beta[j]*X[ij],q+1,p-1)+e[i]
# Degree of freedom = n-2, error df = n-c

# Q3
# Because the proportion of undercounted votes to be more variable in smaller counties than larger ones, we can use ballots as weight.
library(faraway)
percunder <- (gavote$ballots - gavote$votes)/gavote$ballots
pergore = gavote$gore/gavote$votes
perbush = gavote$bush/gavote$votes
model1 = lm(percunder ~ equip + econ + rural + perAA + equip:econ +equip:perAA, data = gavote)
model2 = lm(percunder ~ equip + econ + rural + perAA + equip:econ +equip:perAA, data = gavote, weights = ballots)
summary(model1)
summary(model2)
# explanation: The P-value of the weighted model is much smaller than unweighted.

# Q5
modelmax = lm(perm ~ (area+peri+shape)^2, rock)
modelbest = step(modelmax,trace=FALSE)
summary(modelbest)
drop1(modelbest,test="F")
# Final model : perm ~ area + peri + shape + area:shape + peri:shape