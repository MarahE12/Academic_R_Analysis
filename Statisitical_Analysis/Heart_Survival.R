##In this analysis, I was provided "heart_survival.csv" containing survival rates following a heart attack.
##I was tasked to perform survival analysis and identify the estimated survival probability of males and females.

####

#install packages
install.packages("dplyr")
install.packages("survminer")
library(dplyr)
library(survival)
library(survminer)

#Read data in 
heart <- read.csv('heart_survival.csv', header = TRUE)
summary(heart)
head(heart)

#question 2.1#
#fit survival curves separately by sex & produce a p value using surv_pvalue.
sfit <- survfit(Surv(time, event)~sex, data=heart)
sfit
summary(sfit)
surv_pvalue(sfit)

#look at the range of time
range(heart$time)
#1 - 2358 days

#create a sequence of numbers going from one number to another number by increments of yet another number
seq(0, 2358, 100)
summary(sfit, time=seq(0, 2358, 100))

#we have fit survival curve, we can now visualize it with a Kaplan-Meier plot
plot(sfit)
ggsurvplot(sfit, data = heart,
           ggtheme = theme_gray(),
           xlab = 'Time (days)', ylab='Survival Probability')

#females have a greater survival probability but men survive longer

#question 2.2#
#p value
surv_pvalue(sfit)

#question 2.3#
time_points <- c(500, 1000, 1500, 2000)
newdata <- data.frame(sex = c("Male", "Female"))
#summary to get probabilities at specific timepoints
surv_summary <- summary(sfit, times = time_points, newdata = newdata)

#survival probabilities
surv_prob <- surv_summary$surv
summary(surv_prob)
print(surv_prob)
