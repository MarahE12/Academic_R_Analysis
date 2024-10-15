##In this analysis, I was provided "blood_pressure.csv " containing prevalence of raised blood pressure from different countries.
##I was tasked to 'provide a model of the association between Sex, Region and Prevalence of raised blood pressure'.

####

#read dataset in 
blood <- read.csv('blood_pressure.csv', header = TRUE)
head(blood)
summary(blood)

#factoring character variables
blood$Sex<-factor(blood$Sex)
blood$Region<-factor(blood$Region)

#Question 1.1#
#subset dataset to be male only
male_blood <- subset(blood, Sex == 'Men')
#identify country with highest prevalence
highest_country <- male_blood$Country[which.max(male_blood$Prevalence)]
#identify country with lowest prevalence
lowest_country <- male_blood$Country[which.min(male_blood$Prevalence)]
#print results
cat('highest: ', highest_country, 'lowest:', lowest_country)

#Question 1.2#
#plot association between Sex, Region and Prevalence 
library(ggplot2)
#boxplot
ggplot(blood, aes(x = Region, y = Prevalence, fill = Sex)) +
  geom_boxplot(color = "black", outlier.color = "red", width = 0.5) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +  
  labs(x = "Global Region", y = "Prevalence (Proportion)") +
  theme_gray() +
  theme(text = element_text(size = 12))

#Question 1.3#
#model selection
#check prevalence for a normal distribution
hist(blood$Prevalence)
#data shows a normal distribution, so will test for association significance using lm

#maximal model
blood.model <- lm(Prevalence ~ Sex + Region, data = blood)
summary(blood.model)
blood.model2 <- lm(Prevalence ~ Sex * Region, data = blood)
summary(blood.model2)
drop1(blood.model, test ='F')
#visualise
library(lattice)
xyplot(Prevalence ~ Region | Sex, data = blood)

anova(blood.model, blood.model2, test ='F')
#model 2 has a better fit.
#drop1 test
anova(blood.model2, test ='F')


#both sex and region show a significant association with prevalence
