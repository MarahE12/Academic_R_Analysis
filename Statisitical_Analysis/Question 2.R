#set working directory
setwd("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2")
#check files within working directory
dir()

#read in piwi.txt file
piwi<-read.table('piwi.txt', header = TRUE)
#look at size of the data
dim(piwi)
#look at header of the data
head(piwi)
#look at summary of the data
summary(piwi)
#check objects within data
str(piwi)

#create subset of only piwi mutants
piwi.mut<-subset(piwi, line == 'piwi')
dim(piwi.mut) #there are 20 mutants

#create subset of only wild-type
not.piwi<-subset(piwi, line == 'wt')

#find mean and standard deviation of transpositions for both mutants and wild-type
mean(piwi.mut$transpositions) #mean for mutant is 2.25
sd(piwi.mut$transpositions) #standard deviation for mutant is 1.517442
mean(not.piwi$transpositions) #mean for wild-type is 1.45
sd(not.piwi$transpositions) #standard deviation for wild-type is 1.356272

#load ggplot2
library(ggplot2)
#make line a factor to plot
piwi$line.fac<-factor(piwi$line, levels=c('wt', 'piwi'), labels=c('Wild-Type', 'Piwi Mutant'))
#plot transpositions against both wild-type and mutant
plot1<-ggplot(data=piwi, aes(y=transpositions, x=line.fac))+
  geom_boxplot()+ geom_point()+
  labs(title='Comparison of Transposition Distribution in Drosophila subobscura',
  y='Number of Transpositions', x='')+
  theme_classic()+
  theme(plot.title = element_text(size = 11, face = "bold", hjust=0.5))
print(plot1)
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 2 Boxplot.tiff", height=10, width=20, units='cm', compression='lzw', res=300)
plot(plot1)
dev.off()
#there is visual increase in transpositions in piwi mutant than wild-type


##
#conducting a t-test to test for a significant difference between piwi mutant and wild-type in 
#number of transpositions as t-test requires a normal data distribution, produce a histogram plot for 
#each genetic line, separately
#histogram for wild-type
WT_hist<-hist(not.piwi$transpositions, xlab='Number of Transpositions', ylab='Frequency', main='Distribution of Transpositions in Wild-Type')
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 2 WT Hist.tiff", height=10, width=15, units='cm', compression='lzw', res=300)
plot(WT_hist, xlab='Number of Transpositions', ylab='Frequency', main='Distribution of Transpositions in Wild-Type')
dev.off()

##histogram for piwi mutant
PM_hist<-hist(piwi.mut$transpositions, xlab='Number of Transpositions', ylab='Frequency', main='Distribution of Transpositions in Piwi Mutant')
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 2 PM Hist.tiff", height=10, width=15, units='cm', compression='lzw', res=300)
plot(PM_hist, xlab='Number of Transpositions', ylab='Frequency', main='Distribution of Transpositions in Piwi Mutant')
dev.off()

#nether wt or piwi data follows a normal distribution? 

#use a Wilcoxon-Mann-Whitney test as it does not assume normality and is a non-parametric 
#equivalent to t-test
wilcox.test(transpositions~line, data = piwi)

#wilcoxon test showed a no significance difference, to double check, conduct a GLM
stat<-glm(transpositions~line, data = piwi, family='poisson')
drop1(stat, test='Chisq')

      