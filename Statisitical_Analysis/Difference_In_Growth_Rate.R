##In this analysis, I was provided a dataset containing growth rate of water flea Daphnia pulex exposed to pollutants cadmium and glyphosate.
##I was tasked to explore the question: "Is there a difference in the growth rate between the groups?"

####

#set working directory
setwd("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2")
#check files within working directory
dir()

#read in Daphnia.txt file
daphnia<-read.table('Daphnia.txt', header = TRUE)
#look at size of the data
dim(daphnia)
#look at header of the data
head(daphnia)
#look at summary of the data
summary(daphnia)
#check objects within data
str(daphnia)

#use sub-setting, to identify number of individuals in each group
CD_sub<-subset(daphnia, trt == 'CD')
dim(CD_sub) #there are 10 individuals in the CD treatment
CTRL_sub<-subset(daphnia, trt == 'CTRL')
dim(CTRL_sub) #there are 10 individuals int he CTRL treatment
GLY_sub<-subset(daphnia, trt == 'GLY')
dim(GLY_sub) #there are 10 individuals in the GLY treatment

#find median of age at maturity in each treatment group
median(CD_sub$age_mat) #median is 7
median(CTRL_sub$age_mat) #median is 6
median(GLY_sub$age_mat) #median is 6

#find mean of growth rate in each treatment group
mean(CD_sub$GR_pre) #mean is 0.171
mean(CTRL_sub$GR_pre) #mean is 0.2355
mean(GLY_sub$GR_pre) #mean is 0.2457

#find standard deviation of growth rate in each treatment group
sd(CD_sub$GR_pre) #standard deviation is 0.008041559
sd(CTRL_sub$GR_pre) #standard deviation is 0.00643342
sd(GLY_sub$GR_pre) #standard deviation is 0.008743887

#load ggplot2
library(ggplot2)
#make the treatment group a factor to plot on a graph
daphnia$trt.fac<-factor(daphnia$trt, levels=c('CD', 'CTRL', 'GLY'), labels=c('0.05 μg L-1 Cadmium Chloride', 'Clean Water', '100
μg L-1 Glyphosate'))
#plot growth rate vs each treatment group using ggplot
GR_plot<-ggplot(data=daphnia, aes(y=GR_pre, x=trt.fac)) +
  geom_boxplot() + geom_point()
GR_plot<-GR_plot+
  labs(title='Distribution of Growth Rates in Treatment Groups',
       y='Growth rate (mm/day)', x='')+
  theme_classic()+
  theme(plot.title = element_text(size = 13, face = "bold", hjust=0.5))
print(GR_plot)
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 1 Boxplot.tiff", height=10, width=20, units='cm', compression='lzw', res=300)
plot(GR_plot)
dev.off()


#check for normal distribution
library(ggplot2)
GR_log<-log(daphnia$GR_pre)
plot_hist<-ggplot(data=daphnia, aes(x=GR_pre, fill=trt.fac))+
  geom_histogram(alpha=0.5)+
  theme_classic()+
  labs(title='Daphnia pulex Growth Rate in Different Treatments Exposures',
       x='Growth rate (mm/day)', y='Frequency')+
  theme(plot.title = element_text(size = 12, face = "bold", hjust=1.5))+
  scale_fill_manual(values = c("orange","lightblue","darkblue"), name = "Treatment Conditions")+
  facet_wrap(~trt.fac, ncol=4)
print(plot_hist)
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 1 Hist.tiff", height=10, width=20, units='cm', compression='lzw', res=300)
plot(plot_hist)
dev.off()

#carry out an ANOVA test to statistically test for a significance different among the treatment groups
growth.aov<-aov(GR_pre~trt, data=daphnia)
summary(growth.aov)



