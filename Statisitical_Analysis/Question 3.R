#set working directory
setwd("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2")
#check files within working directory
dir()

#read in mushroom_data.csv file
mushroom<-read.table('mushroom_data.csv', header = TRUE, sep=',')
#look at size of the data
dim(mushroom)
#look at header of the data
head(mushroom)
#look at summary of the data
summary(mushroom)
#check objects within data
str(mushroom)

##calculating statistical summary

#create subset of edible mushrooms only
edible<-subset(mushroom, class == 'e')
#create subset of poisonous mushrooms only
poison<-subset(mushroom, class == 'p')

#find mean and standard deviation of cap diameter and stem height of edible mushrooms
mean(edible$cap.diameter) #mean of cap diameter is 7.7 98696
sd(edible$cap.diameter) #standard deviation of cap diameter is 6.373404
mean(edible$stem.height) #mean of stem height is 7.039077
sd(edible$stem.height) #standard deviation of stem height is 3.583433

#find mean and standard deviation of cap diameter and stem height of poisonous mushrooms
mean(poison$cap.diameter) #mean of cap diameter is 5.879763
sd(poison$cap.diameter) #standard deviation of cap diameter is 3.966391
mean(poison$stem.height) #mean of stem height is 6.214554
sd(poison$stem.height) #standard deviation of stem height is 3.140778

###creating visual presentation for the cap diameter
#load ggplot2
library(ggplot2)
#make class a factor to plot
mushroom$class.fac<-factor(mushroom$class, levels=c('e', 'p'), labels=c('Edible', 'Poisonous'))

#convert the cap diameters and stem height into log for better visualisation
mushroom$cap.log<-log(mushroom$cap.diameter)
mushroom$stem.log<-log(mushroom$stem.height)

#generate a boxplot to visualise data
plot_box<-ggplot(data=mushroom, aes(y=cap.diameter, x=class.fac))+
  geom_boxplot()+
  labs(title='Cap Diameter: Edible vs. Poisonous Mushrooms',
       y='Mushroom Cap Diameter (cm)', x='')+
  theme_linedraw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust=0.5))
print(plot_box)
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 3 Cap Boxplot.tiff", height=10, width=15, units='cm', compression='lzw', res=300)
plot(plot_box)
dev.off()

plot_box2<-ggplot(data=mushroom, aes(y=stem.height,x=class.fac))+
  geom_boxplot()+
  labs(title='Stem Height: Edible vs. Poisonous Mushrooms',
       y='Mushroom Stem Height (cm)', x='')+
  theme_linedraw()+
  theme(plot.title = element_text(size = 12, face = "bold", hjust=0.5))
print(plot_box2)
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 3 Stem Boxplot.tiff", height=10, width=15, units='cm', compression='lzw', res=300)
plot(plot_box2)
dev.off()

#plot histogram of cap diameter for both edible and poisonous mushrooms
plot_hist<-ggplot(data=mushroom, aes(x=cap.log, fill=class.fac))+
  geom_histogram(bins=30, alpha=0.5)+
  labs(title='Cap Diameter Comparison: Edible vs. Poisonous Mushrooms',
       x='Mushroom Cap Diameter (cm)', y='Frequency')+
  theme_classic()+
  scale_fill_manual(values = c("grey","red"), name = "Mushroom Type")+
  theme(plot.title = element_text(size = 12, face = "bold"))
print(plot_hist)
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 3 Cap Hist.tiff", height=10, width=20, units='cm', compression='lzw', res=300)
plot(plot_hist)
dev.off()

#square root stem.height data as cannot log to visual normal distribution
stem.sq<-sqrt(mushroom$stem.height)

#plot histogram of stem height for both edible and poisonous mushrooms
plot_hist2<-ggplot(data=mushroom, aes(x=stem.sq, fill=class.fac))+
  geom_histogram(bins=30, alpha=0.5)+
  labs(title='Stem Height Distribution in Edible and Poisonous Mushrooms',
    x='Mushroom Stem Height (cm)', y='Frequency')+
  theme_classic()+
  scale_fill_manual(values = c("grey","red"), name = "Mushroom Type")+
  theme(plot.title = element_text(size = 12, face = "bold"))
print(plot_hist2)
tiff("~/Desktop/UNI/MSc Bioinformatics/R. Studio/Assessment 2/Question 3 Stem Hist.tiff", height=10, width=20, units='cm', compression='lzw', res=300)
plot(plot_hist2)
dev.off()

#perform t-test for cap.diameter and class
t.test(cap.diameter~class, data = mushroom) 
#t = 43.359, df = 43335, p-value < 2.2e-16, this shows a high significant difference

#perform t-test for stem height and class
t.test(stem.height~class, data = mushroom)
#t = 29.84, df = 54422, p-value < 2.2e-16, this shows a high significant difference


