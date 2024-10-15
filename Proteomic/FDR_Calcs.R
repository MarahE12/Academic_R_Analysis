#load libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tibble)

#read data in and ignore first row
df = read.csv("Animal_20200325_TM_HStmtpro_CoV12_ph2_fr1.txt",sep="\t", skip = 1)

#order df by e value
df = df[order(df$e.value),]

#### Identifying false positives and targets ####

##add additional columns to df identifying phosphorylated peptides and
##to count number of FP, Targets and TP to allow FDR calculation and later plotting

#create column to identify phosphorylated peptides
df$PTM <- grepl("\\[79.966", df$modified_peptide)

#create a new column for decoy status in df to identify FP and Target rows
df$Decoy_status <- grepl("DECOY", df$protein)
#create 'FP_count' and 'Target_count' columns and start count at 0
df$FP_count <- cumsum(df$Decoy_status)
df$Target_count <- cumsum(!df$Decoy_status)

#calculate cumulative FDR for each row in df
df$FDR <- (df$FP_count/(df$FP_count + df$Target_count))
df$TP <- (df$Target_count-df$FP_count)

#filter df to only show FP
fp_df <- df[df$Decoy_status == "TRUE",]
#filter df to only show Targets
target_df <- df[df$Decoy_status == "FALSE", ]


#### Identifying significant e_values & calculating FDR ####
#set e value threshold
e_value <- 0.05

#identify how many FPP flag as significant according to e value threshold
significant_fp <- subset(fp_df, `e.value` < e_value)

#identify significant Targets
significant_target <- subset(target_df, `e.value` < e_value)

#calculate FDR of the significant data
significant_fdr = nrow(significant_fp)/(nrow(significant_fp) + nrow(significant_target))

#### Visualising FDR vs TP ####
ggplot(significant_target, aes(x = TP, y = FDR)) +
  geom_point() +
  xlab("True Positive (TP)") +
  ylab("False discovery rate (FDR)") +
  xlim(0, 7800) +
  theme(axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))

#save plot
ggsave("FDR_vs_TP.png", dpi = 300, limitsize = FALSE)

#### Creating a FDR & PSMs count table####
##generate table to show how many PSMs have been identified at 1%, 5%, and 10% FDR
#define the FDR thresholds
fdr_thresholds <- c(0.01, 0.05, 0.10)

#create an empty tibble df for the FDR table
fdr_table <- tibble(`FDR(%)` = character(), PSMs = integer())

#loop through each FDR threshold to find the corresponding TP
for (threshold in fdr_thresholds) {
  #find the max TP value before the FDR exceeds the threshold
  max_tp_at_fdr <- max(df$TP[df$FDR <= threshold])
  
  #append the results to the fdr_table
  fdr_table <- bind_rows(fdr_table, tibble(`FDR(%)` = paste0((threshold * 100), "%"), 
                                           PSMs = max_tp_at_fdr))}

#print the FDR table
print(paste("The FDR Table:\n", fdr_table))

#save the FDR table as a CSV file
write.csv(fdr_table, "FDR_Table.csv", row.names = FALSE, quote = FALSE)


#### Filtering for significant PTMs ####
#identify significant phosphopeptides starting serine residue
ptm_df <- subset(significant_target, PTM == "TRUE")

#save final main df results as csv (this will be submitted)
write.csv(df, "main_dataframe.csv", row.names = FALSE)
#save the df of significant phosphopeptides for further investigations of PTMs
write.csv(ptm_df, "significant_PTM.csv", row.names = FALSE)





