#libraries
library(tidyverse)
library(dplyr)
library(ggplot2)

#open dataframe
df <- read.csv('reshuffled_filtered_mod_peptides.csv', sep=',', header = TRUE)
head(df)
summary(df)

####Transform Data####
#address the issue of having '0' values before performing the natural log transformation
#replace those '0' values with '-1' 
#list of column names with 'ion'
ion_columns <- colnames(df)[grep("_ion\\d+", colnames(df))]
#replace '0' with '-1' in 'ion' columns
df[ion_columns][df[ion_columns] == 0] <- 0.1
#check
summary(df)

#take the natural log of non-negative values in 'ion' columns
df[, ion_columns] <- lapply(df[, ion_columns], function(x) ifelse(x > 0, log(x), x))

#create a boxplot for the logged ion intensities
boxplot(df[, ion_columns], ylab = "Log Transformed Ion Intensity", names = ion_columns)

####Normalize Data####
#calculate the median of the first log '_ion' column
first_ion_median <- median(df[, ion_columns[1]])
#calculate differences between medians
median_differences <- sapply(df[, ion_columns], function(x) median(x) - first_ion_median)
#create normalized versions of every column
normalized_df <- df
for (col in ion_columns) {
  # Check if the column contains numeric data before performing subtraction
  if (is.numeric(df[[col]])) {
    normalized_df[[col]] <- df[[col]] - median_differences[col]
  }}

#create a new boxplot to confirm that we have removed systematic differences between each channel
ion_columns <- colnames(normalized_df)[grep("_ion\\d+", colnames(normalized_df))]
#'_ion\\d+' looks for names with _ion and a number following it.
#boxplot for the normalized ion intensities
boxplot(normalized_df[, ion_columns], ylab = "Log Transformed Ion Intensity", names = ion_columns)


####Statistical Testing####
#perform t-test for every phosphopeptide
#filter for phosphopeptides
phospho_df <- normalized_df[grep("Phospho", normalized_df$Modifications), ]

#specify columns for comparison
ion_columns_S1 <- c("S1_36_ion2", "s1_36_ion10", "s1_36_ion11")
ion_columns_S2 <- c("s2_36_ion4", "s2_36_ion6", "s2_36_ion9")
ion_columns_mock <- c("mock_36_ion3", "mock_36_ion8", "mock_36_ion13")

# Create a function to perform t-tests
perform_t_tests <- function(values_S, values_mock) {
  p_value_S_vs_mock <- t.test(values_S, values_mock)$p.value
  return(p_value_S_vs_mock)
}

#apply t-tests using rowwise and mutate
ttest_results <- phospho_df %>%
  rowwise() %>%
  mutate(
    p_value_S1_vs_mock = perform_t_tests(across(all_of(ion_columns_S1)), across(all_of(ion_columns_mock))),
    p_value_S2_vs_mock = perform_t_tests(across(all_of(ion_columns_S2)), across(all_of(ion_columns_mock)))
  ) %>%
  ungroup()  #remove grouping for subsequent operations

#filter for significant phosphopeptides (p < 0.05)
significant_phospho <- ttest_results %>%
  filter(p_value_S1_vs_mock < 0.05 | p_value_S2_vs_mock < 0.05)


#### visualising differential expression with volcano plot ####
#create a new column for log2 fold change for S1 vs Mock
ttest_results$log2FC_S1_vs_mock <- rowMeans(ttest_results[, c("S1_36_ion2", "s1_36_ion11", "s1_36_ion10")]) - 
  rowMeans(ttest_results[, c("mock_36_ion3", "mock_36_ion13", "mock_36_ion8")])

#create a new column for log2 fold change for S2 vs Mock
ttest_results$log2FC_S2_vs_mock <- rowMeans(ttest_results[, c("s2_36_ion4", "s2_36_ion6", "s2_36_ion9")]) - 
  rowMeans(ttest_results[, c("mock_36_ion3", "mock_36_ion13", "mock_36_ion8")])

#create a volcano plot for S1 vs Mock
ggplot(ttest_results, aes(x = log2FC_S1_vs_mock, y = -log10(p_value_S1_vs_mock))) +
  geom_point(aes(color = factor(ifelse(p_value_S1_vs_mock < 0.05, "Significant", "Not Significant"))), size = 1.5) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (FC)", y = "-log10 (p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme(legend.position = "top") +  # position the legend
  guides(color = guide_legend(title = NULL)) +  # remove legend title
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "ivory4")) +  # set custom colors
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5))

#save plot
ggsave("S1_Vs_Mock_Volcano.png", dpi = 300, limitsize = FALSE)

#create a volcano plot for S2 vs Mock
ggplot(ttest_results, aes(x = log2FC_S2_vs_mock, y = -log10(p_value_S2_vs_mock))) +
  geom_point(aes(color = factor(ifelse(p_value_S2_vs_mock < 0.05, "Significant", "Not Significant"))), size = 1.5) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (FC)", y = "-log10 (p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme(legend.position = "top") +  # position the legend
  guides(color = guide_legend(title = NULL)) +  # remove legend title
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "ivory4")) +  # set custom colors
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5))

#save plot
ggsave("S2_Vs_Mock_Volcano.png", dpi = 300, limitsize = FALSE)

##maybe look at volcano plot by removing outliers?
# Filter out outliers based on log2 fold change threshold
S1_filtered_results <- subset(ttest_results, log2FC_S1_vs_mock > -2 & log2FC_S1_vs_mock < 2)
S2_filtered_results <- subset(ttest_results, log2FC_S2_vs_mock > -2 & log2FC_S2_vs_mock < 2)

#create a scatter plot for S1 vs Mock with filtered data
ggplot(S1_filtered_results, aes(x = log2FC_S1_vs_mock, y = -log10(p_value_S1_vs_mock))) +
  geom_point(aes(color = factor(ifelse(p_value_S1_vs_mock < 0.05, "Significant", "Not Significant"))), size = 1.5) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (FC)", y = "-log10 (p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme(legend.position = "top") +  # position the legend
  guides(color = guide_legend(title = NULL)) +  # remove legend title
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "ivory4")) +  # set custom colors
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5))
#save plot
ggsave("Filtered_S1_Vs_Mock_Volcano.png", dpi = 300, limitsize = FALSE)

#create a scatter plot for S2 vs Mock with filtered data
ggplot(S2_filtered_results, aes(x = log2FC_S2_vs_mock, y = -log10(p_value_S2_vs_mock))) +
  geom_point(aes(color = factor(ifelse(p_value_S2_vs_mock < 0.05, "Significant", "Not Significant"))), size = 1.5) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (FC)", y = "-log10 (p-value)") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme(legend.position = "top") + 
  guides(color = guide_legend(title = NULL)) + 
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "ivory4")) + 
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5))

#save plot
ggsave("Filtered_S2_Vs_Mock_Volcano.png", dpi = 300, limitsize = FALSE)
####Identify differentially expressed proteins####
# Split the Proteins column based on ";"
significant_phospho <- significant_phospho %>%
  separate_rows(Proteins, sep = ";", convert = TRUE)

#remove duplicate rows based on the Proteins column
unique_proteins <- significant_phospho %>%
  distinct(Proteins, .keep_all = TRUE)

##same for the phospho df
new_phospho_df <- phospho_df %>%
  separate_rows(Proteins, sep = ";", convert = TRUE)

new_phospho_df <- new_phospho_df %>%
  distinct(Proteins, .keep_all = TRUE)

#save the 'unique_proteins' and 'phospho' df for further analysis
write.csv(new_phospho_df, "phospho_df.csv", row.names = FALSE)
write.csv(unique_proteins, "unique_proteins.csv", row.names = FALSE)

##figure editting ##


