# Load Packages
library(tidyverse)
library(data.table)
library(DESeq2)

#Sample Information
sample_info <- read_csv("/project/bf528/project_3/toxgroups/toxgroup_3_rna_info.csv")
sample_info_control <- sample_info %>% filter(mode_of_action == "Control")

# Read Count Tables
SRR1177981 <- fread("featureCounts/SRR1177981Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177982 <- fread("featureCounts/SRR1177982Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1177983 <- fread("featureCounts/SRR1177983Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178008 <- fread("featureCounts/SRR1178008Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178009 <- fread("featureCounts/SRR1178009Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178010 <- fread("featureCounts/SRR1178010Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178014 <- fread("featureCounts/SRR1178014Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178021 <- fread("featureCounts/SRR1178021Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)
SRR1178047 <- fread("featureCounts/SRR1178047Aligned.sortedByCoord.out.bam.txt", skip = 1, header = T)

# Only using the "Geneid" and count columns (columns 1 and 7 for each of the count files)
SRR1177981 <- SRR1177981[,c(1,7)] %>% dplyr::rename(SRR1177981 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1177981Aligned.sortedByCoord.out.bam")
SRR1177982 <- SRR1177982[,c(1,7)] %>% dplyr::rename(SRR1177982 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1177982Aligned.sortedByCoord.out.bam")
SRR1177983 <- SRR1177983[,c(1,7)] %>% dplyr::rename(SRR1177983 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1177983Aligned.sortedByCoord.out.bam")
SRR1178008 <- SRR1178008[,c(1,7)] %>% dplyr::rename(SRR1178008 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1178008Aligned.sortedByCoord.out.bam")
SRR1178009 <- SRR1178009[,c(1,7)] %>% dplyr::rename(SRR1178009 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1178009Aligned.sortedByCoord.out.bam")
SRR1178010 <- SRR1178010[,c(1,7)] %>% dplyr::rename(SRR1178010 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1178010Aligned.sortedByCoord.out.bam")
SRR1178014 <- SRR1178014[,c(1,7)] %>% dplyr::rename(SRR1178014 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1178014Aligned.sortedByCoord.out.bam")
SRR1178021 <- SRR1178021[,c(1,7)] %>% dplyr::rename(SRR1178021 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1178021Aligned.sortedByCoord.out.bam")
SRR1178047 <- SRR1178047[,c(1,7)] %>% dplyr::rename(SRR1178047 = "/projectnb/bf528/users/lava_lamp/project_3/bams/SRR1178047Aligned.sortedByCoord.out.bam")

# Merge Count Tables
merged <- Reduce(function(...) merge(..., by = "Geneid", all=TRUE), list(SRR1177981, SRR1177982, SRR1177983, SRR1178008, SRR1178009, SRR1178010, SRR1178014, SRR1178021, SRR1178047))
rm(SRR1177981, SRR1177982, SRR1177983, SRR1178008, SRR1178009, SRR1178010, SRR1178014, SRR1178021, SRR1178047)

# Boxplot of distribution
merged_bp <- merged %>% pivot_longer(!Geneid, names_to = "Sample", values_to = "Count")
png(filename = "figures/bp_count_distribution_zeros_limits.png")
merged_bp %>% ggplot(mapping = aes(x = Sample, y = Count)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Boxplot of Sample Count Distributions") +
  ylim(0, 1000)
dev.off()

# Control
control_counts <- read.csv("/project/bf528/project_3/samples/control_counts.csv") %>%
  select(Geneid, as.character(sample_info_control[[1]]))

# Merge controls with count matrices
merged_controls <- merge(merged, control_counts, by = "Geneid") %>% column_to_rownames(var="Geneid")

##### DESeq #####

# Filter out rows that have all zeros
merged_filtered <- subset(merged_controls, rowSums(merged_controls == 0) == 0)

# Boxplot of count distributions, with zeros removed
png(filename = "figures/bp_count_distribution_nozeros_limits.png")
merged_filtered %>% select(c(1:9)) %>% pivot_longer(c(1:9), names_to = "Sample", values_to = "Count") %>%
  ggplot(mapping = aes(x = Sample, y = Count)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Boxplot of Sample Count Distributions, Zero Counts Removed") +
  ylim(0, 1000)
dev.off()
  
# Subsetting runs to be put through DESeq individually
ahr_leflunomide <- sample_info %>% filter(mode_of_action == "AhR" | mode_of_action == "Control" & vehicle == "CORN_OIL_100_%") %>% select(Run)
car_fluconazole <- sample_info %>% filter(mode_of_action == "CAR/PXR" | mode_of_action == "Control" & vehicle == "CORN_OIL_100_%") %>% select(Run)
dna_ifosfamide <- sample_info %>% filter(mode_of_action == "DNA_Damage" | mode_of_action == "Control" & vehicle == "SALINE_100_%") %>% select(Run)

# Create the DESeq objects
#AhR - Leflunomide
dds_ahr <- DESeqDataSetFromMatrix(
  countData = merged_filtered %>% select(ahr_leflunomide[[1]]),
  colData = sample_info[c(1:3, 10:12),],
  design = ~ mode_of_action
)
#CAR/PXR - Fluconazole
dds_car <- DESeqDataSetFromMatrix(
  countData = merged_filtered %>% select(car_fluconazole[[1]]),
  colData = sample_info[c(4:6, 10:12),],
  design = ~ mode_of_action
)
#DNA Damage - Isosfamide
dds_dna <- DESeqDataSetFromMatrix(
  countData = merged_filtered %>% select(dna_ifosfamide[[1]]),
  colData = sample_info[c(7:9, 13:15),],
  design = ~ mode_of_action
)

# Relevel mode_of_action as factor for each DESeq object
dds_ahr$mode_of_action <- relevel(dds_ahr$mode_of_action, ref='Control')
dds_car$mode_of_action <- relevel(dds_car$mode_of_action, ref='Control')
dds_dna$mode_of_action <- relevel(dds_dna$mode_of_action, ref='Control')

# Run DESeq
dds_ahr <- DESeq(dds_ahr)
dds_car <- DESeq(dds_car)
dds_dna <- DESeq(dds_dna)

# DESeq Results
res_ahr <- results(dds_ahr, contrast=c('mode_of_action','AhR','Control'))
res_ahr <- lfcShrink(dds_ahr, coef=2)
res_car <- results(dds_car, contrast=c('mode_of_action','CAR/PXR','Control'))
res_car <- lfcShrink(dds_car, coef=2)
res_dna <- results(dds_dna, contrast=c('mode_of_action','DNA_Damage','Control'))
res_dna <- lfcShrink(dds_dna, coef=2)

# Write DE results
write.csv(res_ahr,'deseq_results_ahr.csv')
write.csv(res_car,'deseq_results_car.csv')
write.csv(res_dna,'deseq_results_dna.csv')

# Write Matrix of Normalized Counts
write.csv(counts(dds_ahr,normalized=TRUE),'deseq_norm_counts_ahr.csv')
write.csv(counts(dds_car,normalized=TRUE),'deseq_norm_counts_car.csv')
write.csv(counts(dds_dna,normalized=TRUE),'deseq_norm_counts_dna.csv')

# Further analysis with DESeq data 
counts_ahr <- read.csv("/projectnb2/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_results_ahr.csv", header = TRUE, row.names = 1)
counts_car <- read.csv("/projectnb2/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_results_car.csv", header = TRUE, row.names = 1) 
counts_dna <- read.csv("/projectnb2/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_results_dna.csv", header = TRUE, row.names = 1) 

# Number of genes with padj < 0.05
counts_ahr %>% filter(padj < 0.05) %>% summarize(count = n()) #1389
counts_car %>% filter(padj < 0.05) %>% summarize(count = n()) #3499
counts_dna %>% filter(padj < 0.05) %>% summarize(count = n()) #91

# Top 10 DE Genes from each analysis
counts_ahr %>% top_n(-10, padj) %>% write.csv("topDE/topDE_ahr.csv", row.names = TRUE)
counts_car %>% top_n(-10, padj) %>% write.csv("topDE/topDE_car.csv", row.names = TRUE)
counts_dna %>% top_n(-10, padj) %>% write.csv("topDE/topDE_dna.csv", row.names = TRUE)

# Histograms of FC Values
png("figures/fc_ahr.png")
counts_ahr %>% filter(padj < 0.05) %>%
  ggplot(mapping = aes(x = log2FoldChange)) +
  geom_histogram() +
  labs(title = "Fold Change Histogram, AhR", x = "Log2 Fold Change", y = "Frequency")
dev.off()
png("figures/fc_car.png")
counts_car %>% filter(padj < 0.05) %>%
  ggplot(mapping = aes(x = log2FoldChange)) +
  geom_histogram() +
  labs(title = "Fold Change Histogram, CAR/PTX", x = "Log2 Fold Change", y = "Frequency")
dev.off()
png("figures/fc_dna.png")
counts_dna %>% filter(padj < 0.05) %>%
  ggplot(mapping = aes(x = log2FoldChange)) +
  geom_histogram() +
  labs(title = "Fold Change Histogram, DNA Damage", x = "Log2 Fold Change", y = "Frequency")
dev.off()

# Volcano Plots
library(EnhancedVolcano)

EnhancedVolcano(res_ahr,
                lab = rownames(res_ahr),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'AhR Differentially Expressed Genes',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                legendPosition = 'bottom')

counts_ahr %>%
  ggplot(mapping = aes(x = log2FoldChange, y = -log(pvalue))) +
  geom_point()