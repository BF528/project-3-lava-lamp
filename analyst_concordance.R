library(dplyr)
library(ggplot2)
#read data in 

#map for calculating observed overlap
map <- read.csv("/project/bf528/project_3/refseq_affy_map.csv")

#RNAseq data 
lef_rna <- read.csv("deseq_rnaseq/deseq_results_ahr.csv")
flu_rna <- read.csv("deseq_rnaseq/deseq_results_car.csv")
ifo_rna <- read.csv("deseq_rnaseq/deseq_results_dna.csv")

#microarray data 
lef_micro <- read.csv('LEFLUNOMIDE_limma_results.csv')
flu_micro <- read.csv('FLUCONAZOLE_limma_results.csv') 
ifo_micro <- read.csv('IFOSFAMIDE_limma_results.csv')

#Data------for each MOA 
#get deferentially expressed genes
total_DEGs_micro_l <- nrow(filter(lef_micro, adj.P.Val < 0.05 & abs(logFC) > log2(1.5)))
total_DEGs_RNA_l <- nrow(filter(lef_rna, padj < 0.05))
#subset-DEG
DEGs_micro_l <- as.data.frame(filter(lef_micro, adj.P.Val < 0.05 & abs(logFC) > log2(1.5)))
DEGs_RNA_l <-  as.data.frame(filter(lef_rna, padj < 0.05))

#get deferentially expressed genes, not duplicated, flu
total_DEGs_micro_f <- nrow(filter(flu_micro, adj.P.Val < 0.05 & abs(logFC) > log2(1.5)))
total_DEGs_RNA_f <- nrow(filter(flu_rna, padj < 0.05))
DEGs_micro_f <- as.data.frame(filter(flu_micro, adj.P.Val < 0.05 & abs(logFC) > log2(1.5)))
DEGs_RNA_f <-  as.data.frame(filter(flu_rna, padj < 0.05))

#get deferentially expressed genes, 
total_DEGs_micro_i <- nrow(filter(ifo_micro, adj.P.Val < 0.05 & abs(logFC) > log2(1.5)))
total_DEGs_RNA_i <- nrow(filter(ifo_rna, padj < 0.05))
DEGs_micro_i <- as.data.frame(filter(ifo_micro, adj.P.Val < 0.05 & abs(logFC) > log2(1.5)))
DEGs_RNA_i <-  as.data.frame(filter(ifo_rna, padj < 0.05))

#---------------------------COMPUTE CONCORDANCE-----------------------------

#observed overlap
n0_calculation <- function(overlap_mic, overlap_rn, map, above, below) {
  #merge each with map some will be duplicated values because they match to different IDs
  overlap_micro <- merge(x = map, y = overlap_mic, by.x = "PROBEID", by.y = 'X')
  overlap_rna <- merge(x = map, y = overlap_rn, by.x = "REFSEQ", by.y = 'X')
  overlap_df <- merge(x = overlap_micro, y =overlap_rna)
  med <- median(overlap_df$baseMean)
  
  #for calculating above and below mean, if 'n' for both skip these ifs. 
  if (above == 'y'){
    overlap_df <-filter(overlap_df, baseMean >= med)
  }
  #merge both dataframes 
  if (below == 'y'){
    overlap_df <-filter(overlap_df, baseMean < med)
  }
  
  #check if log fc direction is the same, and only keep rows that have same direction
  dir_mic <- overlap_df$logFC > 0
  dir_rna <- overlap_df$log2FoldChange > 0
  overlap_df$logfc_dir_micro <- ifelse(dir_mic,'Positive','Negative')  
  overlap_df$logfc_dir_rna <- ifelse(dir_rna,'Positive','Negative')  
  overlap_df <- overlap_df[overlap_df$logfc_dir_micro==overlap_df$logfc_dir_rna, ]
  
  #get total observed overlaps
  n0 <- length(overlap_df$SYMBOL)
  return (n0)
}

#observed overlaps stored in variables
n0_l<- n0_calculation(overlap_mic =DEGs_micro_l, overlap_rn =DEGs_RNA_l, map = map, above = 'n', below = 'n')
n0_f<- n0_calculation(overlap_mic =DEGs_micro_f, overlap_rn =DEGs_RNA_f, map = map, above = 'n', below = 'n')
n0_i<- n0_calculation(overlap_mic =DEGs_micro_i, overlap_rn =DEGs_RNA_i, map = map, above = 'n', below = 'n')


#compute concordance
concord <- function (n0, total_DEG_mi, total_DEG_rn){
  n1 <- total_DEG_mi
  n2 <- total_DEG_rn
  #N: total possible detected genes from microarray and rna-seq platforms 
  N <- 13079 
  #background corrected
  bg_intersect<- (N*n0 - n1*n2)/(n0+N-n1-n2)
  concordance <- (2*abs(bg_intersect))/(n1+n2)
  return (concordance)
}

#plot concordance against MICRO for each treatment------
png("concord.png")
x <- c(total_DEGs_micro_l, total_DEGs_micro_f, total_DEGs_micro_i)
y <- c(concord(n0_l, total_DEGs_micro_l, total_DEGs_RNA_l), 
       concord(n0_f, total_DEGs_micro_f, total_DEGs_RNA_f), 
       concord(n0_i, total_DEGs_micro_i, total_DEGs_RNA_i))

plot(x, y, xlab = 'Micro Array genes', ylab = 'Concordance')
text(x,y-0.005,labels=c('lef','flu', 'ifo'))
dev.off()

#plot concordance against RNA-seq for each treatment--------
png("concord1.png")
x <- c(total_DEGs_RNA_l, total_DEGs_RNA_f, total_DEGs_RNA_i)
y <- c(concord(n0_l, total_DEGs_micro_l, total_DEGs_RNA_l), 
       concord(n0_f, total_DEGs_micro_f, total_DEGs_RNA_f), 
       concord(n0_i, total_DEGs_micro_i, total_DEGs_RNA_i))
plot(x, y, xlab = 'RNA seq genes', ylab = '')
text(x,y-0.005,labels=c('lef','flu', 'ifo'))
dev.off()

#Calculating Above and Below--------------------------------
#above median and below leflu, n0 (observed overlaps)
a<- n0_calculation(overlap_mic =DEGs_micro_l, overlap_rn =DEGs_RNA_l, map = map, above = 'y', below = 'n')
b<- n0_calculation(overlap_mic =DEGs_micro_l, overlap_rn =DEGs_RNA_l, map = map, above = 'n', below = 'y')

#above median and below fluco, n0 (observed overlaps)
c<- n0_calculation(overlap_mic =DEGs_micro_f, overlap_rn =DEGs_RNA_f, map = map, above = 'y', below = 'n')
d<- n0_calculation(overlap_mic =DEGs_micro_f, overlap_rn =DEGs_RNA_f, map = map, above = 'n', below = 'y')

#join above, overlap, below into one variable for each treatment
concord_l = c(concord(a, total_DEGs_micro_l, total_DEGs_RNA_l), 
              concord(n0_l, total_DEGs_micro_l, total_DEGs_RNA_l), 
              concord(b, total_DEGs_micro_l, total_DEGs_RNA_l))

concord_f <-c(concord(c, total_DEGs_micro_f, total_DEGs_RNA_f), 
              concord(n0_f, total_DEGs_micro_f, total_DEGs_RNA_f), 
              concord(d, total_DEGs_micro_f, total_DEGs_RNA_f))

concord_i <-c(0,0,0)

#combine into vector and then dataframe for barplot:
concordances <- c(concord_l, concord_f,concord_i)
chemical <- rep(c('Leflunomide/AhR', 'Fluconazole/CAR', 'Ifosfamide/DNA'), each =3)
group <- rep(c('above', 'overall', 'below'), 3)

df <- data.frame(concordances, chemical, group)

#plot above, below, overall concordance for each treatment 
p <- ggplot(data=df, aes(x=chemical, y=concordances, fill=group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() + ggtitle('Above, Overall, and Below Median Subsets Concardances per Treatment')

#save in png
png('above_below.png')
p
dev.off()
