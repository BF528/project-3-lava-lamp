library(limma)
library(dplyr)
library(gt)
library(gridExtra)
library(ggplot2)
# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/toxgroups/toxgroup_3_mic_info.csv',as.is=TRUE)

rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt', as.is=TRUE, header=TRUE, sep = "\t", row.names=1)

#different comparisons
leflu <- filter(samples, chemical == 'LEFLUNOMIDE' | chemical == 'Control' & vehicle == 'CORN_OIL_100_%') 
flu <- filter(samples, chemical == 'FLUCONAZOLE'| chemical == 'Control' & vehicle == 'CORN_OIL_100_%')
Ifo <- filter(samples, chemical == 'IFOSFAMIDE'| chemical == 'Control' & vehicle == 'SALINE_100_%')

samps <- list(leflu,flu, Ifo)

for (i in samps) {
# subset the full expression matrix to just those in this comparison
  rma.subset <- rma[paste0('X',i$array_id)]
  name <- i$chemical[1]
  # construct a design matrix modeling treatment vs control for use by limma
  design <- model.matrix(
    ~factor(
      i$chemical,
      levels=c('Control', name)
    )
  )
  colnames(design) <- c('Intercept',name)
  
  # run limma
  fit <- lmFit(rma.subset, design)
  fit <- eBayes(fit)
  t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')
  
  # write out the results to file
  write.csv(t,paste0(name,'_limma_results.csv')) 
}

l<- read.csv('LEFLUNOMIDE_limma_results.csv')
f<- read.csv('FLUCONAZOLE_limma_results.csv')
i<- read.csv('IFOSFAMIDE_limma_results.csv')


#number of significant genes, save into variables as well
nrow(filter(l, adj.P.Val < 0.05 & abs(logFC) > log2(1.5))) #183
nrow(filter(f, adj.P.Val < 0.05 & abs(logFC) > log2(1.5))) #726
nrow(filter(i, adj.P.Val < 0.05 & abs(logFC) > log2(1.5))) #0

l_sig <- filter(l, adj.P.Val < 0.05 & abs(logFC) > log2(1.5)) 
f_sig <-filter(f, adj.P.Val < 0.05 & abs(logFC) > log2(1.5)) 
i_sig <-filter(i, adj.P.Val < 0.05 & abs(logFC) > log2(1.5))

map <- read.csv("/project/bf528/project_3/refseq_affy_map.csv")
#top 10 into tables

l_10 <- merge(map, l_sig, by.x = 'PROBEID', by.y = 'X') %>% 
  select(PROBEID, SYMBOL, logFC, t,P.Value,adj.P.Val)  %>%
  arrange(adj.P.Val)%>%
  slice(1:10) %>% gt() %>%
  tab_header(title = md("This is the top 10 genes adj.p-val <0.05 LEFLUNOMIDE")) %>%
  tab_options(heading.background.color = "#EFFBFC",
              table_body.hlines.color = "#989898",
              table_body.border.top.color = "#989898")

l_10 %>%
  gtsave("leflu.png", expand = 10)

f_10 <- merge(map, f_sig, by.x = 'PROBEID', by.y = 'X') %>% 
  select(PROBEID, SYMBOL, logFC, t,P.Value,adj.P.Val)  %>%
  arrange(adj.P.Val)%>%
  slice(1:10) %>% gt() %>% 
  tab_header(title = md("This is the top 10 genes adj.p-val <0.05 FLUCONAZOLE"))  %>% 
  tab_options(heading.background.color = "#EFFBFC",
              table_body.hlines.color = "#989898",
              table_body.border.top.color = "#989898")

f_10 %>%
  gtsave("fluco.png", expand = 10)

#NO significant genes for Ifosfamide

#Histograms of fold change values

png("fc_histograms.png")
par(mfrow =c(1,2))
a <-hist(l_sig$logFC, main = 'LEFLUNOMIDE', xlab = "fold change", ylab = "Frequency", breaks = 20)
b <-hist(f_sig$logFC, main = 'FLUCONAZOLE', xlab = "fold change", ylab = "", breaks = 20)
c <- hist(i_sig$logFC, main = 'IFOSFAMIDE', xlab = "", ylab = "", breaks = 20) # NO values
dev.off() 
#scatter plots
plot1 <- ggplot(l, aes(x=logFC, y= -log10(P.Value)))+
  geom_point(size=2, shape=23) + xlab("")+ggtitle('LEFLUNOMIDE') +
  geom_vline(xintercept=log2(1.5), color ='blue')+
  geom_vline(xintercept=-log2(1.5), color ='blue')
plot2 <- ggplot(f, aes(x=logFC, y= -log10(P.Value))) + 
  geom_point(size=2, shape=23) + ylab("")+ xlab("") + ggtitle('FLUCONAZOLE') +
  geom_vline(xintercept=log2(1.5), color ='blue')+
  geom_vline(xintercept=-log2(1.5), color ='blue')

#plot3 <- ggplot(i_sig, aes(x=logFC, y= -log10(P.Value))) +
  #geom_point(size=2, shape=23) + xlab("") + ylab("")
png("fc_scatter.png")
grid.arrange(plot1, plot2, ncol=2, bottom="logFC")
dev.off() 
