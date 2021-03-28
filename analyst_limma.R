library(limma)
library(dplyr)
library(gt)
library(gridExtra)
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
nrow(filter(l, adj.P.Val < 0.05)) #466
nrow(filter(f, adj.P.Val < 0.05)) #1997
nrow(filter(i, adj.P.Val < 0.05)) #0

l_sig <- filter(l, adj.P.Val < 0.05) 
f_sig <-filter(f, adj.P.Val < 0.05) 
i_sig <-filter(i, adj.P.Val < 0.05)

#top 10 into tables
l_10 <- l %>% slice(1:10)%>% gt() %>% 
  tab_header(title = md("This is the top 10 genes adj.p-val <0.05 LEFLUNOMIDE")) %>%
  tab_options(heading.background.color = "#EFFBFC",
              table_body.hlines.color = "#989898",
              table_body.border.top.color = "#989898")

f_10 <- f %>% slice(1:10) %>% gt() %>% 
  tab_header(title = md("This is the top 10 genes adj.p-val <0.05 FLUCONAZOLE"))  %>% 
  tab_options(heading.background.color = "#EFFBFC",
              table_body.hlines.color = "#989898",
              table_body.border.top.color = "#989898")

i_10 <- i %>% slice(1:10) %>% gt() %>% 
  tab_header(title = md("This is the top 10 genes adj.p-val <0.05 IFOSFAMIDE")) %>% 
  tab_options(heading.background.color = "#EFFBFC",
              table_body.hlines.color = "#989898",
              table_body.border.top.color = "#989898")

#Histograms of fold change values
a <-hist(l$logFC, main = 'LEFLUNOMIDE', xlab = "", ylab = "Frequency", breaks = 20)
b <-hist(f$logFC, main = 'FLUCONAZOLE', xlab = "fold change", ylab = "", breaks = 20)
c <- hist(i$logFC, main = 'IFOSFAMIDE', xlab = "", ylab = "", breaks = 20)

#scatter plots
par(mfrow =c(3,1))
plot1 <- ggplot(l_sig, aes(x=logFC, y= -log10(P.Value)))+
  geom_point(size=2, shape=23) + xlab("")
plot2 <- ggplot(f_sig, aes(x=logFC, y= -log10(P.Value))) +
  geom_point(size=2, shape=23) + ylab("")
plot3 <- ggplot(i_sig, aes(x=logFC, y= -log10(P.Value))) +
  geom_point(size=2, shape=23) + xlab("") + ylab("")
grid.arrange(plot1, plot2, plot3, ncol=3)
