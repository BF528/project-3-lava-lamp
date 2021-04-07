#Clear environment
rm(list = ls())
gc()

setwd('/projectnb/bf528/users/lava_lamp/project_3/esaake_pr3')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install(c("RDAVIDWebService", "RGSEA", "SeqGSEA","GOstats", "AnnotationDbi"))

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("gplots")

library("RDAVIDWebService")

library("RGSEA")

library("GOstats")

library("dplyr")

library("ggplot2")
#library("ggarrange")
library("ggpubr")
library("gplots")

#clear environment


#BiocManager::available("GOstats")


deseqcount_ahr_path<-"/projectnb/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_norm_counts_ahr.csv"
deseqcount_car_path<-"/projectnb/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_norm_counts_car.csv"
deseqcount_dna_path<-"/projectnb/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_norm_counts_dna.csv"

DEgAhrrespath_RNAseq<-"/projectnb/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_results_ahr.csv"
DEgCarrespath_RNAseq<-"/projectnb/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_results_car.csv"
DEgDnarespath_RNAseq<-"/projectnb/bf528/users/lava_lamp/project_3/deseq_rnaseq/deseq_results_dna.csv"

DEgFlurespath_microarray<-"/projectnb/bf528/users/lava_lamp/project_3/FLUCONAZOLE_limma_results.csv"
DEgIforespath_microarray<-"/projectnb/bf528/users/lava_lamp/project_3/IFOSFAMIDE_limma_results.csv"
DEgLefrespath_microarray<-"/projectnb/bf528/users/lava_lamp/project_3/LEFLUNOMIDE_limma_results.csv"

datacount_ahr<-read.csv(deseqcount_ahr_path, row.names=1)
datacount_car<-read.csv(deseqcount_car_path, row.names=1)
datacount_dna<-read.csv(deseqcount_dna_path, row.names=1)


RNAseq_result_ahr<-read.csv(DEgAhrrespath_RNAseq, row.names=1)
RNAseq_result_car<-read.csv(DEgCarrespath_RNAseq, row.names=1)
RNAseq_result_dna<-read.csv(DEgDnarespath_RNAseq, row.names=1)

m_array_result_flu<-read.csv(DEgFlurespath_microarray, row.names=1)
m_array_result_ifo<-read.csv(DEgIforespath_microarray, row.names=1)
m_array_result_lef<-read.csv(DEgLefrespath_microarray, row.names=1)


#CLUSTERING
-----------
  
  #Setting seed for reproducibility
  set.seed(100)


#function to set color for each moa type
#---------------------------------------
  #takes a two arguments:
  #flag to know which moa type to generate specific color
  #data --moa matrix for either ahr,car or dna

setmoacolor<-function(moadata, flag){
  
  #creating an empty string for color
  color='';
  
  newcolnames<-c(names(moadata));
  
  if(flag=='ahr') {color='chartreuse'}
  if(flag=='car') {color='blue'}
  if(flag=='dna') {color='brown1'}
  
  #Generating color code for each moa type
  MoA_color<-rep(color,ncol(moadata))
  
  tmoadata<-t(moadata)
  ##transpose data and append color code 
  datacount<- mutate(as.data.frame(tmoadata), MOAcolor= MoA_color)
  
  #newset of headernames
  newrownames<-names(datacount)
  
  #Resetting the transposed coloumname
  datacount<-t(datacount);
  row.names(datacount)<-newrownames;
  colnames(datacount)<-c(names(moadata)); #new colnames set above
  
  return (datacount)
  
}

bio_datacount_ahr<-setmoacolor(datacount_ahr,'ahr')  
bio_datacount_car<-setmoacolor(datacount_car,'car') 
bio_datacount_dna<-setmoacolor(datacount_dna,'dna') 

#saving row names of data for later use. Rownames are the same for all there files so , saving one is enough
rownm<-row.names(bio_datacount_ahr)

#creating a list pool from all three norm_count_results
combined_norm_count_moa<-cbind.data.frame(bio_datacount_ahr, bio_datacount_car, bio_datacount_dna)

#After transformation, resetting the rownames with previously copied
row.names(combined_norm_count_moa)<-rownm

#getting a list of my colors representing each moa sample
colside_colors<-combined_norm_count_moa['MOAcolor',]
colside_colors<-t(colside_colors)
nrow(colside_colors)

#Coefficient of variation codes
  #Function execute_CoVtest
  #Input: Dataframe or matrix
  #Output: List of row-indices

#Transposing data to allow for clustering users
#removing the row with MOAcolor to avoid it from interfering with clustering
  #Gexpdata<- as.data.frame(t(combined_norm_count_moa[!row.names(combined_norm_count_moa)=='MOAcolor',]));
  bioanal_data<- as.data.frame(t(combined_norm_count_moa[-10696,]), row.names=1);
  

#testing for null values
any(is.na(bioanal_data))#To check for missing or null values
  
#function to clear unused memory space due to error
 gc()
 
#Scaling values
#unlisting resolves error with as.munerical error
data_scaled <- scale(data.matrix(bioanal_data))

#Finding distance between data points by using Euclidean
distanced_data<- dist(data_scaled, method="euclidean")
any(is.na(distanced_data))#checking if any is null

#computing hierarchical cluster average
cluster_average<-hclust(distanced_data, method='average')
plot(cluster_average)


#creating a cutree object to specify the number of clusters.
#because we are targeting AHR, CAR & DNA, I will limit with k=3 ie. 3 clusters
clust_trim<-cutree(cluster_average, k=3)

#storing the cluster return as dataframe
clust_df<-data.frame(clust_trim)


#Lets redraw the cluster and include three cut lines
  plot(cluster_average)
  rect.hclust(cluster_average, k=3, border=2:4)
  abline(h=NULL, col='red');#

#Draw the dendogram in different colors
  suppressPackageStartupMessages(library(dendextend)); #required package
  dendo_obj_avg<-as.dendrogram(cluster_average) #converting cluster average to dendogram object
  colored_dendo_of_dendo_obj_avg<-color_branches(dendo_obj_avg, k=3)

#opening an image object for saving
#jpeg('resultsf/clusterplot.jpg', width = 600, height = 600)
plot(colored_dendo_of_dendo_obj_avg)
#dev.off()


#Creating a cluster data to pass to heatmap by combining data with cutree object
bioanal_data_cluster<-cbind(bioanal_data, cluster=clust_trim)
table(bioanal_data_cluster$cluster)


#must be cluster #solution to Error: 'x' must be a numeric matrix #heatmap(as.matrix(dataset[, -1]))
#bioanal_data_cluster[,-1]

#setting color for colsidecolors #ColSideColors only takes vectors
rowheaders<-row.names(bioanal_data);
columnames<-names(bioanal_data);


#saving heatmap to jpeg
jpeg('results/heatmapnew5.jpg',  width = 900, height = 600 )
  #making sure colside color matches data for heatmap
    #colside_colors<-colside_colors[row.names(as.data.frame(colside_colors)) %in% c(rownames(clust_df)), ]
  #heatmap(data.matrix(data.frame(t(bioanal_data_cluster))), ColSideColors=colside_colors, scale="column", main="Gene Expression Heatmap",ylab= "Genes", xlab="Samples")
  heatmap(data.matrix(data.frame(bioanal_data_cluster)), RowSideColors=colside_colors, scale="column", main="Gene Expression Heatmap",ylab= "Samples", xlab="Genes")
  #legend for rowside colors
  legend("topleft", title = "MOA side color",legend=c("AHR","CAR", "DNA"), 
         fill=c("chartreuse","blue","brown1"), cex=0.8, box.lty=0)

  
dev.off()

#Coefficiet of variation filtering
execute_CoVtest<-function(givendata){
  #Variable Initialization
  numofrows<-nrow(givendata)
  numofcolumns<-ncol(givendata)
  filterpassers<-c();
  count<-0;
  list_index<-1;
  CoVthreshold<-0.186
  
  #loop to traverse rows
  for (i in 1:numofrows){
    gexpmean<-mean(as.numeric(givendata[i,]))
    std_dev<-sd(givendata[i,])
    CoV<-std_dev/gexpmean
    #test for condition of cov>CoVthreshold
    if (CoV>CoVthreshold){
      filterpassers[list_index]<-i
      list_index=list_index+1
    }
  }
  return(filterpassers)
}

#removing the color row to avoid interference in analysis
lastitem<-nrow(combined_norm_count_moa)

filtered_indices<-execute_CoVtest(combined_norm_count_moa[-lastitem,])
filtered_biodata<-combined_norm_count_moa[filtered_indices,]

#trying with filtered data
#filtering resulted in 10639 number of genes which is very close to before filtering 10695
#filtering won't be used because clustering already done well
bioanal_data_filt<- as.data.frame(t(filtered_biodata), row.names=1);


#Filtering for significant genes
#BiocManager::install("")

#The list provided is unfiltered and unordered
  #filtering for only Genebank accession(Genes) with p-adust < 0.05 and orderind in descending
    RNAseq_result_ahr <- RNAseq_result_ahr %>% filter(padj < 0.05) %>%  arrange(desc(padj))
    RNAseq_result_car <- RNAseq_result_car %>% filter(padj < 0.05) %>%  arrange(desc(padj))
    RNAseq_result_dna <- RNAseq_result_dna %>% filter(padj < 0.05) %>%  arrange(desc(padj))

    #writing the qualifying genes/Genebank accession to csv to be used on DAVID website
    #https://david.ncifcrf.gov/
    
     #write.csv(row.names(RNAseq_result_ahr), 'results/ahr_list_for_pathway.csv')
     #write.csv(row.names(RNAseq_result_car),'results/car_list_for_pathway.csv')
     #write.csv(row.names(RNAseq_result_dna), 'results/dna_list_for_pathway.csv')

#USING LIST FROM MICROARRAY
    #The list provided is unfiltered and unordered
    #filtering for only Affymetrix probe_IDs with p-adust < 0.05 and orderind in descending
      m_array_result_flu <- m_array_result_flu %>% filter(adj.P.Val < 0.05) %>%  arrange(desc(adj.P.Val))
      m_array_result_ifo <- m_array_result_ifo %>% filter(adj.P.Val < 0.05) %>%  arrange(desc(adj.P.Val))
      m_array_result_lef <-m_array_result_lef %>% filter(adj.P.Val < 0.05) %>%  arrange(desc(adj.P.Val))
     
     
     #commented out to avoid overwriting by mistake
     #writing the qualifying genes/Genebank accession to csv to be used on DAVID website
      #https://david.ncifcrf.gov/
      #write.csv(row.names(m_array_result_flu), 'results/m_flu_car_list_forpathway.csv')
      #write.csv(row.names(m_array_result_ifo),'results/m_ifo_dna_list_forpathway.csv')
      #write.csv(row.names(m_array_result_lef), 'results/m_lef_Ahr_list_forpathway.csv')

