# Concordance of Microarray and RNA-Seq Differential Gene Expression
Both microarray and RNA-Sequencing measure genome-wide abundance of mRNA molecules. However, the concordance between the results obtained from these two technologies have rarely been compared directly, across multiple treatment groups, and while controlling for a variety of biological factors. Wang et al present a large study comparing microarray and RNA-Seq gene expression data from a set of toxicological treatments with known mechanism of action (MOA) measured in rat liver. The goal of the study was to characterize the concordance of differential gene expression across platforms, test and compare how effective each platform is at detecting expected pathway-level effects based on each treatment’s MOA, and assess the MOA prediction accuracy of each platform using a test set.

The goal of our project is to characterize the concordance of differential gene expression across platforms, test and compare how effective each platform is at detecting expected pathway-level effects based on each treatment’s MOA, and assess the MOA prediction accuracy of each platform using a test set. Therefore, we will reproduce the results from Figure 2A and Figures 3B and C, as well as compare the pathway enrichment results reported in the paper.

__Original Analysis:__ 
+ Wang, Charles, Binsheng Gong, Pierre R. Bushel, Jean Thierry-Mieg, Danielle Thierry-Mieg, Joshua Xu, Hong Fang, et al. 2014. “A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data” Nature Biotechnology 32 (9): 926–32. PMID: 4243706

__Goals of This Project:__
1. Align short reads to the rat genome using STAR and quantify expression using a read counting strategy
2. Perform differential expression analysis of RNA-Seq data with DESeq2
3. Perform differential expression analysis of pre-normalized microarray expression data using limma
4. Map between Affymetrix and refSeq identifier systems

### Contributors

+ Daisy Wenyan Han daisyhan@bu.edu
+ Divya Sundaresan divyas3@bu.edu
+ Alec Jacobsen aggjacob@bu.edu
+ Emmanuel Saake esaake@bu.edu

### Repository Contents and Suggested Workflow
1. `p3_install_packages.sh` and `snakefile` - Run through snakemake, these files make use of the STAR and MultiQC tools to align the pre-processed FASTQ files used in the original analysis to a reference genome. Quality control is performed on these samples to ensure an accurate alignment. This generates nine BAM files, to be used in further analysis. MultiQC is performed on these files to ensure adequate quality.
2. `featureCounts.qsub`, `multiqc_report.html` and `programmer.R` - Count matrices are generated from the nine BAM files generated above using the featureCounts algorithm, available in the `subread` package. Again, MultiQC is run on the summary files to best visualize the quality metrics of these matrices, with a full report documented in `multiqc_report.html`. Differential gene expression analysis is then carried out using DESeq2 in R.
3. `analyst_limma.R` and `analyst_concordance.R` - Differential gene expression analysis is performed on the same treatment samples as above, analyzed using microarray technology rather than RNA-Sequencing. Limma is used to perform this analysis. Concordance between these two analyses is then computed. 
4. 'Bio_analysis.R' - Filter out significant DE genes for gene set enrichment analysis with DAVID (https://david.ncifcrf.gov/). Generate a clustered heatmap with the MOA based  normalized count  matrices from DESeq2 object. 
