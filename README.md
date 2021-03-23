# Concordance of microarray and RNA-Seq differential gene expression

The goal of the study was to characterize the concordance of differential gene expression across platforms, test and compare how effective each platform is at detecting expected pathway-level effects based on each treatment’s MOA, and assess the MOA prediction accuracy of each platform using a test set. In this project, we will reproduce the results from Figure 2a and Figure 3b+c, as well as compare the pathway enrichment results reported in the paper.

Original Analysis to be reproduced: 
+ Wang, Charles, Binsheng Gong, Pierre R. Bushel, Jean Thierry-Mieg, Danielle Thierry-Mieg, Joshua Xu, Hong Fang, et al. 2014. “A comprehensive study design reveals treatment- and transcript abundance–dependent concordance between RNA-seq and microarray data” Nature Biotechnology 32 (9): 926–32. PMID: 4243706

Goals of this project:
1. Align short reads to the rat genome using STAR and quantify expression using a read counting strategy
2. Perform differential expression analysis of RNA-Seq data with DESeq2
3. Perform differential expression analysis of pre-normalized microarray expression data using limma
4. Map between Affymetrix and refSeq identifier systems
# Contributors

+ Daisy Wenyan Han daisyhan@bu.edu
+ Divya Sundaresan divyas3@bu.edu
+ Alec Jacobsen aggjacob@bu.edu
+ Emmanuel Saake esaake@bu.edu

# Repository Contents

Provide a brief description of each script/code file in this repo, what it does, and how to execute it
