#Package
library(SeqVarTools)
library(SeqArray)
library(plyr)
library(tidyverse)
library(VariantAnnotation)
library(GenomicFeatures)
library(ape)
library(viridis)
library(readxl)
library(ggpubr)
library(SNPRelate)
library(scales)
library(flextable)
library(officer)
library(ggrepel)
library(MASS)
library(cowplot)
library(pheatmap)
library(FactoMineR)
library(reshape)
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(Hmisc)
library(cluster)
library(dplyr)
library(ggplot2)
library(readr)
library(Rtsne)
library(fpc)
library(MASS)
library(ComplexHeatmap)
library(randomcoloR)
library(olsrr)
library(ggExtra)
library(vcfR)
library(ROCR)
#Plasmodium vivax SNP data
#SNPs were called using the standard best practice from Genome Analysis Toolkit (GATK) version 4.0.12.0 implemented in nextflow (https://github.com/gatk-workflows/gatk4-germline-snps-indels).
#########################################SNP PCA value PC1-PC10#########################################
#Get the genotype from the SNP genotype data and perform principal component analysis of SNP genotype data using prcomp function and get the first ten PC values
#SNP PCA PC1 ~ PC10 value
dist_SNP1 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP1)[2] <- "PCvalue"
dist_SNP2 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP2)[2] <- "PCvalue"
dist_SNP3 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP3)[2] <- "PCvalue"
dist_SNP4 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP4)[2] <- "PCvalue"
dist_SNP5 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP5)[2] <- "PCvalue"
dist_SNP6 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP6)[2] <- "PCvalue"
dist_SNP7 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP7)[2] <- "PCvalue"
dist_SNP8 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP8)[2] <- "PCvalue"
dist_SNP9 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP9)[2] <- "PCvalue"
dist_SNP10 <- read.csv("SNP/SNP_PCAvalue/PC1.csv")
colnames(dist_SNP10)[2] <- "PCvalue"
#########################################STR Motif count for each sample each loci#########################################
#Get the STR genotype data from the HipSTR output and qccording to the genotype allele and the reference motif get the motif count for each sample each loci.
#23146 STR loci   
#174 samples
motif <- read.csv("STR/Motif_count.csv")
rownames(motif) <- motif$STRid
motif <- motif[,-1]
motif <- t(motif)
motif_new <- motif
#########################################Test the correlation of STR motif count with SNP PC1 ~ PC10 value#########################################
#TRUE correlation
#Test the correlation for each STR loci in only PC1
#PC1
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  test1
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_1_T.csv",row.names = FALSE)

#Test the correlation for each STR loci in both PC1 and PC2, and choose the maximum one
#PC1-PC2
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  test <- max(test1,test2,na.rm = TRUE)
  test
})


correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_2_T.csv",row.names = FALSE)


#PC1-PC3
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  
  
  test <- max(test1,test2,test3,na.rm = TRUE)
  test
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_3_T.csv",row.names = FALSE)


#PC1-PC4
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  
  test <- max(test1,test2,test3,test4,na.rm = TRUE)
  test
})


correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_4_T.csv",row.names = FALSE)




#PC1-PC5
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  test <- max(test1,test2,test3,test4,test5,na.rm = TRUE)
  test
})


correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_5_T.csv",row.names = FALSE)




#PC1-PC6
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  
  test <- max(test1,test2,test3,test4,test5,test6,na.rm = TRUE)
  test
})


correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_6_T.csv",row.names = FALSE)


#PC1-PC7
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  dist_SNP_STR_7 <- dist_SNP7 %>% 
    left_join(dist_STR,by="Sample")
  test7 <- cor(dist_SNP_STR_7$PCvalue, dist_SNP_STR_7$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test7 <- test7*test7
  
  test <- max(test1,test2,test3,test4,test5,test6,test7,na.rm = TRUE)
  test
})


correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_7_T.csv",row.names = FALSE)




#PC1-PC8
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  dist_SNP_STR_7 <- dist_SNP7 %>% 
    left_join(dist_STR,by="Sample")
  test7 <- cor(dist_SNP_STR_7$PCvalue, dist_SNP_STR_7$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test7 <- test7*test7
  
  dist_SNP_STR_8 <- dist_SNP8 %>% 
    left_join(dist_STR,by="Sample")
  test8 <- cor(dist_SNP_STR_8$PCvalue, dist_SNP_STR_8$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test8 <- test8*test8
  
  
  test <- max(test1,test2,test3,test4,test5,test6,test7,test8,na.rm = TRUE)
  test
})


correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_8_T.csv",row.names = FALSE)


#PC1-PC9
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  dist_SNP_STR_7 <- dist_SNP7 %>% 
    left_join(dist_STR,by="Sample")
  test7 <- cor(dist_SNP_STR_7$PCvalue, dist_SNP_STR_7$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test7 <- test7*test7
  
  dist_SNP_STR_8 <- dist_SNP8 %>% 
    left_join(dist_STR,by="Sample")
  test8 <- cor(dist_SNP_STR_8$PCvalue, dist_SNP_STR_8$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test8 <- test8*test8
  
  dist_SNP_STR_9 <- dist_SNP9 %>% 
    left_join(dist_STR,by="Sample")
  test9 <- cor(dist_SNP_STR_9$PCvalue, dist_SNP_STR_9$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test9 <- test9*test9
  
  test <- max(test1,test2,test3,test4,test5,test6,test7,test8,test9,na.rm = TRUE)
  test
})


correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_9_T.csv",row.names = FALSE)



#PC1-PC10
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  dist_SNP_STR_7 <- dist_SNP7 %>% 
    left_join(dist_STR,by="Sample")
  test7 <- cor(dist_SNP_STR_7$PCvalue, dist_SNP_STR_7$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test7 <- test7*test7
  
  dist_SNP_STR_8 <- dist_SNP8 %>% 
    left_join(dist_STR,by="Sample")
  test8 <- cor(dist_SNP_STR_8$PCvalue, dist_SNP_STR_8$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test8 <- test8*test8
  
  dist_SNP_STR_9 <- dist_SNP9 %>% 
    left_join(dist_STR,by="Sample")
  test9 <- cor(dist_SNP_STR_9$PCvalue, dist_SNP_STR_9$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test9 <- test9*test9
  
  dist_SNP_STR_10 <- dist_SNP10 %>% 
    left_join(dist_STR,by="Sample")
  test10 <- cor(dist_SNP_STR_10$PCvalue, dist_SNP_STR_10$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test10 <- test10*test10
  
  test <- max(test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,na.rm = TRUE)
  test
})

correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_TRUE/correlation_STR_10_T.csv",row.names = FALSE)
####FALSE correlation (permutation test where the population labels were permuted between samples to derive a null distribution to determine an appropriate correlation threshold that maximised the difference between high and low quality STRs. )
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  test1
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_1_F.csv",row.names = FALSE)

#PC2
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  test <- max(test1,test2,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_2_F.csv",row.names = FALSE)

#PC3
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  #PC2
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  #PC3
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  test <- max(test1,test2,test3,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_3_F.csv",row.names = FALSE)




#PC4
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  #PC2
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  #PC3
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  #PC4
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  test <- max(test1,test2,test3,test4,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_4_F.csv",row.names = FALSE)





#PC5
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  #PC2
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  #PC3
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  #PC4
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  
  #PC5
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  test <- max(test1,test2,test3,test4,test5,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_5_F.csv",row.names = FALSE)


#PC6
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  #PC2
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  #PC3
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  #PC4
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  
  #PC5
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  #PC6
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  test <- max(test1,test2,test3,test4,test5,test6,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_6_F.csv",row.names = FALSE)





#PC7
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  #PC2
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  #PC3
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  #PC4
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  
  #PC5
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  #PC6
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  #PC7
  dist_SNP_STR_7 <- dist_SNP7 %>% 
    left_join(dist_STR,by="Sample")
  test7 <- cor(dist_SNP_STR_7$PCvalue, dist_SNP_STR_7$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test7 <- test7*test7
  
  
  test <- max(test1,test2,test3,test4,test5,test6,test7,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_7_F.csv",row.names = FALSE)



#PC8
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  #PC2
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  #PC3
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  #PC4
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  
  #PC5
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  #PC6
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  #PC7
  dist_SNP_STR_7 <- dist_SNP7 %>% 
    left_join(dist_STR,by="Sample")
  test7 <- cor(dist_SNP_STR_7$PCvalue, dist_SNP_STR_7$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test7 <- test7*test7
  
  #PC8
  dist_SNP_STR_8 <- dist_SNP8 %>% 
    left_join(dist_STR,by="Sample")
  test8 <- cor(dist_SNP_STR_8$PCvalue, dist_SNP_STR_8$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test8 <- test8*test8
  
  
  test <- max(test1,test2,test3,test4,test5,test6,test7,test8,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_8_F.csv",row.names = FALSE)




#PC9
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  #PC2
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  #PC3
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  #PC4
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  
  #PC5
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  #PC6
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  #PC7
  dist_SNP_STR_7 <- dist_SNP7 %>% 
    left_join(dist_STR,by="Sample")
  test7 <- cor(dist_SNP_STR_7$PCvalue, dist_SNP_STR_7$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test7 <- test7*test7
  
  #PC8
  dist_SNP_STR_8 <- dist_SNP8 %>% 
    left_join(dist_STR,by="Sample")
  test8 <- cor(dist_SNP_STR_8$PCvalue, dist_SNP_STR_8$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test8 <- test8*test8
  
  #PC9
  dist_SNP_STR_9 <- dist_SNP9 %>% 
    left_join(dist_STR,by="Sample")
  test9 <- cor(dist_SNP_STR_9$PCvalue, dist_SNP_STR_9$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test9 <- test9*test9
  
  
  
  test <- max(test1,test2,test3,test4,test5,test6,test7,test8,test9,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_9_F.csv",row.names = FALSE)

#PC10
correlation <- lapply(1:ncol(motif), function(i){
  dist_STR <- data.frame(Sample=rownames(motif),STRRepeat=motif[,i])
  reorder_STR <- dist_STR[sample(nrow(dist_STR)),1]
  dist_STR <- data.frame(Sample=reorder_STR,STRRepeat=dist_STR$STRRepeat)
  dist_SNP_STR <- dist_SNP1 %>% 
    left_join(dist_STR,by="Sample")
  test1 <- cor(dist_SNP_STR$PCvalue, dist_SNP_STR$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test1 <- test1*test1
  #PC2
  dist_SNP_STR_2 <- dist_SNP2 %>% 
    left_join(dist_STR,by="Sample")
  test2 <- cor(dist_SNP_STR_2$PCvalue, dist_SNP_STR_2$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test2 <- test2*test2
  #PC3
  dist_SNP_STR_3 <- dist_SNP3 %>% 
    left_join(dist_STR,by="Sample")
  test3 <- cor(dist_SNP_STR_3$PCvalue, dist_SNP_STR_3$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test3 <- test3*test3
  #PC4
  dist_SNP_STR_4 <- dist_SNP4 %>% 
    left_join(dist_STR,by="Sample")
  test4 <- cor(dist_SNP_STR_4$PCvalue, dist_SNP_STR_4$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test4 <- test4*test4
  
  
  #PC5
  dist_SNP_STR_5 <- dist_SNP5 %>% 
    left_join(dist_STR,by="Sample")
  test5 <- cor(dist_SNP_STR_5$PCvalue, dist_SNP_STR_5$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test5 <- test5*test5
  
  #PC6
  dist_SNP_STR_6 <- dist_SNP6 %>% 
    left_join(dist_STR,by="Sample")
  test6 <- cor(dist_SNP_STR_6$PCvalue, dist_SNP_STR_6$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test6 <- test6*test6
  
  #PC7
  dist_SNP_STR_7 <- dist_SNP7 %>% 
    left_join(dist_STR,by="Sample")
  test7 <- cor(dist_SNP_STR_7$PCvalue, dist_SNP_STR_7$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test7 <- test7*test7
  
  #PC8
  dist_SNP_STR_8 <- dist_SNP8 %>% 
    left_join(dist_STR,by="Sample")
  test8 <- cor(dist_SNP_STR_8$PCvalue, dist_SNP_STR_8$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test8 <- test8*test8
  
  #PC9
  dist_SNP_STR_9 <- dist_SNP9 %>% 
    left_join(dist_STR,by="Sample")
  test9 <- cor(dist_SNP_STR_9$PCvalue, dist_SNP_STR_9$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test9 <- test9*test9
  
  #PC10
  dist_SNP_STR_10 <- dist_SNP10 %>% 
    left_join(dist_STR,by="Sample")
  test10 <- cor(dist_SNP_STR_10$PCvalue, dist_SNP_STR_10$STRRepeat,method = "spearman",use = "pairwise.complete.obs")
  test10 <- test10*test10
  
  test <- max(test1,test2,test3,test4,test5,test6,test7,test8,test9,test10,na.rm = TRUE)
  test
  
})
correlation <- do.call(rbind,correlation)
correlation_all <- data.frame(STRid=motif_new$STRid,correlation)
head(correlation_all)
write.csv(correlation_all,"SNP_STR_Correlation/Correlation_FALSE_Permutationtest/correlation_STR_10_F.csv",row.names = FALSE)

###################################################Built the model###################################################
#The optimal cut-point of the correlation to distinguish high-quality or low-quality STRs was determined using the resampling permutation test where the population labels were permuted between samples to derive a null distribution to determine an appropriate correlation threshold that maximised the difference between high and low quality STRs. This was determined using the AUC.
#Threshold Mononucleotide 0.048065497, Polynucleotide 0.047691231
#Examination of Spearman's correlation coefficients (R2) suggested that selecting the first five SNP principal components (PCs) were sufficient to capture the signal STRs. 
#Model(mononucleotide STR-Model)
#Set the correlation threhold
correlation_all <- read.csv("SNP_STR_Correlation/Correlation_TRUE/correlation_STR_5_T.csv")
#Choose the right threhold
correlation_all <- correlation_all %>% 
  mutate(STRQuality=ifelse(correlation>0.048065497, "Good","Bad"))
correlation_all$STRQuality[is.na(correlation_all$STRQuality)==TRUE] <- "Bad"
#change to factor
data <- read.csv("STR/STR_Feature.csv")
data <- correlation_all %>% 
  left_join(data,by="STRid")
data$STRQuality <- ifelse(data$STRQuality=="Good",1,0)
data <- data %>% na.omit(data)
data <- data %>% 
  filter(Length==1)
#Scale
X =data[,4:ncol(data)]
X.scaled = scale(X, center= TRUE, scale=TRUE)
colMeans(X.scaled)
var(X.scaled)
data <- data.frame(data[,1:3],X.scaled)
#######################################glm###############################################
#First build the model using all the STR feature
model <- glm(STRQuality ~ complexity+ GC_flank+GC_STR + GC_diff+  MaxHe + MeanHe + MinHe + Heterozygosity+JostD + Repeatunits+ mean_indel+mean_Posterior + mean_stutter+ Missingness, data = data,family = binomial)
summary(model)
#Model selection in the multivariable regression models was employed for stepwise regression analysis based on the Akaike information criterion (AIC) using the R package MASS (Version704 7.3-51.5)
#Make predictions
step.model <- model %>% stepAIC(trace = FALSE)
coef(step.model)
#Significantly associated STR feature
#Rebuild the model using the selected STR feature
model <- glm(STRQuality ~ GC_flank+GC_STR + GC_diff + MaxHe + MeanHe + MinHe + Heterozygosity+ Repeatunits+ mean_Posterior + mean_stutter, data = data,family = binomial)
#Calculate the 95% confidence level
confint(model)
#Test the model performance 
#The effectiveness of the prediction was evaluated by calculating the AUC on the ROC curve. 
train <- data
#test data
test <- data
summary(fit <- glm(STRQuality ~ GC_flank+GC_STR + GC_diff + MaxHe + MeanHe + MinHe + Heterozygosity+ Repeatunits+ mean_Posterior + mean_stutter, data = train,family = binomial))
prob <- predict(fit, newdata=test, type="response")
pred <- prediction(prob, test$STRQuality)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc
#Predict the STR quailty
predicted.data <- data.frame(probability.of.hd=model$fitted.values,STRQuality=data$STRQuality,STRid=data$STRid)
predicted.data <- predicted.data[order(predicted.data$probability.of.hd,decreasing=FALSE),]
predicted.data$rank <- 1:nrow(predicted.data)
predicted.data$STRQuality <- ifelse(predicted.data$STRQuality== 1,"Good","Bad")
#Model(polynucleotide STR-Model)
#The polynucleotide STR-Model build step is same as the mononucleotide STR-Model
######################################Select the high quality STRs##########################################
#To select the optimal cut-off value to remove low-quality STRs, the predicted probabilities are sorted into five bins ([0, 0.2), [0.2, 0.4), [0.4, 0.6), [0.6, 0.8), [0.8, 1]). For each bin, we randomly selected 500 STRs and calculated the correlation of sample pairwise distances of STR and SNP based on PCA analysis, and repeated this 100 times to calculate the mean value of correlation.