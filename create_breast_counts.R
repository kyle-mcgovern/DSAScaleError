library(biomformat)

healthy_metadata <- read.table("GSE86354_GTEx_FeatureCounts.samples.txt")
breast_healthy_names <- gsub("-", ".", healthy_metadata[healthy_metadata[,2]=="Breast",1])
tumor_metadata <- read.table("GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt")
breast_tumor_names <- gsub("-", ".", tumor_metadata[tumor_metadata[,2]=="BRCA",1])
healthy_counts <- read.table("GSE86354_GTEx_FeatureCounts.txt")
healthy_breast_counts <- healthy_counts[,breast_healthy_names]
tumor_counts <- read.table("GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt")
set.seed(1234)
breast_tumor_names <- sample(breast_tumor_names, 92)
set.seed(NULL)
tumor_breast_counts <- tumor_counts[,breast_tumor_names]

healthy_breast_counts_filt <- healthy_breast_counts[rowSums(healthy_breast_counts==0)<=89,]
tumor_breast_counts_filt <- tumor_breast_counts[rowSums(tumor_breast_counts==0)<=89,]
merged_counts <- merge(healthy_breast_counts_filt, tumor_breast_counts_filt, by="row.names", all=FALSE)
row.names(merged_counts) <- merged_counts[,1]
merged_counts <- merged_counts[,-1]
dim(merged_counts)
merged_counts[1:5,1:5]
write.table(merged_counts+0.1, "breast_tumor_counts.txt")
metadata <- data.frame(sampleid=c(breast_healthy_names, breast_tumor_names), condition=c(rep("Healthy", length(breast_healthy_names)), rep("Tumor", length(breast_tumor_names))))
write.table(metadata, "breast_tumor_metadata.txt", sep="\t", row.names=F, quote=F)

