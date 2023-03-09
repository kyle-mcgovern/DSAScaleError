library(biomformat)

healthy_metadata <- read.table("GSE86354_GTEx_FeatureCounts.samples.txt")
thyroid_healthy_names <- gsub("-", ".", healthy_metadata[healthy_metadata[,2]=="Thyroid",1])
NAT_metadata <- read.table("GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt")
thyroid_NAT_names <- gsub("-", ".", NAT_metadata[NAT_metadata[,2]=="THCA",1])
healthy_counts <- read.table("GSE86354_GTEx_FeatureCounts.txt")
healthy_thyroid_counts <- healthy_counts[,thyroid_healthy_names]
NAT_counts <- read.table("GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt")
NAT_thyroid_counts <- NAT_counts[,thyroid_NAT_names]
paper_genes <- read.table("paper_genes.txt")[,1]
healthy_thyroid_counts_filt <- healthy_thyroid_counts[paper_genes,]
NAT_thyroid_counts_filt <- NAT_thyroid_counts[paper_genes,]
merged_counts <- merge(healthy_thyroid_counts_filt, NAT_thyroid_counts_filt, by="row.names", all=TRUE)
row.names(merged_counts) <- merged_counts[,1]
merged_counts <- merged_counts[,-1]
biom_counts <- make_biom(merged_counts)
write_biom(biom_counts, "thyroid_counts.biom")
metadata <- data.frame(sampleid=c(thyroid_healthy_names,thyroid_NAT_names), condition=c(rep("Healthy", length(thyroid_healthy_names)), rep("NAT", length(thyroid_NAT_names))))
write.table(metadata, "thyroid_metadata.txt", sep="\t", row.names=F, quote=F)
system("songbird multinomial --input-biom thyroid_counts.biom --metadata-file thyroid_metadata.txt --formula \"condition\" --epochs 5000 --summary-dir throid_songbird")
lfc <- read.table("./throid_songbird/differentials.tsv", header=T)
write.table(lfc[,c(1,3)], "thyroid_lfc.rnk", col.names=F, quote=F, sep="\t", row.names=F)
for(epsilon in seq(-1.2, 0.6, 0.05)) {
  epsilon <- round(epsilon, 2)
  print(epsilon)
  if (epsilon < 0) {
    epsilon_n <- paste0("n", abs(epsilon))
  } else {
    epsilon_n <- as.character(epsilon)
  }
  
  rnk_file <- paste0("thyroid_lfc_", epsilon_n)
  cat("#extra line \n", file=paste0(rnk_file, ".rnk"))
  tmp_lfc <- lfc[,c(1,3)]
  tmp_lfc[,2] <- tmp_lfc[,2] + epsilon
  write.table(tmp_lfc, paste0(rnk_file, ".rnk"), col.names=F, quote=F, sep="\t", row.names=F, append=T)
  cmd <- paste("~/GSEA_4.2.3/gsea-cli.sh GSEAPreranked -rnk", paste0(rnk_file, ".rnk"), "-gmx ~/data/input/h.all.v7.4.symbols.gmt -nperm 5000 -rpt_label", rnk_file, "-out thyroid_gsea \n")
  system(cmd)
}
all_data <- data.frame()
for(epsilon in seq(-1.2, 0.6, 0.05)) {
  epsilon <- round(epsilon, 2)
  if (epsilon < 0) {
    epsilon_n <- paste0("n", abs(epsilon))
  } else {
    epsilon_n <- as.character(epsilon)
  }
  dir_name <- list.files(path="./thyroid_gsea", pattern=paste0("thyroid_lfc_", epsilon_n, ".GseaPreranked", "[0-9]*"), include.dirs=T)
  dir <- paste0("./thyroid_gsea/", dir_name, "/")
  tsv_files <- list.files(path=dir, pattern="gsea_report_for.*\\.tsv")
  for (tsv_file in tsv_files) {
    df <- read.table(paste0(dir, tsv_file), sep="\t", header=T)[,c("NAME", "NES", "FDR.q.val")]
    df[df=="---"] <- NA
    df$NES <- as.numeric(df$NES)
    df$FDR.q.val <- as.numeric(df$FDR.q.val)
    df$epsilon <- epsilon
    all_data <- rbind(all_data, df)
  }
}
write.table(all_data, "thyroid_gsea_results.tsv", row.names=F, sep="\t", quote=F)


