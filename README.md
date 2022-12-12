# DSAScaleError

This repo contains code used to generate figures/perform analyses in the paper "Addressing Erroneous Scale Assumptions in Gene and Microbe Enrichment Analysis"

## Files

gsea.R is the core file, containing the code to run GSEA-LFC and GSEA-CW parallelized.
There exists a simple set of tests for this file.

## Requirements

* Songbird: https://github.com/biocore/songbird
* GSEA CLI tool: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html

## Data Access

The data files for the analyses are too big to include in this repo. Instead they can be found at the following links (warning, these links download directly and can be fairly large: a gig or more).

### Gene Sets

* http://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/7.4/msigdb_v7.4_files_to_download_locally.zip

### Thyroid Tissue Data

* GSE86354_GTEx_FeatureCounts.samples.txt: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86354&format=file&file=GSE86354%5FGTEx%5FFeatureCounts%2Esamples%2Etxt%2Egz
* GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE62944&format=file&file=GSE62944%5F06%5F01%5F15%5FTCGA%5F24%5FNormal%5FCancerType%5FSamples%2Etxt%2Egz
* GSE86354_GTEx_FeatureCounts.txt: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86354&format=file&file=GSE86354%5FGTEx%5FFeatureCounts%2Etxt%2Egz
* GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1697009&format=file&file=GSM1697009%5F06%5F01%5F15%5FTCGA%5F24%2Enormal%5FRsubread%5FFeatureCounts%2Etxt%2Egz

### Breast Tissue Data

* GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE62944&format=file&file=GSE62944%5F06%5F01%5F15%5FTCGA%5F24%5FCancerType%5FSamples%2Etxt%2Egz
* GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1536837&format=file&file=GSM1536837%5F06%5F01%5F15%5FTCGA%5F24%2Etumor%5FRsubread%5FFeatureCounts%2Etxt%2Egz

