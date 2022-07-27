source("gsea.R")
source("utils.R")

### TEST GSEA ###

# Simulate dataset
print("TESTING gsea")
set.seed(1056)
lfc <- sort(round(rnorm(1000), 5), decreasing=T)
names <- rand_strs(1000)
inds_A <- c(50,56,67,87,99,123,153,177,200,201,220,410,601,700,706,809)
inds_B <- c(54,156,896,999,876,789,567,983,653,902,306,402,782,503,740,799)
set.seed(NULL)

# Write required files for GSEA

# RNK File
cat("# Random 1000 data \n",file="GSEA_Test_1.rnk")
write.table(cbind(names, lfc), "GSEA_Test_1.rnk", sep="\t", append=T, col.names=F, row.names=F, quote=F)

# GMT File
write.table(rbind(c("Test_A", "None", names[inds_A]), c("Test_B", "None", names[inds_B])), "GSEA_Test_1.gmt", sep="\t", col.names=F, row.names=F, quote=F)

# Code to run GSEA-LFC to generate test data
# ~/GSEA_4.2.3/gsea-cli.sh GSEAPreranked -rnk GSEA_Test_1.rnk -gmx GSEA_Test_1.gmt -nperm 100000 -rpt_label GSEA_test_1 -out GSEA_test_1

# Results from GSEA pipeline
GSEA_test_1_A_es <- 0.636
GSEA_test_1_B_es <- -0.482
GSEA_test_1_A_pval <- 0.002
GSEA_test_1_B_pval <- 0.065

# Test score and p-value
print("EXPECTED")
print(c(0.636, 0.002))
print(c(-0.482, 0.065))
print("OBSERVED")
gsea(lfc, inds_A, 100000)
gsea(lfc, inds_B, 100000)
# Test score and p-value with null_es precomputed
null_es <- c()
for(i in 1:100000) {
	shuffled_inds <- sample(1:length(lfc), length(inds_A), replace=F)
	null_es <- c(null_es, gsea_base(lfc, shuffled_inds))
}
print("EXPECTED")
print(c(0.636, 0.002))
print(c(-0.482, 0.065))
print("OBSERVED")
gsea(lfc, inds_A, null_scores=null_es)
gsea(lfc, inds_B, null_scores=null_es)

print("OBSERVED gsea_parallel")
gsea_parallel(lfc, list(inds_A=inds_A, inds_B=inds_B), 100000, 6)


### TEST MEAN EUCLIDEAN WEIGHT FUNCTION ###

print("TESTING score_mean_Euclidean_weights")

test_mean_euc <- c(5,4,3,-1,-2)
expected_mean_euc <- c(16/4,13/4,12/4,16/4)
print("OBSERVED")
score_mean_Euclidean_weights(c(5,4,3,-1), test_mean_euc)
print("EXPECTED")
expected_mean_euc

### TEST GSEA PARALLEL FOR COLUMN SHUFFLING ###

# ~/GSEA_4.2.3/gsea-cli.sh GSEA -res GSEA_2_test_table_clrs.gct -gmx GSEA_Test_2.gmt -cls GSEA_2_test_table_clrs.cls -nperm 10000 -permute phenotype -metric Diff_of_Classes -rpt_label GSEA_test_2 -out GSEA_test_2
# Results GSEA_COL_TEST_1 ES -0.27 p-value 0.918 GSEA_COL_TEST_2 ES -0.56 p-value 0.651
print ("TESTING gsea_parallel_matrix")
S <- read.table("GSEA_2_test_table_clrs.txt", row.names=1)
X <- cbind(rep(1,20), c(rep(1,10), rep(0,10)))
GSEA_COL_TEST_1 <- match(c("ZNF205", "UQCR11", "UBE2E2", "DDX11L2", "CD82", "C9orf102", "C7orf10", "C11orf24", "ARL4C", "ARID1B", "AHR", "DIAPH2", "EPB41L4A-AS1", "EPHA5", "FBLN7", "FCRL2"), row.names(S))
GSEA_COL_TEST_2 <- match(c("JUNB", "CXCL2", "ATF3", "NFKBIA", "PTGS2", "CXCL1", "PLK2", "IER3", "CD83", "CCL20", "CXCL3", "MAFF", "NFKB2", "TNFAIP2", "KLF6", "BIRC3", "PLAUR", "ZFP36", "ICAM1", "JUN", "EGR3", "IL1B"), row.names(S))
print("EXPECTED")
print(c("GSEA_COL_TEST_1, -0.27, 0.918"))
print(c("GSEA_COL_TEST_2, -0.56, 0.651"))
print("OBSERVED")
gsea_parallel_matrix(S, X, list(GSEA_COL_TEST_1=GSEA_COL_TEST_1, GSEA_COL_TEST_2=GSEA_COL_TEST_2), 2, iterations=10000, cores=6)


