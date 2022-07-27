library(fgsea)
source("gsea.R")

cores <- 7
niter <- 2000

pathways <- gmtPathways("c2.all.v7.4.symbols.gmt")
Y <- as.matrix(read.table("breast_tumor_counts.txt", row.names=1))
metadata <- read.table("breast_tumor_metadata.txt", header=T)
metadata[,2] <- factor(metadata[,2], level=c("Healthy", "Tumor"))
row.names(metadata) <- NULL
X <- model.matrix(~condition, metadata)

path_inds <- list()
path_names <- list()
for (p in names(pathways)) {
	if(!grepl("REACTOME_", p)) {
		next
	}
	inds <- match(pathways[[p]], row.names(Y))
	inds <- inds[!is.na(inds)]
	if((length(inds)<15)|(length(inds)>500)) {
		next
	}
	path_inds[p] <- list(inds)
	path_names[[p]] <- pathways[[p]][pathways[[p]]%in%row.names(Y)]
}

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

calculate_W <- function(Y, X, lfc_column, epsilon) {
	W_par <- t(t(Y)/colSums(Y))
	W_perp <- log(apply(W_par, 2, function(col) 1/gm_mean(col)))
	error <- epsilon * X[,lfc_column]
	W <- log(W_par)
	W <- sweep(W, 2, W_perp, "+")
	W <- sweep(W, 2, error, "+")
	W
}

W <- calculate_W(Y, X, 2, 0)

gsea_unweighted_res <- gsea_parallel_matrix(W, X, path_inds, 2, iterations=niter, cores=cores, unweighted=T)
gsea_cw_res <- gsea_parallel_matrix(W, X, path_inds, 2, iterations=niter, cores=cores, score_func=score_mean_Euclidean_weights)
gsea_res <- gsea_parallel_matrix(W, X, path_inds, 2, iterations=niter, cores=cores)

saveRDS(gsea_unweighted_res, "fig_s2_t1_data_gsea_unweighted.RDS")
saveRDS(gsea_cw_res, "fig_s2_t1_data_gsea_cw.RDS")
saveRDS(gsea_res, "fig_s2_t1_data_gsea.RDS")

