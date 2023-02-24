source("gsea.R")

set.seed(200)

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

Y <- readRDS("./data/smoking_vs_non_smoking.RDS")
metadata <- readRDS("./data/smoking_vs_non_smoking_categories.RDS")
microbe_sets <- readRDS("./data/smoking_vs_non_smoking_microbe_sets.RDS")

X <- cbind(1, metadata)
W <- calculate_W(Y+0.1, X, 2, 0)
lfc <- ols_solution(W, X, 2)
# Process Microbe Sets
path_inds <- list()
path_names <- list()
for (p in names(microbe_sets)) {
	inds <- match(microbe_sets[[p]], row.names(Y))
	inds <- inds[!is.na(inds)]
	path_inds[p] <- list(inds)
	path_names[[p]] <- microbe_sets[[p]][microbe_sets[[p]]%in%row.names(Y)]
}

# GSEA-LFC
obs_scores <- c()
obs_p_vals <- c()
error_range <- seq(-1,1,0.05)
for(error in error_range) {
	print(error)
	gsea_res <- gsea_parallel(lfc+error, path_inds, 5000, 4)
	obs_scores <- cbind(obs_scores, gsea_res$obs_scores)
	obs_p_vals <- cbind(obs_p_vals, gsea_res$p_values)
}
row.names(obs_scores) <- names(path_inds)
row.names(obs_p_vals) <- names(path_inds)
colnames(obs_scores) <- error_range
colnames(obs_p_vals) <- error_range
write.table(obs_scores, "./data/figure_s5_obs_scores_lfc.txt")
write.table(obs_p_vals, "./data/figure_s5_obs_p_vals_lfc.txt")

# GSEA-LFC-S
obs_scores <- c()
obs_p_vals <- c()
for(error in error_range) {
	print(error)
	noise <- X[,2] * error
	noise_adj_W <- sweep(W, 2, noise, "+")
	gsea_res <- gsea_parallel_matrix(noise_adj_W, X, path_inds, 2, 5000, cores=4)
	obs_scores <- cbind(obs_scores, gsea_res$obs_scores)
	obs_p_vals <- cbind(obs_p_vals, gsea_res$p_values)
}
row.names(obs_scores) <- names(path_inds)
row.names(obs_p_vals) <- names(path_inds)
colnames(obs_scores) <- seq(-1,1,0.05)
colnames(obs_p_vals) <- seq(-1,1,0.05)
write.table(obs_scores, "./data/figure_s5_obs_scores_lfc_s.txt")
write.table(obs_p_vals, "./data/figure_s5_obs_p_vals_lfc_s.txt")


