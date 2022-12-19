library(fgsea)
source("gsea.R")

pathways <- gmtPathways("~/data/input/c2.all.v7.4.symbols.gmt")
Y <- as.matrix(read.table("./data/breast_tumor_counts.txt", row.names=1))
metadata <- read.table("./data/breast_tumor_metadata.txt", header=T)
metadata[,2] <- factor(metadata[,2], level=c("Healthy", "Tumor"))
head(metadata)
row.names(metadata) <- NULL
X <- model.matrix(~condition, metadata)
path_inds <- list()
path_names <- list()
for (p in names(pathways)) {
	if(!grepl("KEGG_", p)) {
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

length(path_inds)

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
obs_scores <- c()
obs_p_vals <- c()
for(e in seq(-0.25,0.25,0.025)) {
	print(e)
	noise <-c(rep(0, 92), rep(e, 92))
	noise_adj_W <- sweep(W, 2, noise, "+")
	gsea_res <- gsea_parallel_matrix(noise_adj_W, X, path_inds, 2, iterations=5000, cores=6)
	obs_scores <- cbind(obs_scores, gsea_res$obs_scores)
	obs_p_vals <- cbind(obs_p_vals, gsea_res$p_values)
}
row.names(obs_scores) <- names(path_inds)
row.names(obs_p_vals) <- names(path_inds)
write.table(obs_scores, "figure_s4_obs_scores.txt")
write.table(obs_p_vals, "figure_s4_obs_p_vals.txt")
