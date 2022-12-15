#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(fgsea)
library(MASS)
library(reshape2)
library(limma)
library(foreach)
library(doSNOW)
source("gsea.R")

pathways <- gmtPathways("c2.all.v7.4.symbols.gmt")
Y <- as.matrix(read.table("breast_tumor_counts.txt", row.names=1))
metadata <- read.table("breast_tumor_metadata.txt", row.names=1, header=1)
metadata[,1] <- factor(metadata[,1], level=c("Healthy", "Tumor"))
row.names(metadata) <- NULL
X <- model.matrix(~condition, metadata)

path_inds <- list()
path_names <- list()
for (p in names(pathways)) {
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

get_corr <- function(Ym, X, inds) {
	        sY <- Ym[inds,]
        B <- ginv(t(X)%*%X)%*%t(X)%*%t(sY)
	        Y_pred <- t(B)%*%t(X)
	        corr_mat <- cor(t(sY-Y_pred))
		        return(list(mean_corr=mean(corr_mat[lower.tri(corr_mat)]), corr_mat=corr_mat))
}

simple_camera <- function(LFC, path_names, X, lfc_column, inter.gene.cor) {
	  p_values <- c()
  count <- 1
    for (path_name in names(path_names)) {
	        inds <- path_names[[path_name]]
      G <- length(LFC)
          m1 <- length(inds)
          m2 <- G-m1
	      delta <- ( mean(LFC[inds]) - mean(LFC) ) * G / m2
	      s2 <- var(LFC)
	          s2p <- ( ((G-1)*s2)-(delta^2*m1*m2/G) ) / (G-2)
	          VIF <- 1 + ( (m1-1) * (inter.gene.cor[count]) )
		      correlation <- (VIF-1)/(m1-1)
		      T_stat <- delta / (sqrt(s2p * ((VIF/m1) + (1/m2))) )
		          dfs <- G-2
		          p_values <- c(p_values, 2*pmin(pt(T_stat, df=dfs), pt(T_stat, df=dfs,lower.tail=F)))
			      count <- count + 1
			    }
    p_values
}

full_cov <- function(W_par, W_perp) {
	  corr_mat <- matrix(0, nrow=nrow(W_par), ncol=nrow(W_par))
  for(i in 1:nrow(W_par)) {
	      for(j in 1:nrow(W_par)) {
		            corr_mat[i,j] <- cov(W_par[i,], W_par[j,]) + cov(W_par[i,], W_perp) + cov(W_perp, W_par[j,]) + cov(W_perp, W_perp)
      }
    }
    row.names(corr_mat) <- row.names(W_par)
    colnames(corr_mat) <- row.names(W_par)
      cov2cor(corr_mat)
}


W <- calculate_W(Y+0.1, X, 2, 0)

camera_lfc_results <- c()
GSEA_results <- c()

set.seed(100)
permutations <- replicate(5000, unname(sample(X[,2], replace=F)))
set.seed(NULL)

noise_val <- 0.225

#LFC <- ols_solution(W, X, 2)
#inter_gene_cors <- unname(unlist((lapply(path_inds, function(ind) get_corr(W, X, ind)$mean_corr))))
#GSEA_results <- gsea_parallel_matrix(W, X, path_inds, 2, permutation_matrix=permutations, cores=24, cw=F)$p_value
#saveRDS(GSEA_results, "GSEA_LFC_C2_0.RDS")
#camera_lfc_results <- simple_camera(LFC, path_names, X, 2, inter_gene_cors)
#saveRDS(camera_lfc_results, "CAMERA_LFC_C2_0.RDS")

for(i in 1:50) {
	noise_A <- runif(ncol(W)/2, -noise_val, noise_val)
	noise_A <- noise_A - mean(noise_A)
	noise_B <- runif(ncol(W)/2, -noise_val, noise_val)
	noise_B <- noise_B - mean(noise_B)
	noise <- c(noise_A, noise_B)
	noise_adj_W <- sweep(W, 2, noise, "+")
	LFC <- ols_solution(noise_adj_W, X, 2)
	inter_gene_cors <- unname(unlist((lapply(path_inds, function(ind) get_corr(noise_adj_W, X, ind)$mean_corr))))
	GSEA_results <- cbind(GSEA_results, gsea_parallel_matrix(noise_adj_W, X, path_inds, 2, permutation_matrix=permutations, cores=23)$p_value)
	camera_lfc_results <- cbind(camera_lfc_results, simple_camera(LFC, path_names, X, 2, inter_gene_cors))
}

row.names(GSEA_results) <- names(path_names)
write.table(GSEA_results, paste0("GSEA_LFC_C2_", noise_val, "_01.txt"), quote=F, sep="\t")
row.names(camera_lfc_results) <- names(path_names)
write.table(camera_lfc_results, paste0("CAMERA_LFC_C2_", noise_val, "_01.txt"), quote=F, sep="\t")

