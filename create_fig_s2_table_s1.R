library(fgsea)
library(ggplot2)
library(ggpubr)
source("gsea.R")

theme_Publication <- function(base_size=14, base_family="helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

gsea_res <- readRDS("fig_s2_t1_data_gsea.RDS")
gsea_cw_res <- readRDS("fig_s2_t1_data_gsea_cw.RDS")
gsea_unw_res <- readRDS("fig_s2_t1_data_gsea_unweighted.RDS")

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

# WP_DNA_IRDOUBLE_STRAND_BREAKS_DSBS_AND_CELLULAR_RESPONSE_VIA_ATM
# WP_G1_TO_S_CELL_CYCLE_CONTROL

gsea_p_vals <- gsea_res$p_value
gsea_cw_p_vals <- gsea_cw_res$p_value
gsea_unw_p_vals <- gsea_unw_res$p_value
names(gsea_p_vals) <- names(path_names)
names(gsea_cw_p_vals) <- names(path_names)
names(gsea_unw_p_vals) <- names(path_names)
res <- cbind(gsea_p_vals, gsea_cw_p_vals, gsea_unw_p_vals)
res[(res[,1]<0.05)&(res[,2]<0.05)&(res[,3]>0.05),]
res[(res[,1]>0.05)&(res[,2]>0.05)&(res[,3]<0.05),]
# REACTOME_AMINE_LIGAND_BINDING_RECEPTORS

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

res[c("REACTOME_AMINE_LIGAND_BINDING_RECEPTORS", "REACTOME_CYCLIN_A_B1_B2_ASSOCIATED_EVENTS_DURING_G2_M_TRANSITION", "REACTOME_MITOTIC_PROPHASE", "REACTOME_INTERCONVERSION_OF_NUCLEOTIDE_DI_AND_TRIPHOSPHATES", "REACTOME_INTRINSIC_PATHWAY_FOR_APOPTOSIS", "REACTOME_SUMOYLATION_OF_DNA_REPLICATION_PROTEINS", "REACTOME_OLFACTORY_SIGNALING_PATHWAY"),]

plot_one <- function(p_name) {
	W <- calculate_W(Y, X, 2, 0)
	lfc <- ols_solution(W, X, 2)
	inds <- path_inds[[p_name]]
	sorted_obj <- sort(lfc, decreasing=T, index.return=T)
	lfc <- sorted_obj$x
	idx <- sort(which(sorted_obj$ix%in%(inds)))

	fill <- rep("Not In Set", length(lfc))
	fill[idx] <- "In Set"
	df <- data.frame(x=1:length(lfc), y=lfc, fill=fill)
	g <- ggplot(df, aes(x=x, y=y, color=fill, shape=fill)) + geom_point(size=1.25, stroke=1.25) + theme_bw()
	g <- g + geom_point(data=data.frame(x=idx,y=lfc[idx]), aes(x=x, y=y), size=1.25, stroke=1.25,shape=3, color="#000000")
	g <- g + xlab("Entity List Rank") + ylab("Log Fold Change")
	g <- g + scale_color_manual(values=c("#000000", "#CCCCCC"))
	g <- g + scale_shape_manual(values=c(3, 1))
	g <- g + geom_segment(data=data.frame(x=df$x[idx]), aes(x=x, xend=x, y=-Inf, yend=min(lfc)+0.001), inherit.aes=F, color="red")
	g <- g + theme(legend.title=element_blank())
	g
}

g1 <- plot_one("REACTOME_CYCLIN_A_B1_B2_ASSOCIATED_EVENTS_DURING_G2_M_TRANSITION")
g2 <- plot_one("REACTOME_OLFACTORY_SIGNALING_PATHWAY")

g1 <- g1 + ggtitle("Cyclin A B1 B2 AEDGMT")
g2 <- g2 + ggtitle("Olfactory Signaling Pathway")
gar <- ggarrange(g1, g2, nrow=2)
ggsave("~/data/output/figure_s2.png", gar, units="in", height=7, width=10)


