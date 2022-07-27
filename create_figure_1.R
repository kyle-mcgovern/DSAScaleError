source("gsea.R")
library(ggplot2)
library(ggpubr)


n <- 60
set.seed(12)
lfc <- sort(rnorm(n), decreasing=T)
set.seed(NULL)
inds <- c(3,5,7,12,15,17,20,25,40,45,48,55)
in_set <- rep("Not In Set", n)
in_set[inds] <- "In Set"

running_sum <- function(lfc, inds) {
	gsea_hit <- rep(0, length(lfc))
	gsea_hit[inds] <- abs(lfc[inds])/sum(abs(lfc[inds]))
	gsea_miss <- rep(0, length(lfc))
	gsea_miss[-inds] <- 1/(length(lfc)-length(inds))
	gsea_running_sum <- cumsum(gsea_hit-gsea_miss)
	list(running=c(0, gsea_running_sum), sup=gsea_running_sum[which.max(abs(gsea_running_sum))], ind=which.max(abs(gsea_running_sum)))
}

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

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}


print("STARTING PLOTS")
### PART A ###

df <- data.frame(x=1:length(lfc), y=lfc, in_set=in_set, error="0")
df <- rbind(df, data.frame(x=1:length(lfc), y=lfc+1, in_set=in_set, error="1"))
df <- rbind(df, data.frame(x=1:length(lfc), y=lfc-1, in_set=in_set, error="-1"))
g1 <- ggplot(df, aes(x=x, y=y, shape=in_set, color=error))
g1 <- g1 + geom_hline(yintercept=0, color="grey", linetype="dashed", size=1) + geom_point(size=1.25, stroke=1.25)
g1 <- g1 + theme_Publication()
g1 <- g1 + xlab("Entity List Rank")
g1 <- g1 + ylab("Log Fold Change")
g1 <- g1 + scale_shape_manual(values=c(3, 21))
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1 + scale_colour_manual(labels=c("ε=-1", "ε=0", "ε=1"), values=c("#E69F00", "#000000", "#0072B2"))
g1 <- g1 + theme(legend.position="right") + theme(legend.direction='vertical')
g1 <- g1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print("PART A DONE")

### PART B ###

res_0 <- gsea(lfc, inds, iters=5000)
res_1 <- gsea(lfc+1, inds, iters=5000)
res_n1 <- gsea(lfc-1, inds, iters=5000)

df <- data.frame(x=0:length(lfc), y=running_sum(lfc, inds)$running, error="0")
df <- rbind(df, data.frame(x=0:length(lfc), y=running_sum(lfc+1, inds)$running, error="1"))
df <- rbind(df, data.frame(x=0:length(lfc), y=running_sum(lfc-1, inds)$running, error="-1"))
df2 <- data.frame(x=c(running_sum(lfc, inds)$ind, running_sum(lfc+1, inds)$ind, running_sum(lfc-1, inds)$ind),
                  y=c(running_sum(lfc, inds)$sup, running_sum(lfc+1, inds)$sup, running_sum(lfc-1, inds)$sup),
                  color=c("0", "1", "-1"))

g2 <- ggplot(df, aes(x=x, y=y, color=error)) + geom_hline(yintercept=0, color="grey", linetype="dashed", size=1) + geom_line(size=1.75)
g2 <- g2 + geom_segment(data=df2, aes(x=x, xend=x, y=0, yend=y, color=color), inherit.aes=F, linetype="dotted", size=1.5)
g2 <- g2 + theme_Publication()
g2 <- g2 + xlab("Entity List Rank")
g2 <- g2 + ylab("Running Sum")
g2 <- g2 + scale_colour_manual(values=c("#E69F00", "#000000", "#0072B2"))
g2 <- g2 + theme(legend.position = "none")
g2 <- g2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
l <- paste0(round(res_1$obs_score, 2), " (", round(res_1$p_value, 2), ")")
g2 <- g2 + annotate(geom="text", x=12, y=res_1$obs_score+0.05, label=l, color="#0072B2", hjust=0, size=3.5)
l <- paste0(round(res_n1$obs_score, 2), " (", round(res_n1$p_value, 2), ")")
g2 <- g2 + annotate(geom="text", x=40, y=res_n1$obs_score-0.01, label=l, color="#E69F00", hjust=0, size=3.5)
l <- paste0(round(res_0$obs_score, 2), " (", round(res_0$p_value, 2), ")")
g2 <- g2 + annotate(geom="text", x=0, y=res_0$obs_score+0.09, label=l, color="#000000", hjust=0, size=3.5)
g2 <- g2 + geom_segment(aes(x=25.5, y=res_1$obs_score+0.05, xend=27.5, yend=res_1$obs_score+0.05), arrow = arrow(length = unit(0.1, "cm")), color="#808080")
g2 <- g2 + annotate(geom="text", x=28, y=res_1$obs_score+0.05, label="Enrichment Score (p-value)", color="#808080", hjust=0, size=3.5)
print("PART B DONE")

### PART C ###

perms <- replicate(10000, sample(1:length(lfc), length(inds), replace=F))
perms_null <- replicate(5000, sample(1:length(lfc), length(inds), replace=F))

p_mat <- c()
score_mat <- c()
for (i in seq(-2, 2, 0.1)) {
	print(i)
	p_vals <- c()
	obs_scores <- c()
	null_scores <- c()
	for(k in 1:ncol(perms_null)) {
		null_scores <- c(null_scores, gsea_base(lfc+i, perms_null[,k]))
	}
	for (j in 1:ncol(perms)) {
		gsea_res <- gsea(lfc+i, perms[,j], null_scores=null_scores)
		p_vals <- c(p_vals, gsea_res$p_value)
		obs_scores <- c(obs_scores, gsea_res$obs_score)
	}
	p_mat <- cbind(p_mat, p_vals)
	score_mat <- cbind(score_mat, obs_scores)
}
colnames(p_mat) <- round(seq(-2, 2, 0.1),1)
colnames(score_mat) <- round(seq(-2, 2, 0.1),1)


zero_p_vals <- p_mat[,21]
zero_scores <- score_mat[,21]
ppvs <- c()
for (i in 1:ncol(p_mat)) {
  tp <- sum( (zero_p_vals<=0.05)&(p_mat[,i]<=0.05)&(sign(zero_p_vals)==sign(p_mat[,i])) )
  fp <- sum( (zero_p_vals<=0.05)&(p_mat[,i]>0.05) )
  ppv <- tp/(tp+fp)
  ppvs <- c(ppvs, ppv)
}
ppvs
g3 <- ggplot(data.frame(x=seq(-2, 2, 0.1), y=round(ppvs*100, 2)), aes(x=x,y=y)) + geom_line(size=2)
g3 <- g3 + theme_Publication()
g3 <- g3 + xlab(expression(bold(Error (epsilon))))
g3 <- g3 + ylab("PPV (%)")
g3 <- g3 + geom_vline(xintercept=1, color="#0072B2", linetype="dashed", size=1.5)
g3 <- g3 + geom_vline(xintercept=0, color="#000000", linetype="dashed", size=1.5)
g3 <- g3 + geom_vline(xintercept=-1, color="#E69F00", linetype="dashed", size=1.5)

g4 <- ggarrange(g1, ggarrange(g2, g3, ncol=2, labels=c("b", "c"), widths = c(1, 0.75)), labels="a", nrow=2)
ggsave("~/data/output/Figure_1.png", g4, units="in", width=7.5, height=5.5, dpi=400)
