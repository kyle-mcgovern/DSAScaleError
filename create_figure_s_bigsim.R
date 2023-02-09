library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(extrafont)
library(latex2exp)
loadfonts()

# True positive (\epsilon=0) p-value <=0.05 (\epsilon=X) p-value <=0.05, ES same
# False positive (\epsilon=0) p-value <=0.05
# cols: 1=dist_type, 2=ds, 3=d, 4=epsilon, 5=index, 6=p-value, 7=es
create_ppv_matrix <- function(m) {
	epsilons <- unique(m[,4])
	D <- unique(m[,3])
	Ds <- unique(m[,2])
	dist_types <- unique(m[,1])
	ppvs <- c()
	for(d in D) {
		for(ds in Ds) {
			for(dist_type in dist_types) {
				cm <- m[(m[,3]==d)&(m[,2]==ds)&(m[,1]==dist_type),]
				assumed_pos_up_cm <- cm[(cm[,4]==0)&(cm[,6]<=0.05)&(cm[,7]>0),]
				assumed_pos_down_cm <- cm[(cm[,4]==0)&(cm[,6]<=0.05)&(cm[,7]<0),]
				assumed_all_cm <- rbind(assumed_pos_up_cm, assumed_pos_down_cm)
				for(epsilon in epsilons) {
					true_pos <- nrow(cm[(cm[,4]==epsilon)&(cm[,6]<=0.05)&(cm[,7]>0)&(cm[,5]%in%assumed_pos_up_cm[,5]),])
					true_pos <- true_pos + nrow(cm[(cm[,4]==epsilon)&(cm[,6]<=0.05)&(cm[,7]<0)&(cm[,5]%in%assumed_pos_down_cm[,5]),])
					false_pos <- nrow(cm[(cm[,4]==epsilon)&(cm[,6]>0.05)&(cm[,5]%in%assumed_all_cm[,5]),])
					ppv <- (true_pos / (true_pos+false_pos)) * 100
					ppvs <- rbind(ppvs, c(d, ds, dist_type, epsilon, ppv))
				}
			}
		}
	}
	colnames(ppvs) <- c("D", "Ds", "dist_type", "epsilon", "ppv")
	ppvs
}

#m <- read.table("./data/figure_s1_data.txt")
#head(m)
#ppv_mat <- create_ppv_matrix(m)
#write.table(ppv_mat, "ppv_mat.txt")

ppv_mat <- read.table("ppv_mat.txt")
ppv_mat

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
               panel.grid.major = element_blank(),
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

create_fig_2_line <- function(data) {
	data$dist_type <- factor(data$dist_type, levels=c("norm", "uniform", "left_skew", "right_skew"))
	g <- ggplot(data, aes(x=epsilon, y=ppv, color=dist_type, fill=dist_type, shape=dist_type)) + geom_line() + geom_point(size=2, color="black")
	g <- g + theme_Publication()
	print(head(data))
	g <- g + facet_grid(facets=D~Ds, labeller=label_bquote(cols=D[S]:.(Ds), rows=D:.(D)))
	g <- g + theme(strip.text.x=element_text(family="Times", face="bold", size=13), strip.text.y=element_text(family="Times", face="bold", size=13))
	g <- g + xlab(TeX("$ϵ^\\perp$"))
	g <- g + theme(legend.title=element_blank())
	g <- g + theme(panel.grid.major = element_line(color="#dfdfdf", size=0.5))
	g <- g + scale_shape_manual(values=c(21, 22, 23, 24), labels=c("Normal", "Uniform", "Left-skew", "Right-skew"))
	g <- g + scale_fill_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"), labels=c("Normal", "Uniform", "Left-skew", "Right-skew"))
	g <- g + scale_color_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"), labels=c("Normal", "Uniform", "Left-skew", "Right-skew"))
	g <- g + scale_x_continuous(breaks=seq(-2,2,1))
	g <- g + ylab("Positive Predictive Value (PPV)")
	return(g)

}

create_fig_2_samples <- function(RDS_path, shape, color, title) {
	curr_data <- data.frame(sort(readRDS(RDS_path), decreasing=T))
	curr_data$gene_index <- 1:nrow(curr_data)
	colnames(curr_data) <- c("lfc", "glr")
	g <- ggplot(curr_data, aes(x=glr, y=lfc)) + geom_point(shape=shape, fill=color, color="black")
	g <- g + theme_Publication()
	g <- g + theme(plot.title=element_text(size=12, face="bold"))
	g <- g + theme(axis.title=element_text(size=8))   
	g <- g + theme(axis.text=element_text(size=7))
	g <- g + scale_x_continuous(breaks=c(1,250,500))
	g <- g + xlab("Entity List Rank")
	g <- g + ylab("Log Fold Change")
	g <- g + theme(panel.grid.major = element_line(color="#dfdfdf", linetype=1))
	g <- g + ggtitle(title)
	g
}

g1 <- create_fig_2_line(ppv_mat)
g2 <- create_fig_2_samples("150_500_norm.RDS", 21, "#a6cee3", "Normal")
g3 <- create_fig_2_samples("150_500_uniform.RDS", 22, "#1f78b4", "Uniform")
g4 <- create_fig_2_samples("150_500_left_skew.RDS", 23, "#b2df8a", "Left-skew")
g5 <- create_fig_2_samples("150_500_right_skew.RDS", 24, "#33a02c", "Right-skew")

#ggsave("~/data/output/figure_s1.png", g)
g6 <- ggarrange(ggarrange(g2, g3, g4, g5, ncol=4), g1, labels=c("a", "b"), nrow=2, heights=c(0.8,2))
ggsave("~/data/output/s_larger_sim.png", g6, units="in", width=9, height=8, dpi=400)


