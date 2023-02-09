library(ggplot2)
library(reshape2)
library(ggthemes)
library(latex2exp)
library(ggpubr)

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

interleave <- function(x,y){
  lx <- length(x)
  ly <- length(y)
  n <- max(lx,ly)
  as.vector(rbind(rep(x, length.out=n), rep(y, length.out=n)))
}

GSEA_constant <- read.table("./data/GSEA_LFC_C2_constant.txt")
CAMERA_constant <- read.table("./data/CAMERA_LFC_C2_constant.txt")
GSEA_constant
colnames(GSEA_constant) <- seq(-0.25, 0.25, 0.05)
colnames(CAMERA_constant) <- seq(-0.25, 0.25, 0.05)

ppv_camera <- apply(CAMERA_constant, 2, function(col) sum(col<=0.05&CAMERA_constant[,6]<=0.05)/sum(CAMERA_constant[,6]<=0.05))
ppv_gsea <- apply(GSEA_constant, 2, function(col) sum(col<=0.05&GSEA_constant[,6]<=0.05)/sum(GSEA_constant[,6]<=0.05))
constant_df <- cbind(ppv_camera, ppv_gsea, 1) * 100
colnames(constant_df) <- c("CAMERA", "GSEA-LFC-S", "GSEA-CW-S")
constant_df

all_df <- data.frame()
zeros_gsea <- read.table("./data/GSEA_LFC_C2_0_01.txt")
zeros_camera <- read.table("./data/CAMERA_LFC_C2_0_01.txt")
for(i in seq(0, 0.25, 0.025)) {
	gsea_df <- read.table(paste0("./data/GSEA_LFC_C2_", i, "_01.txt"))
	gsea_ppvs <- apply(gsea_df, 2, function(col) sum(col<=0.05&zeros_gsea<=0.05)/sum(zeros_gsea<=0.05))
	camera_df <- read.table(paste0("./data/CAMERA_LFC_C2_", i, "_01.txt"))
	camera_ppvs <- apply(camera_df, 2, function(col) sum(col<=0.05&zeros_camera<=0.05)/sum(zeros_camera<=0.05))
	all_df <- rbind(all_df, data.frame(epsilon=i, ppv=gsea_ppvs*100, method="GSEA-LFC-S"))
	all_df <- rbind(all_df, data.frame(epsilon=i, ppv=camera_ppvs*100, method="CAMERA"))
	all_df <- rbind(all_df, data.frame(epsilon=i, ppv=100, method="GSEA-CW-S"))
}

g1 <- ggplot(melt(constant_df), aes(x=Var1, y=value, color=Var2, shape=Var2)) + geom_point(size=2, stroke=1.5)
g1 <- g1 + scale_shape_manual(values=c(3, 4, 5))
g1 <- g1 + theme_bw()
g1 <- g1 + theme(legend.title=element_blank())
g1 <- g1+ xlab(TeX("Error in $\\hat{\\theta}^{\\perp}$ ($δ^\\perp$)"))
g1 <- g1 + ylab("Positive Predictive Value (%)")
g1 <- g1 + ylim(50, 100)
g1 <- g1 + scale_x_continuous(breaks=seq(-.25, .25, 0.05), labels=c("", -0.2, "", -0.1, "", 0, "", 0.1, "", 0.2, ""))
g1 <- g1 + theme(legend.position="bottom")
g1 <- g1 + scale_colour_manual(values=c("#000000", "#E69F00", "#0072B2"))
g1 <- g1 + theme(legend.margin=margin(0,0,0,0), legend.box.margin=margin(-8,-8,-8,-8))

g2 <- ggplot(all_df, aes(x=as.factor(epsilon), y=ppv, color=factor(method, levels=c("CAMERA", "GSEA-LFC-S", "GSEA-CW-S")))) + geom_boxplot(outlier.size=0.25)
g2 <- g2 + theme_bw()
g2 <- g2 + scale_colour_manual(values=c("#000000", "#E69F00", "#0072B2"))
g2 <- g2 + theme(legend.title=element_blank())
g2 <- g2 + ylab("Positive Predictive Value (%)")
g2 <- g2 + xlab(TeX("Error in $\\hat{W}^{\\perp}$ ($𝛾^\\perp$)"))
g2 <- g2 + ylim(50, 100)
labels <- c("0", "", "0.05", "", "0.1", "", "0.15", "", "0.2", "", "0.25")
g2 <- g2 + scale_x_discrete(breaks=seq(0, 0.25, 0.025), labels=labels) + theme(legend.position="bottom")
g2 <- g2 + theme(legend.margin=margin(0,0,0,0), legend.box.margin=margin(-8,-8,-8,-8))

g <- ggarrange(g1, g2, labels=c("a", "b"), ncol=2, nrow=1)
g
ggsave("~/data/output/Scale_Sensitivity_Analysis.png", g, units="in", width=7.25, height=3.25, dpi=400)

