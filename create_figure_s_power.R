library(ggplot2)
library(latex2exp)

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

df <- read.table("./data/figure_s3_table.txt", row.names=1, header=T)
dim(df)
length(seq(-2,2,0.2))
colnames(df) <- c("lfc_sens_testing", seq(-2,2,0.2))
df[is.na(df)] <- T

power <- c()
for(i in seq(-2,2,0.2)) {
	sig_epsilon <- nrow(df[(df[,as.character(i)]),])
	sig_both <- nrow(df[(df[,"lfc_sens_testing"])&(df[,as.character(i)]),])
	power <- c(power, (sig_both/sig_epsilon) * 100)
}
power
power_df <- data.frame(x=seq(-2,2,0.2), y=power)
sum(df[,"lfc_sens_testing"])
g <- ggplot(power_df, aes(x=x,y=y)) + geom_point() + geom_line()
g <- g + ylab("Power (%)") + xlab("ϵ")
g <- g + theme_Publication()
g <- g + ylim(5,15)
ggsave("~/data/output/figure_power.png", g, units="in", width=6, height=4)

