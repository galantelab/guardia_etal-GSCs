#Make barplots for RNA-Seq results (with statistical test for differences) ############################################

library("ggplot2")
library(ggpubr)

data = read.table("MES_events.txt", header=T)
data$sample <- factor(data$sample, levels = as.character(data$sample))

pdf("MES_plot.pdf", width=7, height=7)
ggplot(data, aes(x=gene, y=delta.psi, fill=sample)) + geom_bar(stat="identity", color="black", position=position_dodge()) + theme(axis.title.x = element_blank(), legend.title=element_blank()) + ylab("Delta PSI") + scale_y_continuous(limits=c(0,1.1), expand = c(0,0)) + stat_compare_means(aes(group = group), label="p.format",label.y=1.05, method="wilcox.test")
dev.off()

data = read.table("PRO_events.txt", header=T)
data$sample <- factor(data$sample, levels = as.character(data$sample))

pdf("PRO_plot.pdf", width=9, height=7)
ggplot(data, aes(x=gene, y=delta.psi, fill=sample)) + geom_bar(stat="identity", color="black", position=position_dodge()) + theme(axis.title.x = element_blank(), legend.title=element_blank()) + ylab("Delta PSI") + scale_y_continuous(limits=c(0,1.1), expand = c(0,0)) + stat_compare_means(aes(group = group), label="p.format",label.y=1.05, method="wilcox.test")
dev.off()
