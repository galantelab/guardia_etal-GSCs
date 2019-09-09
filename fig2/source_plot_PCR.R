#Make barplots for PCR results ############################################

library("ggplot2")

data = read.table("MES_events.txt", header=T)
data$sample <- factor(data$sample, levels = as.character(data$sample))

pdf("MES_plot.pdf", width=7, height=7)
ggplot(data, aes(x=gene, y=fold.change, fill=sample)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=fold.change, ymax=fold.change+standard.deviation), width=.2, position=position_dodge(.9)) + theme(axis.title.x = element_blank(), legend.title=element_blank()) + ylab("Fold change") + scale_y_continuous(limits=c(0,5), expand = c(0,0))
dev.off()

data = read.table("PRO_events.txt", header=T)
data$sample <- factor(data$sample, levels = as.character(data$sample))

#Make plot
pdf("PRO_plot.pdf", width=9, height=7)
ggplot(data, aes(x=gene, y=fold.change, fill=sample)) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=fold.change, ymax=fold.change+standard.deviation), width=.2, position=position_dodge(.9)) + theme(axis.title.x = element_blank(), legend.title=element_blank()) + ylab("Fold change") + scale_y_continuous(limits=c(0,1.5), expand = c(0,0))
dev.off()
