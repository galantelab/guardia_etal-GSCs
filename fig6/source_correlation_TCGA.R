#Calculate Spearman correlation among selected RBPs on TCGA GBM data and create correlation plots #########################################

library(corrplot)
library(Hmisc)

mes = read.table("TCGA_GBM_MES.txt", header=T, row.names=1, check.names=F)
pro = read.table("TCGA_GBM_PN.txt", header=T, row.names=1, check.names=F)

group_mes1 = c("SAP18", "CLK4", "SMG6", "PPARGC1A")
mes1 = mes[,as.vector(na.omit(pmatch(group_mes1,colnames(mes))))]

group_mes2 = c("SMN2", "RBM11", "SLU7", "ZNF326", "HNRNPA2B1", "PTBP2", "SFPQ", "TIA1", "SF3B6", "THOC1", "CCNT2")
mes2 = mes[,as.vector(na.omit(pmatch(group_mes2,colnames(mes))))]

group_pro1 = c("SCAF4", "SFPQ", "TIA1", "THOC1", "GPATCH1", "LUC7L", "CLK1")
pro1 = pro[,as.vector(na.omit(pmatch(group_pro1,colnames(pro))))]

group_pro2 = c("UPF3A", "SLU7", "ZNF326", "U2AF1", "SON", "SRSF12", "PSRPK2", "HNRNPA2B1", "PTBP3", "QKI", "CCNT2")
pro2 = pro[,as.vector(na.omit(pmatch(group_pro2,colnames(pro))))]

#Calculate pairwise correlations
mes1_cor = rcorr(as.matrix(mes1),type="spearman")
mes2_cor = rcorr(as.matrix(mes2),type="spearman")

pro1_cor = rcorr(as.matrix(pro1),type="spearman")
pro2_cor = rcorr(as.matrix(pro2),type="spearman")

#Make correlation plots - MES
pdf("plot_MES1.pdf")
corrplot(mes1_cor$r, p.mat = mes1_cor$P, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, method = "circle", tl.cex = 0.8, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()

pdf("plot_MES2.pdf")
corrplot(mes2_cor$r, p.mat = mes2_cor$P, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, method = "circle", tl.cex = 0.8, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()

#Make correlation plots - PRO
pdf("plot_PRO1.pdf")
corrplot(pro1_cor$r, p.mat = pro1_cor$P, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, method = "circle", tl.cex = 0.8, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()

pdf("plot_PRO2.pdf")
corrplot(pro2_cor$r, p.mat = pro2_cor$P, insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = .9, method = "circle", tl.cex = 0.8, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()
