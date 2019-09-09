library(corrplot)

mes = read.table("MES_GBM_data.txt",header=T,row.names=1)
pro = read.table("PRO_GBM_data.txt",header=T,row.names=1)

group_mes1 = c("SAP18", "CLK4", "SMG6", "PPARGC1A")
mes1 = mes[,as.vector(na.omit(pmatch(group_mes1,colnames(mes))))]
group_mes2 = c("HNRNPA2B1", "PTBP2", "SFPQ", "TIA1", "THOC1", "CCNT2")
mes2 = mes[,as.vector(na.omit(pmatch(group_mes2,colnames(mes))))]

group_pro1 = c("SCAF4", "SFPQ", "TIA1", "THOC1", "GPATCH1", "LUC7L")
pro1 = pro[,as.vector(na.omit(pmatch(group_pro1,colnames(pro))))]
group_pro2 = c("UPF3A", "U2AF1", "SON", "HNRNPA2B1", "PTBP3", "QKI", "CCNT2")
pro2 = pro[,as.vector(na.omit(pmatch(group_pro2,colnames(pro))))]

mes1_cor = cor(mes1,method="spearman")
mes2_cor = cor(mes2,method="spearman")
pro1_cor = cor(pro1,method="spearman")
pro2_cor = cor(pro2,method="spearman")

pdf("plot_MES_GBM_1.pdf")
corrplot(mes1_cor, method = "circle", tl.cex = 0.8, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()

pdf("plot_MES_GBM_2.pdf")
corrplot(mes2_cor, method = "circle", tl.cex = 0.8, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()

pdf("plot_PRO_GBM_1.pdf")
corrplot(pro1_cor, method = "circle", tl.cex = 0.8, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()

pdf("plot_PRO_GBM_2.pdf")
corrplot(pro2_cor, method = "circle", tl.cex = 0.8, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()
