library(corrplot)

mes_data = read.table("MES_GSCs_data.txt",header=T)
pro_data = read.table("PRO_GSCs_data.txt",header=T)

df_mes = data.frame(gene1=character(), gene2=character(), coefficient=numeric(), p.value=numeric())
df_pro = data.frame(gene1=character(), gene2=character(), coefficient=numeric(), p.value=numeric())

for(i in 1:ncol(mes_data)){
  for(j in i:ncol(mes_data)){

    res = suppressWarnings(cor.test(mes_data[,i], mes_data[,j], alternative="two.sided", method="spearman", exact=F))
    df_mes = rbind(df_mes, data.frame(gene1=colnames(mes_data)[i], gene2=colnames(mes_data)[j], coefficient=as.numeric(res$estimate), p.value=res$p.value, stringsAsFactors=F))

  }
}

for(i in 1:ncol(pro_data)){
  for(j in i:ncol(pro_data)){

    res = suppressWarnings(cor.test(pro_data[,i], pro_data[,j], alternative="two.sided", method="spearman", exact=F))
    df_pro = rbind(df_pro, data.frame(gene1=colnames(pro_data)[i], gene2=colnames(pro_data)[j], coefficient=as.numeric(res$estimate), p.value=res$p.value,stringsAsFactors=F))

  }
}

df_mes = df_mes[!is.na(df_mes$coefficient) & (df_mes$coefficient>0.8 | df_mes$coefficient<(-0.8)),]
df_pro = df_pro[!is.na(df_pro$coefficient) & (df_pro$coefficient>0.8 | df_pro$coefficient<(-0.8)),]

write.table(df_mes,"corr_matrix_MES_GSCs.txt",row.names=F,col.names=T,quote=F,sep="\t")
write.table(df_pro,"corr_matrix_PRO_GSCs.txt",row.names=F,col.names=T,quote=F,sep="\t")

mes_cor = suppressWarnings(cor(mes_data,method="spearman"))
pdf("plot_MES_GSCs.pdf")
corrplot(mes_cor, method = "circle", tl.cex = 0.5, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()

pro_cor = suppressWarnings(cor(pro_data,method="spearman"))
pro_cor[,49] = 0
pro_cor[49,] = 0
pro_cor[49,49] = 1

pdf("plot_PRO_GSCs.pdf")
corrplot(pro_cor, method = "circle", tl.cex = 0.5, tl.col = "black", order = "hclust", type = "lower",tl.srt = 45)
dev.off()
