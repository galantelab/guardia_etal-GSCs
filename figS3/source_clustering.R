library("dendextend")

data = read.table("diff_exp_lncRNAs.txt", header=T)
sample_distances = dist(t(data))
clusters = hclust(sample_distances,method="ward.D")
dend = as.dendrogram(clusters)
groupCodes = c(rep("MES",3),rep("PN",3))
colorCodes = c(MES="red", PN="blue")
labels_colors(dend) = colorCodes[groupCodes][order.dendrogram(dend)]
pdf("cluster_lncRNA.pdf") 
plot(dend,main="lncRNAs expression",ylab="distance",ylim=c(0,400))
dev.off()
