library("dendextend")

data = read.table("GSCs_psi.txt",header=T)
rownames(data) = paste(as.character(data$ID),as.character(data$GeneID),sep=".")

data = data[,c(3:ncol(data))]
data = as.data.frame(t(data))

sample_distances = dist(data)
clusters = hclust(sample_distances,method="ward.D")
dend = as.dendrogram(clusters)
groupCodes = c(rep("MES",3),rep("PN",3))
colorCodes = c(MES="red", PN="blue")
labels_colors(dend) = colorCodes[groupCodes][order.dendrogram(dend)]

pdf("cluster_psi.pdf")
plot(dend,main="Splicing events",ylab="distance",ylim=c(0,50))
dev.off()

