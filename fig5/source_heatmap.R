#Plot heatmap of normalized values (ordered by up/down and gene symbol) ####################################

#Load libraries
library(gplots)

#Read DE RBPs
up = read.table("up_multi-exon_rbps.txt",header=T,row.names=1)
down = read.table("down_multi-exon_rbps.txt",header=T,row.names=1)

#Read normalized expression data
data = read.table("normalized_counts.txt")

#Select only up/down-regulated RBPs from normalized expression data
up_data = merge(up,data,by=0)
down_data = merge(down,data,by=0)

#Reorder dataframes by gene names
up_data = up_data[order(up_data$name),]
down_data = down_data[order(down_data$name),]

#Merge dataframes and process them
data = rbind(down_data, up_data)
rownames(data) = data$name
data = data[,9:ncol(data)]
colnames(data) = c("MES-83","MES-326","MES-1123","PN-19","PN-157","PN-528")

pdf("heatmap.pdf")
heatmap.2(as.matrix(data),density.info="none",trace="none",keysize=0.8,col=redblue(80),scale="row",dendrogram="none",cexRow=0.5,key.title = "",key.xlab="",srtCol=0, Rowv=FALSE,Colv=FALSE,labCol="")
dev.off()
