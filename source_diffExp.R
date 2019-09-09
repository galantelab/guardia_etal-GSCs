library("DESeq2")

data = read.table("counts_GSCs.txt", header=T, row.names=1)
data = data[1:(dim(data)[1]-5),]

condition = c(rep("MES",3), rep("PN",3))
col_data = as.data.frame(condition)
rownames(col_data) = colnames(data)

dataset = DESeqDataSetFromMatrix(countData = data, colData = col_data, design=~condition)
dataset = dataset[rowSums(counts(dataset))>1, ]

dds = DESeq(dataset)
res = results(dds,contrast=c("condition", "MES", "PN"))

res = as.data.frame(res[order(res$padj),])
res = res[!is.na(res$padj),]
res = res[(res$padj<0.05),]
res = res[(abs(res$log2FoldChange)>1),]
up = res[(res$log2FoldChange>1),]
down = res[(res$log2FoldChange<(-1)),]

write.table(up,file="up_genes.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
write.table(down,file="down_genes.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
