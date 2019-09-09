#Perform survival analyses on TCGA GBM data

library(survival)
library(survminer)

options(scipen = 999)

types = c("MES","PRO")
gene_info = read.table("genes_DE_AS.txt", header=T)
df_mes = data.frame(gene=character(),MES.pvalue=character(),MES.prognosis=character())
df_pro = data.frame(gene=character(),PN.pvalue=character(),PN.prognosis=character())

for (type in types){

  data = read.table(paste(type,"_data.txt",sep=""),header=T,row.names=1)

  for (i in 3:ncol(data)){

    res = coxph(Surv(OS.time, OS)~., data=data[,c(1:2,i)])
    coef = coef(summary(res))[,1]

    categorical_data = data
    categorical_data[,i] = ifelse(categorical_data[,i] > median(categorical_data[,i]),"high", "low")

    tryCatch({
      diff = survdiff(as.formula(paste("Surv(OS.time, OS)", paste(colnames(categorical_data)[i], collapse=" + "), sep=" ~ ")),data=categorical_data);
      p.val = round(1 - pchisq(diff$chisq, length(diff$n) - 1),4)
    }, error = function(e){ p.val = NA})  

    if(as.numeric(p.val)<0.05){

      if(type=="MES" & coef>0){
          df_mes = rbind(df_mes,data.frame(gene=colnames(categorical_data)[i], MES.pvalue=p.val,MES.prognosis="poor"))
      } else if(type=="MES" & coef<0){
          df_mes = rbind(df_mes,data.frame(gene=colnames(categorical_data)[i], MES.pvalue=p.val,MES.prognosis="better"))
      } else if(type=="PRO" & coef>0){
          df_pro = rbind(df_pro,data.frame(gene=colnames(categorical_data)[i], PN.pvalue=p.val,PN.prognosis="poor"))
      } else{
          df_pro = rbind(df_pro,data.frame(gene=colnames(categorical_data)[i], PN.pvalue=p.val,PN.prognosis="better"))
      }

    }

  }

}

res = merge(df_mes,df_pro,by="gene",all=T)
res = merge(gene_info,res,by="gene")
res_mes = res[(!is.na(res$MES.pvalue) & !is.na(res$PN.pvalue)) | (!is.na(res$MES.pvalue) & is.na(res$PN.pvalue) & res$higher.expression=="MES"),]
res_pro = res[(!is.na(res$MES.pvalue) & !is.na(res$PN.pvalue)) | (!is.na(res$PN.pvalue) & is.na(res$MES.pvalue) & res$higher.expression=="PN"),]

write.table(res_mes,"results_univLogRank_MES.txt", row.names=F, col.names=T, quote=F, sep="\t")
write.table(res_pro,"results_univLogRank_PRO.txt", row.names=F, col.names=T, quote=F, sep="\t")

for (type in types){

  data = read.table(paste(type,"_data.txt",sep=""),header=T,row.names=1)
  if(type=="MES"){genes = as.character(res_mes$gene)} else{genes = as.character(res_pro$gene)}
  genes = c("OS", "OS.time", genes)
  genes = paste(genes,collapse="|")
  data = data[,grep(genes,colnames(data))]
  clin = read.delim("clinical.txt", header=T, row.names=1)[,c(1,2,6:8,10:11)]
  clin$age.diagnosis.new <- NA
  clin[(!is.na(clin$age.diagnosis) & clin$age.diagnosis<=55),]$age.diagnosis.new <- "18-55 years"
  clin[(!is.na(clin$age.diagnosis) & clin$age.diagnosis>55),]$age.diagnosis.new <- "56+ years"
  clin$age.diagnosis.new = as.factor(clin$age.diagnosis.new)
  clin = clin[,c(8,2:7)]
  colnames(clin)[1] = "age.diagnosis"
  data = merge(clin,data,by=0)
  rownames(data) = data$Row.names
  data = data[,c(9:10,2:8,11:ncol(data))]

  for(i in 10:ncol(data)){

    covariates = c(colnames(data)[i], "age.diagnosis", "gender", "CIMP.status", "IDH1.status", "MGMT.status", "chr.19.20.co.gain", "chr.7.gain.chr.10.loss")
    res = coxph(formula = formula(paste('Surv(OS.time, OS)~', paste(covariates,collapse=" + "))) , data = data)
    tmp = as.data.frame(coef(summary(res)))
    tmp$variable = rownames(tmp)
    tmp = tmp[,c(6,1,4:5)]
    colnames(tmp) = c("variable", "beta.coef", "Wald.statistic", "Cox.pvalue")
    write.table(tmp, paste("result_",type,"_",colnames(data)[i],".txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")

  }

}

genes = c("CC2D2A","EMILIN2","WBSCR22","LRRFIP1")
for (type in types){
  for (gene in genes){

    myGene = paste("^",gene,"$",sep="")
    data = read.table(paste(type,"_data.txt",sep=""),header=T,row.names=1)
    categorical_data = data[,grep(paste("OS", "OS.time", myGene, sep="|"),colnames(data)),drop=F]
    categorical_data[,3] = ifelse(categorical_data[,3] > median(categorical_data[,3]),"high", "low")
    tryCatch({
      diff = survdiff(as.formula(paste("Surv(OS.time, OS)", paste(gene, collapse=" + "), sep=" ~ ")),data=categorical_data);
      p.val = round(1 - pchisq(diff$chisq, length(diff$n) - 1),4)
    }, error = function(e){ p.val = NA})

    pdf(paste("plot_",type,"_",gene,".pdf",sep=""))
    pp = ggsurvplot(survfit(as.formula(paste("Surv(OS.time, OS)", paste(gene, collapse=" + "), sep=" ~ ")), data=categorical_data), pval=paste("p =",format(p.val,scientific=T),sep=" "),font.legend=20, font.x=22, font.y=22, font.tickslab=18, pval.size=8, pval.coord=c(0,0.05), title= "", legend = c(0.7, 0.9), legend.title="", break.time.by=400, censor=T, legend.labs = c(paste(gene,"high",sep="="),paste(gene,"low",sep="="))) + xlab("Survival time (days)")
    print(pp$plot, newpage = FALSE)
    dev.off()

  }
}
