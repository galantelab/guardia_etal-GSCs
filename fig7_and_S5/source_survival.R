#Perform survival analysis on TCGA GBM data

library(survival)
library(survminer)

options(scipen = 999)

data = read.table("surv_lncRNAs.txt",header=T,row.names=1)
df = data.frame(gene=character(), beta.coef=numeric(),Wald.statistic=numeric(),Cox.pvalue=numeric(),logRank.pvalue=numeric())

for (i in 3:ncol(data)){

  result = coxph(Surv(survival_time, status)~., data=data[,c(1:2,i)])
  result = as.data.frame(coef(summary(result)))
  result$id = rownames(result)
  result = result[,c(6,1,4:5)]

  categorical_data = data
  categorical_data[,i] = ifelse(categorical_data[,i] > median(categorical_data[,i]),"high", "low")

  tryCatch({
    diff = survdiff(as.formula(paste("Surv(survival_time, status)", paste(colnames(categorical_data)[i], collapse=" + "), sep=" ~ ")),data=categorical_data);
    p.val = round(1 - pchisq(diff$chisq, length(diff$n) - 1),4)
  }, error = function(e){ p.val = NA})

  df = rbind(df,data.frame(gene=as.character(result[1]), beta.coef=as.numeric(result[2]),Wald.statistic=as.numeric(result[3]),Cox.pvalue=as.numeric(result[4]), logRank.pvalue=as.numeric(p.val)))

  pdf(paste(colnames(categorical_data)[i],".pdf",sep=""))
  pp = ggsurvplot(survfit(as.formula(paste("Surv(survival_time, status)", paste(colnames(categorical_data)[i], collapse=" + "), sep=" ~ ")), data=categorical_data), pval=paste("p =",format(p.val,scientific=T),sep=" "),font.legend=20, font.x=22, font.y=22, font.tickslab=18, pval.size=8, pval.coord=c(0,0.05), title= "", legend = c(0.7, 0.9), legend.title="", break.time.by=400, censor=T, legend.labs = c(paste(colnames(categorical_data)[i],"high",sep="="),paste(colnames(categorical_data)[i],"low",sep="="))) + xlab("Survival time (days)")
  print(pp$plot, newpage = FALSE)
  dev.off()

}

write.table(df,"results.txt",quote=F,row.names=F,col.names=T,sep="\t")

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

  res = coxph(formula = formula(paste('Surv(survival_time, status)~', paste(covariates,collapse=" + "))) , data = data)
  tmp = as.data.frame(coef(summary(res)))
  tmp$variable = rownames(tmp)
  tmp = tmp[,c(6,1,4:5)]
  colnames(tmp) = c("variable", "beta.coef", "Wald.statistic", "Cox.pvalue")

  write.table(tmp, paste("multiCox_",colnames(data)[i],".txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")

}

