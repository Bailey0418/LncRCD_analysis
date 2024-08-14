###计算每个样本的IC50值###
#准备表达谱数据
setwd("C:/Users/Dell/Desktop/drug/UCEC")
e <- read.csv("F:/RCD/cancer/UCEC/rcdlnc/sig_pairs.csv")
lnc <- unique(e$lncRNA)
data <- read.csv("F:/RCD/cancer/UCEC/fpkm/TCGA_UCEC_lncRNA_tum_fpkm.csv",row.names = 1)
data <- data[lnc,]
#计算药物敏感性
CTRP2_Expr <- readRDS("F:/RCD/课题/数据备份/分析/oncopredict/CTRP2_Expr (TPM, not log transformed).rds")
CTRP2_Res <- readRDS("F:/RCD/课题/数据备份/分析/oncopredict/CTRP2_Res.rds")
library(oncoPredict)
calcPhenotype(trainingExprData = CTRP2_Expr,
              trainingPtype = CTRP2_Res,
              testExprData = as.matrix(data),#需要matrix
              batchCorrect = 'eb',  
              #IC50是对数转换的，所以表达矩阵也用对数转换过的
              powerTransformPhenotype = F,
              minNumSamples = 20,
              printOutput = T,
              removeLowVaryingGenes = 0.2,
              removeLowVaringGenesFrom = "homogenizeData")
###构建lncRNA——药物敏感性关系###
setwd("C:/Users/Dell/Desktop/drug/UCEC")
e <- read.csv("F:/RCD/cancer/UCEC/rcdlnc/sig_pairs.csv")
lnc <- unique(e$lncRNA)
data <- read.csv("F:/RCD/cancer/UCEC/fpkm/TCGA_UCEC_lncRNA_tum_fpkm.csv",row.names = 1)
sample_Riskscore <- data[lnc,]

sample_HighLow <- as.data.frame(apply(sample_Riskscore, 1, function(x) {
  ifelse(x > median(x), "High", "Low")
}))

sample_drug_data <- read.csv("calcPhenotype_Output/DrugPredictions.csv",row.names = 1)
sample_Riskscore <- as.data.frame(t(sample_Riskscore))
library(doParallel)
library(foreach)
library(dplyr)
library(psych)
library(tidyverse)
ceRNA_names <- colnames(sample_HighLow)
DE_drug_data <- data.frame()

# 找到核心个数# 构建核心类，并不是越多越好
cl <- makeCluster(detectCores() - 1)
# 调用核心
registerDoParallel(detectCores())

for(i in 1:length(ceRNA_names)){
  #2.1 在每个ceRNA中将高低风险组样本根据进行差异药物反应分析 DE
  sample_drug_data_temp <- cbind(group=sample_HighLow[rownames(sample_drug_data),ceRNA_names[i]],
                                 sample_drug_data)
  drug_names <- colnames(sample_drug_data_temp)[-1]
  
  # 由于第一组是分组变量，# foreach（for循环） - 耗时0.03 s
  #FC
  group_mean <- aggregate(sample_drug_data_temp[,-1], list(sample_drug_data_temp$group), mean)
  rownames(group_mean) <- as.character(group_mean[,1])
  group_mean <- group_mean[,-1]
  
  log2FC <- log2(as.numeric(group_mean["High",])/as.numeric(group_mean["Low",]))
  
  #检验#https://www.jianshu.com/p/938c6ee56dd9
  pvalue <- foreach(j = 2 : ncol(sample_drug_data_temp), .combine = "c") %dopar% {wilcox.test( sample_drug_data_temp[,j] ~ group, data = sample_drug_data_temp)$p.value}
  pvalue <- cbind(drug_names, cbind(log2FC=as.numeric(t(log2FC)),pvalue))
  
  DE_drug_data <- rbind(DE_drug_data,
                        cbind(data.frame(ceRNA=ceRNA_names[i],method="Wilcoxon rank-sum test"),
                              pvalue))
  
}
# 退出
stopImplicitCluster()

#2.2 在每个ceRNA中,riskcore——样本药物反应进行斯皮尔曼相关
corr_ce_drug <- corr.test(sample_Riskscore,sample_drug_data,use="pairwise",method="spearman",adjust = "fdr",ci=FALSE)
corr_ce_drug_r <- as.data.frame(cbind(ceRNA=rownames(corr_ce_drug$r),corr_ce_drug$r))
corr_ce_drug_r <- gather(corr_ce_drug_r, key = "drug",value = "scc_r",-`ceRNA`)
corr_ce_drug_p <- as.data.frame(cbind(ceRNA=rownames(corr_ce_drug$p),corr_ce_drug$p))
corr_ce_drug_p <- gather(corr_ce_drug_p, key = "drug",value = "pvalue",-`ceRNA`)
corr_ce_drug_fdr <- as.data.frame(cbind(ceRNA=rownames(corr_ce_drug$p.adj),corr_ce_drug$p.adj))
corr_ce_drug_fdr <- gather(corr_ce_drug_fdr, key = "drug",value = "fdr",-`ceRNA`)
corr_ce_drug_r$pvalue <- corr_ce_drug_p[,3]
corr_ce_drug_r$fdr <- corr_ce_drug_fdr[,3]
#汇总
SCC_drug_data <- cbind(data.frame(ceRNA=corr_ce_drug_r[,1],method="spearman"),
                       corr_ce_drug_r[,2:5])
DE_drug_data_temp <- DE_drug_data
SCC_drug_data_temp <- SCC_drug_data
DE_drug_data_temp <- cbind(DE_drug_data_temp[,1:3],as.data.frame(lapply(DE_drug_data_temp[,4:5],as.numeric)))
SCC_drug_data_temp <- cbind(SCC_drug_data_temp[,1:3],as.data.frame(lapply(SCC_drug_data_temp[,4:6],as.numeric)))
#合并 cancer ceRNA drug train_database method1 log2FC pvalue method2 scc_r pvalue fdr
colnames(SCC_drug_data_temp) <- c("lncRNA","method_scc","drug","scc_r","scc_pvalue","fdr")
colnames(DE_drug_data_temp) <- c("lncRNA","method_DE","drug","log2FC","DE_pvalue")
cancer_drug_data <- right_join(DE_drug_data_temp, SCC_drug_data_temp, by = c("lncRNA","drug"))  
#cancer_drug_data <- cancer_drug_data[,c(1:2,5,3:4,6:11)]
#write.csv(cancer_drug_data,"KICH_drugsensitity_data.csv")
#卡阈值#设置阈值 DE(log2FC>0.1) Scc(r<-0.3)
cancer_drug_data_sig <- filter(cancer_drug_data,abs(log2FC) > 0.25 & scc_r < (- 0.3) & DE_pvalue <0.05 & scc_pvalue <0.05)
table(cancer_drug_data_sig$drug)
write.csv(cancer_drug_data_sig,"drugsensitity_sigdata.csv")
###绘制药物在lncRNA高低分组敏感性箱线图###
library(ggplot2)#绘图包
library(ggpubr)#基于ggplot2的可视化包，主要用于绘制符合出版要求的图形
library(ggsignif)#用于P值计算和显著性标记
library(tidyverse)
library(ggsci)
d1 <- as.data.frame(cbind(sample_HighLow[,c("RP1-63M2.7")],sample_drug_data[,c("PF.3758309")]))
names(d1) <- c("group","IC50")
d1$IC50 <- as.numeric(d1$IC50)

p <- ggplot(d1, aes(group, IC50))+
  geom_boxplot(aes(fill = group ))+
  geom_signif(comparisons = list(c("High", "Low")), test = wilcox.test, step_increase = 0.2)+
  theme_bw() + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  scale_fill_npg()
p
