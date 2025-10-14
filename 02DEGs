#求差异lnc并提取表达谱
library(DESeq2)
library(limma)
rm(list = ls())
#求差异表达lnc
counts <- read.csv("F:/RCD/cancer/BLCA/count/TCGA_BLCA_lnc_count.csv",sep = ',',quote ="",row.names=1)
row.names(counts) <- counts$gene_name
counts = 2^counts-1
counts = cbind(gene_name = row.names(counts), counts)
row.names(counts)=NULL
metadata <- data.frame(TCGA_id =colnames(counts)[-1])
table(substring(metadata$TCGA_id,14,15))
sample <- ifelse(substring(metadata$TCGA_id,14,15)=="11","normal","cancer")
metadata$sample <- as.factor(sample)
mycounts <- as.data.frame(avereps(counts[,-1],ID = counts$gene_name) )
dds <-DESeqDataSetFromMatrix(countData=round(mycounts), 
                             colData=metadata, 
                             design=~sample,
                             tidy=FALSE)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
#用Deseq2内置的主成分分析来看一下样本分布
plotPCA(vsd, "sample")
#获取标准化后的数据,这一步会自动过滤掉不符合规定的基因
mRNA_exprSet <- as.data.frame(assay(vsd))
#write.csv(mRNA_exprSet, file = "KICH_pcg_mRNA_exprset_norm.csv",sep=",", row.names = T,quote = F)
#如果需要差异分析
expr_for_diff <- results(dds, tidy=TRUE)
expr_for_diff = na.omit(expr_for_diff)

#单基因差异分析,判断上下调的正反关系
plot_data<-data.frame(mRNA_exprSet["CA9",])
plot_data<-t(plot_data)
plot_data<-data.frame(plot_data)
plot_data <- cbind(plot_data,sample=metadata$sample)
#单基因表达可视化
library(ggplot2)
p<-ggplot(plot_data,aes(x=sample,y=CA9,fill=sample))+
  geom_boxplot()+
  theme_classic()+
  ggpubr::stat_compare_means(color="red")+
  ggsci::scale_fill_jco()+
  theme(legend.position = "none")+
  ylab("CA9 counts")
print(p)

#基因上下调
logFC_cutoff <- 2
k1 = (expr_for_diff$padj < 0.05)&(expr_for_diff$log2FoldChange < -logFC_cutoff)
k2 = (expr_for_diff$padj < 0.05)&(expr_for_diff$log2FoldChange > logFC_cutoff)
#expr_for_diff$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
expr_for_diff$change = ifelse(k1,"UP",ifelse(k2,"DOWN","NOT"))#FC正负与上下调相反
table(expr_for_diff$change)
head(expr_for_diff)
#del <- which(expr_for_diff$change=="NOT")
#data <- expr_for_diff[-del,]
expr_for_diff$log2FoldChange <- -expr_for_diff$log2FoldChange
write.csv(expr_for_diff,file = "LUSC_geneset_change.csv")
