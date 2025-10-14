###GSEA判断通路是否在癌症中差异表达###
#1.首先计算基因的差异倍数，以差异倍数作为排序依据
library(DESeq2)
library(limma)
counts <- read.csv("F:/RCD/cancer/COAD/count/TCGA_COAD_pcg_count.csv",sep = ',',quote ="",row.names=1)
counts = 2^counts-1#deseq2输入的矩阵为未log、整数
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
plot_data<-data.frame(mRNA_exprSet["TSPAN6",])
plot_data<-t(plot_data)
plot_data<-data.frame(plot_data)
plot_data <- cbind(plot_data,sample=metadata$sample)
#单基因表达可视化
library(ggplot2)
p<-ggplot(plot_data,aes(x=sample,y=	TSPAN6,fill=sample))+
  geom_boxplot()+
  theme_classic()+
  ggpubr::stat_compare_means(color="red")+
  ggsci::scale_fill_jco()+
  theme(legend.position = "none")+
  ylab("TSPAN6 counts")
print(p)

expr_for_diff$log2FoldChange <- -expr_for_diff$log2FoldChange
rownames(expr_for_diff) <- expr_for_diff$row
expr_for_diff <- expr_for_diff[,-1]
FCgenelist <- expr_for_diff$log2FoldChange
names(FCgenelist) <- row.names(expr_for_diff) 
FCgenelist <- sort(FCgenelist, decreasing = T)
#2.GSEA富集
library(fgsea)
pathway <- readRDS("C:/Users/Dell/Desktop/ImmuLnc/pathway.Rdata")
set.seed(123)#FGSEA函数记得都要设置seed
fgseaRes <- fgsea(pathways = pathway, 
                  stats = FCgenelist,
                  minSize=1, 
                  maxSize=5000, 
                  nperm=1000)
fgseaRes1 <- fgseaRes[which( fgseaRes$pval<0.05&fgseaRes$padj<0.25),]
write.csv(fgseaRes1,"C:/Users/Dell/Desktop/pathway/COAD_GSEA.csv")
#3.富集结果可视化
library(GseaVis)
gseaNb(object = fgseaRes1,
       geneSetID = 'GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS')
###单细胞GSEA###
result <- FindAllMarkers(data,only.pos =T, min.pct = 0.5, logfc.threshold = 0.5)
write.csv(result,"C:/Users/Dell/Desktop/26.大类及上皮细胞亚群rds/Main/deg.csv")
result <- read.csv("C:/Users/Dell/Desktop/26.大类及上皮细胞亚群rds/Main/deg.csv")
TandNK  <- result[which(result$cluster== "TandNK"),]
write.csv(TandNK,"C:/Users/Dell/Desktop/26.大类及上皮细胞亚群rds/EpithelialCells/TandNK_deg.csv")
###小鼠基因GSEA###
Apoptosis <- read.csv("F:/RCD/cancer/####RCD-pathway/Apoptosis.csv")
Apoptosis <- homologene(Apoptosis$Apoptosis, inTax = 9606, outTax = 10116)
Apoptosis <- as.data.frame(Apoptosis[2])
names(Apoptosis) <- "gene"
Autophage <- read.csv("F:/RCD/cancer/####RCD-pathway/Autophagy.csv")
Autophage <- homologene(Autophage$Autophagy, inTax = 9606, outTax = 10116)
Autophage <- as.data.frame(Autophage[2])
names(Autophage) <- "gene"
Cuproptosis <- read.csv("F:/RCD/cancer/####RCD-pathway/Cuproptosis.csv")
Cuproptosis <- homologene(Cuproptosis$Cuproptosis, inTax = 9606, outTax = 10116)
Cuproptosis <- as.data.frame(Cuproptosis[2])
names(Cuproptosis) <- "gene"
Ferroptosis <- read.csv("F:/RCD/cancer/####RCD-pathway/Ferroptosis.csv")
Ferroptosis <- homologene(Ferroptosis$Ferroptosis, inTax = 9606, outTax = 10116)
Ferroptosis <- as.data.frame(Ferroptosis[2])
names(Ferroptosis) <- "gene"
Necroptosis <- read.csv("F:/RCD/cancer/####RCD-pathway/Necroptosis.csv")
Necroptosis <- homologene(Necroptosis$Necroptosis, inTax = 9606, outTax = 10116)
Necroptosis <- as.data.frame(Necroptosis[2])
names(Necroptosis) <- "gene"
Pyroptosis <- read.csv("F:/RCD/cancer/####RCD-pathway/Pyroptosis.csv")
Pyroptosis <- homologene(Pyroptosis$Pyroptosis, inTax = 9606, outTax = 10116)
Pyroptosis <- as.data.frame(Pyroptosis[2])
names(Pyroptosis) <- "gene"
pathway <- list(Apoptosis$gene,Autophage$gene,Cuproptosis$gene,Ferroptosis$gene,Necroptosis$gene,Pyroptosis$gene)
names(pathway) <- c("Apoptosis","Autophage","Cuproptosis","Ferroptosis","Necroptosis","Pyroptosis")

geneList= TandNK$avg_log2FC
names(geneList)= toupper(TandNK$X)
geneList=sort(geneList,decreasing = T)
head(geneList)
library(fgsea)
set.seed(123)#FGSEA函数记得都要设置seed
fgseaRes <- fgsea(pathways = pathway, 
                  stats = geneList,
                  minSize=1, 
                  maxSize=5000, 
                  nperm=1000)
fgseaRes1 <- fgseaRes[which( fgseaRes$pval<0.05&fgseaRes$padj<0.25),]
