setwd("D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC亚型")
#####KIRC immune subtype
lnc <- read.csv("D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC预后模型构建/KIRC_单因素COX分析结果.csv")
data <- read.csv("D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC预后模型构建/KIRC_cox_data.csv",row.names = 1)
expr <- data[,-c(1:3)]
expr <- t(expr)
expr <- as.data.frame(expr)
lnc$gene <- gsub("-", ".", lnc$gene)
expr <- expr[rownames(expr) %in% lnc$gene,]

library(NMF)
ranks <- 2:8
f = "rank.rdata"
if(!file.exists(f)){
  result = nmf(expr,ranks,nrun = 150,seed=123)
  save(result,file = f)
}
load(f)
plot(result)

rank = 4
result2 <- nmf(expr,rank = rank,nrun = 150)
index <- extractFeatures(result2,"max") # 能给出选了哪些关键的基因，可以用这些基因再次聚类
fpkm = expr[unlist(index),]
fpkm <- na.omit(fpkm)
result2 <- nmf(fpkm,rank = rank,nrun = 150)
group <- predict(result2) # 提出亚型
table(group)
my_palette <- c("#98d09d","#fbf398","#e77381","#9b8191","#7EABCA")
annColors <- list(group = setNames(my_palette, unique(group)))
basismap(dp,
         cexCol = 1.5,
         cexRow = 1,
         annColors = annColors)
dp = fpkm[,order(group)]
write.csv(group,"KIRC_group.csv")

#group 生存
library(survival)
library(survminer)
library(tidyverse)
common_samples <- intersect(rownames(data), names(group))
data <- data[common_samples, ]
group <- group[common_samples]
all(rownames(data) == names(group))
data$group = group
data$OS <- as.numeric(data$OS)
sfit <- survfit(Surv(OS.time, OS) ~ group, data = data)
ggsurvplot(sfit, data = data, pval = TRUE, palette = my_palette)

#group PCA
library(Rtsne)
expr <- t(expr)#矩阵应为样本数*基因数
tsne_out = Rtsne(expr,perplexity = 30)
all(rownames(expr) == names(group))
pdat = data.frame(tsne_out$Y,factor(group))
colnames(pdat) = c("Y1","Y2","group")
head(pdat)
library(ggplot2)
library(paletteer)
ggplot(pdat, aes(x = Y1, y = Y2, fill = group, color = group)) +
  geom_point(shape = 21, color = "black", size = 3) +  # 点使用黑边，填充颜色来自 group
  stat_ellipse(geom = "polygon", alpha = 0.3, linetype = 2) +  # 椭圆透明度
  scale_color_manual(values = my_palette) +  # 线条颜色
  scale_fill_manual(values = my_palette) +   # 填充颜色
  theme_classic() +
  theme(legend.position = "top")

#marker lncRNA热图
library(pheatmap)
library(ComplexHeatmap)
sorted_samples <- names(sort(group))
fpkm_sorted <- fpkm[, sorted_samples]      # 表达矩阵按样本排序
group_sorted <- group[sorted_samples]      # 分组向量对应排序后的样本
fpkm_num <- apply(fpkm_sorted, 2, as.numeric)
rownames(fpkm_num) <- rownames(fpkm_sorted)
annotation_col <- data.frame(Group = factor(group_sorted))
rownames(annotation_col) <- colnames(fpkm_num)
group_levels <- levels(factor(group_sorted))
ann_colors <- list(Group = setNames(my_palette[1:length(group_levels)], group_levels))
pheatmap(fpkm_num,
         cluster_rows = TRUE,        # 聚类基因
         cluster_cols = FALSE,       # 样本按分组排序，不聚类
         scale = "row",              # 按行标准化
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_colnames = FALSE,       # 显示样本名
         show_rownames = TRUE,       # 显示基因名
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_legend = TRUE)

#计算肿瘤纯度与免疫分数
KIRC_count <- read.csv("E:/RCD/data/cancer/KIRC/count/TCGA_KIRC_pcg_tum_count.csv",row.names = 1)

library(estimate)  #定义ESTIMATE包
estimate <- function(pcg_count,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(pcg_count,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   #注意platform
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='KIRC'#记得改名称
scores=estimate(KIRC_count,pro)
scores <- as.data.frame(scores)
scores$TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])
write.csv(scores,"KIRC_immunescores.csv")

#画箱线图
group <- read.csv("KIRC_group.csv",row.names = 1)
rownames(group) <- gsub("-", ".", rownames(group))
library(stringr)
sample_names <- str_sub(rownames(scores), 1, 12)
keep_index <- !duplicated(sample_names)
scores <- scores[keep_index, ]
rownames(scores) <- sample_names[keep_index]
scores1 <- scores[rownames(group),]
identical(rownames(group), rownames(scores1))
data <- cbind(group,scores1)
library(ggplot2)
library(ggpubr)
p = ggboxplot(data,x = "x",y = "StromalScore",color = "x",add = "jitter")+#TumorPurity,StromalScore,ImmuneScore,ESTIMATEScore
  scale_color_manual(values = c("#98d09d","#fbf398","#e77381","#9b8191","#7EABCA"))+
  stat_compare_means()#多组差异
p

#Cibersort计算免疫细胞
setwd("D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC亚型")
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
group <- read.csv("KIRC_group.csv",row.names = 1)
fpkm <- read.csv("E:/RCD/data/cancer/KIRC/fpkm/TCGA_KIRC_pcg_tum_fpkm.csv",row.names = 1)
rownames(group) <- gsub("-", ".", rownames(group))
library(stringr)
fpkm <- as.data.frame(t(fpkm))
sample_names <- str_sub(rownames(fpkm), 1, 12)
keep_index <- !duplicated(sample_names)
fpkm <- fpkm[keep_index, ]
rownames(fpkm) <- sample_names[keep_index]
fpkm = 2^fpkm-1#要求输入未log处理数据
fpkm <- fpkm[!duplicated(rownames(fpkm)),]
fpkm <- as.data.frame(t(fpkm))
fpkm <- cbind(RowName = rownames(fpkm), fpkm)
rownames(fpkm) <- NULL
colnames(fpkm)[1] <- "Gene Symbol"
head(fpkm)
write.csv(fpkm,"cibersort_data.csv",row.names = T)
write.table(fpkm, file = "cibersort_data.txt", sep = "\t", row.names = F, col.names = TRUE, quote = FALSE)

setwd("C:/Users/Dell/Desktop/data/subtype(2)")
source('Cibersort.R')
result <- CIBERSORT('LM22.txt','cibersort_data.txt', perm = 1000, QN = T)
result <- read.table("CIBERSORT-Results.txt",header = T,sep="\t",row.names = 1)

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(reshape2)
b <- read.csv("KIRC_group.csv",row.names = 1)
rownames(b) <- gsub("-", ".", rownames(b))
a <- result[,1:22]
identical(rownames(a),rownames(b))#查看行名是否一致
b <- b[rownames(a), ]
rownames(a) <- substr(rownames(a), 1, 12)
rownames(b) <- substr(rownames(b), 1, 12)
sum(rownames(a) %in% rownames(b))
common <- intersect(rownames(a), rownames(b))
a <- a[common, ]
b <- b[common, , drop = FALSE]
identical(rownames(a), rownames(b))

a$group <- b$x
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=CIBERSORT,value = Proportion,-c(group,sample))
b$group <- as.factor(b$group)
pdf("cibersort.pdf",width = 10,height = 6)
ggboxplot(b, 
          x = "CIBERSORT", 
          y = "Proportion",
          fill = "group") +
  scale_fill_manual(values = c("#98d09d","#fbf398","#e77381","#9b8191","#7EABCA")) +
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args = list(
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("****", "***", "**", "*", "ns"))
  ) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#免疫细胞比例与marker基因的相关性
genes <- c("RP5-1091N2.9","RP11-291B21.2","RP11-327F22.2","LINC01260","RP11-10J5.1","LINC00944","CTD-2023M8.1","USP30.AS1","AC079466.1",
           "RP4-764O22.1","LINC01428","RP11-807H17.1")
a <- read.table("CIBERSORT-Results.txt",sep = "\t",row.names = 1,check.names = F,header = T)
a <- a[,1:22]
expr <- read.csv("E:/RCD/data/cancer/KIRC/fpkm/TCGA_KIRC_lncRNA_tum_fpkm.csv",row.names =1)
genes_expr <- as.data.frame(t(expr[rownames(expr) %in% genes,]))
genes_expr <- genes_expr[rownames(a), ]
rownames(a) <- substr(rownames(a), 1, 12)
rownames(genes_expr) <- substr(rownames(genes_expr), 1, 12)
sum(rownames(a) %in% rownames(genes_expr))
common <- intersect(rownames(a), rownames(genes_expr))
a <- a[common, ]
genes_expr <- genes_expr[common, , drop = FALSE]
identical(rownames(a),rownames(genes_expr))

library(linkET)
library(tidyr)
library(tibble)
library(dplyr)
cor_res <- correlate(genes_expr, a ,method = "spearman")

df_r <- cor_res$r %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-1,names_to = "cell_type", values_to = "correlation")
df_p <- cor_res$p %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-1,names_to = "cell_type", values_to = "pvalue")
df_cor <- df_r %>%
  left_join(df_p) %>%
  # 增加FDR校正
  mutate(
    p_adj = p.adjust(pvalue, method = "BH"),   # Benjamini–Hochberg FDR
    stars = cut(
      p_adj,
      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
      right = FALSE,
      labels = c("***", "**", "*", " ")
    )
  )
df_filtered <- df_cor %>%
  filter(!cell_type %in% c("Correlation", "P-value", "RMSE"))
head(df_filtered)
write.csv(df_filtered,"Gene_cibersort_cor.csv")

df_matrix <- df_filtered %>%
  select(gene, cell_type, correlation) %>%
  pivot_wider(names_from = cell_type, values_from = correlation) %>%
  column_to_rownames("gene")
gene_dist <- dist(df_matrix)
gene_clust <- hclust(gene_dist)# 使用层次聚类
df_filtered$gene <- factor(df_filtered$gene, levels = gene_clust$labels[gene_clust$order])# 对 gene 进行聚类后的排序

library(ggplot2)
pdf("cibersort_cor_FDR.pdf",width = 15,height = 8)
ggplot(df_filtered, aes(cell_type, gene)) +
  geom_tile(aes(fill = correlation)) +
  geom_text(aes(label = stars), color = "black", size = 4) +
  scale_fill_gradient2(low = "lightblue", high = "lightcoral", mid = "white", limit = c(-1, 1),
                       name = paste0("*    FDR < 0.05", "\n\n",
                                     "**  FDR < 0.01", "\n\n",
                                     "*** FDR < 0.001", "\n\n",
                                     "Correlation"))+
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.ticks.y = element_blank(),
        panel.background = element_blank())
dev.off()

library(readxl)
data <- read_excel("mmc2.xlsx")
table(data$`TCGA Study`)
KIRC <- data[data$`TCGA Study` == "KIRC",]
write.csv(KIRC,"KIRC_mmc2.csv")

data <- read.csv("KIRC_mmc2.csv",row.names = 1)
group <- read.csv("KIRC_group.csv")
sample <- intersect(data$TCGA.Participant.Barcode,group$X)
data1 <- data[data$TCGA.Participant.Barcode %in% sample,]
group1 <- group[group$X %in% sample,]
library(dplyr)
norm_id <- function(x){
  x <- gsub("\\.", "-", x)
  substr(x, 1, 12)
}
group1$sample <- norm_id(group1$X)
data1$sample  <- norm_id(data1$TCGA.Participant.Barcode)
# 合并标签（只保留两方都有的样本）
df <- group1 %>%
  inner_join(data1, by = "sample") %>%
  dplyr::select(sample, my_subtype = x, immune_subtype = Immune.Subtype, tcga_subtype = TCGA.Subtype)
head(df)
library(vcd)
library(mclust)
tab_mmc2  <- table(df$my_subtype, df$immune_subtype)
tab_tcga  <- table(df$my_subtype, df$tcga_subtype)
# 卡方 + Cramer's V
chisq.test(tab_mmc2)
assocstats(tab_mmc2)$cramer

chisq.test(tab_tcga)
assocstats(tab_tcga)$cramer

adjustedRandIndex(df$my_subtype, df$immune_subtype)
adjustedRandIndex(df$my_subtype, df$tcga_subtype)

library(ggalluvial)
library(ggplot2)
library(grid)
if (!exists("gg_par")) gg_par <- grid::gpar
df$my_subtype <- factor(df$my_subtype, levels = sort(unique(df$my_subtype)))
my_colors <- c("#98d09d","#fbf398","#e77381","#9b8191")
### ---------------- Sankey 1: My → MMC2 ----------------
df_mmc2 <- df %>% 
  filter(!is.na(immune_subtype)) %>%  # 去除 NA
  count(my_subtype, immune_subtype)   # 统计频数

pdf("Sankey_my_vs_mmc2_COUNT.pdf", width = 6, height = 4)
ggplot(df_mmc2,
       aes(axis1 = my_subtype, axis2 = immune_subtype, y = n)) +
  geom_alluvium(aes(fill = my_subtype), width = 1/8, alpha = 0.9) +
  geom_stratum(width = 1/8, color = "grey30") +
  geom_label(stat = "stratum",
             aes(label = paste0(after_stat(stratum), " (n=", after_stat(count), ")")),
             size = 4, fill="white") +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(limits = c("My Subtype", "MMC2 Subtype")) +
  ylab("Sample Count") + xlab("") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=13, face="bold"),
    axis.text.x = element_text(color="black", size=11),
    axis.text.y = element_text(size=10),
    panel.grid = element_blank()
  )
dev.off()
### ---------------- Sankey 2: My → TCGA ----------------
df_tcga <- df %>%
  filter(!is.na(tcga_subtype)) %>%
  count(my_subtype, tcga_subtype)

pdf("Sankey_my_vs_TCGA_COUNT.pdf", width = 6, height = 4)
ggplot(df_tcga,
       aes(axis1 = my_subtype, axis2 = tcga_subtype, y = n)) +
  geom_alluvium(aes(fill = my_subtype), width = 1/8, alpha = 0.9) +
  geom_stratum(width = 1/8, color = "grey30") +
  geom_label(stat = "stratum",
             aes(label = paste0(after_stat(stratum), " (n=", after_stat(count), ")")),
             size = 4, fill="white") +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(limits = c("My Subtype", "TCGA Subtype")) +
  ylab("Sample Count") + xlab("") +
  ggtitle("Sankey: My subtype → TCGA subtype") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size=13, face="bold"),
    axis.text.x = element_text(color="black", size=11),
    axis.text.y = element_text(size=10),
    panel.grid = element_blank()
  )
dev.off()
######分组的差异表达
setwd("D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC亚型/DEGs_enrich")
count <- read.csv("E:/RCD/data/cancer/KIRC/count/TCGA_KIRC_pcg_count.csv",row.names = 1)
count <- as.data.frame(sapply(count, as.numeric))
rownames(count) <- rownames(read.csv("E:/RCD/data/cancer/KIRC/count/TCGA_KIRC_pcg_count.csv", row.names = 1))
count = 2^count-1
count <- round(count)
group <- read.csv("D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC亚型/Immune/KIRC_group.csv")
group1 <- group[which(group$x=="1"),]
group2 <- group[which(group$x=="2"),]
group3 <- group[which(group$x=="3"),]
group4 <- group[which(group$x=="4"),]
library(DESeq2)
match_samples <- function(sample_names, all_colnames) {
  matched <- sapply(sample_names, function(s) {
    s_clean <- gsub("-", ".", s)  # 把 - 替换成 .
    matched_cols <- grep(s_clean, all_colnames, value = TRUE)
    if (length(matched_cols) > 0) matched_cols[1] else NA
  })
  matched <- matched[!is.na(matched)]
  return(matched)
}
group1_samples <- match_samples(group1$X, colnames(count))
group2_samples <- match_samples(group2$X, colnames(count))
group3_samples <- match_samples(group3$X, colnames(count))
group4_samples <- match_samples(group4$X, colnames(count))
tumor_samples_all <- unique(c(group1_samples, group2_samples, group3_samples, group4_samples))
normal_samples <- setdiff(colnames(count), tumor_samples_all)

group_list <- list(group1 = group1_samples, group2 = group2_samples,
                   group3 = group3_samples, group4 = group4_samples)

for (g in names(group_list)) {
  tumor <- group_list[[g]]
  sub_count <- count[, c(tumor, normal_samples)]
  
  condition <- factor(c(rep("tumor", length(tumor)), rep("normal", length(normal_samples))))
  colData <- data.frame(condition = condition, row.names = colnames(sub_count))
  
  dds <- DESeqDataSetFromMatrix(countData = sub_count, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  
  sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 2), ]
  
  write.csv(as.data.frame(sig), file = paste0("DEG_", g, "_vs_normal.csv"))
  cat("✅ 完成:", g, "\n")
}
###分组的富集结果
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)
deg_path <- "D:/Project/03SCI/250710LncRCD修稿/GPB/补充数据/KIRC亚型/DEGs_enrich/"
deg_files <- list.files(deg_path, pattern = "DEG_.*_vs_normal\\.csv$", full.names = TRUE)

for (file in deg_files) {
  
  # 读取差异基因
  gene <- read.csv(file)
  gene <- gene[abs(gene$log2FoldChange) > 2 & gene$padj < 0.05, ]
  gene_symbols <- gene$X
  
  if (length(gene_symbols) == 0) {
    cat("⚠️ 无显著差异基因:", file, "\n")
    next
  }
  
  prefix <- sub("DEG_|_vs_normal\\.csv", "", basename(file))
  
  # ---------------- GO 富集 ----------------
  ego <- enrichGO(
    gene = gene_symbols,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  ego_sig <- ego@result %>% filter( pvalue < 0.05)
  
  if (nrow(ego_sig) > 0) {
    write.csv(ego_sig, paste0("GO_", prefix, "_enrichment.csv"), row.names = FALSE)
    
    pdf(paste0("GO_dotplot_", prefix, ".pdf"))
    print(dotplot(ego, showCategory = 20, title = paste0("GO Enrichment - ", prefix)))
    dev.off()
  }
  
  # ---------------- KEGG 富集 ----------------
  gene_KEGG <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db, drop = TRUE)
  
  if (nrow(gene_KEGG) == 0) {
    cat("⚠️ 无 KEGG ENTREZID:", prefix, "\n")
    next
  }
  
  ekegg <- enrichKEGG(
    gene = gene_KEGG$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    use_internal_data = TRUE
  )
  
  if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
    
    # 使用后缀避免重名，确保 mutate 成功
    kegg_category <- kegg_category %>% mutate(id = paste0("hsa", id))
    KEGG_sig <- ekegg@result %>%
      dplyr::left_join(
        kegg_category %>% dplyr::select(id, name, category, subcategory),
        by = c("ID" = "id"),
        suffix = c("", "_cat")  # 避免重名
      ) %>%
      dplyr::mutate(
        Description = ifelse(!is.na(name) & name != "", name, Description),
        category = ifelse(!is.na(category_cat) & category_cat != "", category_cat, NA_character_),
        subcategory = ifelse(!is.na(subcategory_cat) & subcategory_cat != "", subcategory_cat, NA_character_)
      ) %>%
      dplyr::select(-name, -category_cat, -subcategory_cat) #%>%
      #dplyr::filter(!is.na(p.adjust) & p.adjust < 0.05)
    
    write.csv(KEGG_sig, paste0("KEGG_", prefix, "_enrichment.csv"), row.names = FALSE)
    
    pdf(paste0("KEGG_dotplot_", prefix, ".pdf"))
    print(dotplot(ekegg, showCategory = 20, title = paste0("KEGG Enrichment - ", prefix)))
    dev.off()
  }
  
  cat("✅ 完成 GO + KEGG 富集:", prefix, "\n\n")
}

###功能韦恩图
library(VennDiagram)
group1 <- read.csv("GO_group1_vs_normal.csv_enrichment.csv")
group1 <- group1[,3]
group2 <- read.csv("GO_group2_vs_normal.csv_enrichment.csv")
group2 <- group2[,3]
group3 <- read.csv("GO_group3_vs_normal.csv_enrichment.csv")
group3 <- group3[,3]
group4 <- read.csv("GO_group4_vs_normal.csv_enrichment.csv")
group4 <- group4[,3]

veen1 <- venn.diagram(x=list(group1,group2,group3,group4),
                      scaled = F, # 根据比例显示大小
                      alpha= 0.5, #透明度
                      lwd=1,lty=1,col=c("#98d09d","#fbf398","#e77381","#9b8191"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
                      label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
                      cex = 2, # 数字大小
                      fontface = "bold",  # 字体粗细；加粗bold 
                      fill=c("#98d09d","#fbf398","#e77381","#9b8191"), # 填充色 配色https://www.58pic.com/
                      category.names = c("group1", "group2","group3","group4") , #标签名
                      cat.dist = c(0.2, 0.2, 0.1, 0.1), # 标签距离圆圈的远近
                      cat.pos = c(-20, 20, -20, 20), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
                      cat.cex = 2, #标签字体大小
                      cat.fontface = "bold",  # 标签字体加粗
                      cat.col=c("#98d09d","#fbf398","#e77381","#9b8191"),   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
                      cat.default.pos = "outer",  # 标签位置, outer内;text 外
                      filename=NULL
)
pdf('GO_venn.pdf', width = 10, height = 10)
grid.draw(veen1)
dev.off()

df_inter <- get.venn.partitions(list(group1,group2,group3,group4))
for (i in 1:nrow(df_inter)) df_inter[i,'values'] <- paste(df_inter[[i,'..values..']], collapse = ', ')
df_inter[-c(5, 6)]
df_inter <- df_inter[,-6]
write.csv(df_inter,"GO_venn.csv")

