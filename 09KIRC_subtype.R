setwd("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/KIRC/Subtytpe")
a <- read.csv("E:/RCD/data/cancer/KIRC/Apoptosis.csv")
b <- read.csv("E:/RCD/data/cancer/KIRC/Cuproptosis.csv")
c <- read.csv("E:/RCD/data/cancer/KIRC/Necroptosis.csv")
d <- read.csv("E:/RCD/data/cancer/KIRC/Pyroptosis.csv")
e <- rbind(a,b,c,d)
lnc <- unique(e$lncRNA)
rm(a,b,c,d,e)

library(data.table)
library(dplyr)
library(stringr)
library(tibble)
survival <- fread("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/KIRC/Prognosis/TCGA-KIRC.survival.tsv.gz", data.table = FALSE)
fpkm <- read.csv("E:/RCD/data/cancer/KIRC/fpkm/TCGA_KIRC_lncRNA_tum_fpkm.csv", row.names = 1)
colnames(fpkm) <- gsub("\\.", "-", colnames(fpkm))
surv2 <- survival %>%
  transmute(
    sample = as.character(sample),
    patient_id = as.character(`_PATIENT`),
    time = as.numeric(`OS.time`),
    status = as.numeric(OS)
  ) %>%
  filter(!is.na(sample), !is.na(patient_id), !is.na(time), !is.na(status)) %>%
  filter(time > 0) %>%
  distinct(patient_id, .keep_all = TRUE)
exp <- fpkm[row.names(fpkm) %in% lnc,]

library(NMF)
library(survival)
library(survminer)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(cluster)
library(mclust)
library(factoextra)
exp <- as.data.frame(exp)
exp <- as.matrix(exp)
exp_sample <- substr(colnames(exp), 1, 15)
colnames(exp) <- exp_sample
tumor_idx <- substr(colnames(exp), 14, 15) == "01"
exp_tumor <- exp[, tumor_idx, drop = FALSE]
dim(exp_tumor)
head(colnames(exp_tumor))

surv2 <- as.data.frame(surv2)
surv2$sample15 <- substr(surv2$sample, 1, 15)
surv2$patient15 <- substr(surv2$patient_id, 1, 12)

exp_patient <- substr(colnames(exp_tumor), 1, 12)
colnames(exp_tumor) <- exp_patient
exp_tumor <- exp_tumor[, !duplicated(colnames(exp_tumor)), drop = FALSE]

surv_use <- surv2 %>%
  mutate(patient12 = substr(patient_id, 1, 12)) %>%
  distinct(patient12, .keep_all = TRUE)
common_patients <- intersect(colnames(exp_tumor), surv_use$patient12)
length(common_patients)
exp_use <- exp_tumor[, common_patients, drop = FALSE]
surv_use2 <- surv_use %>%
  filter(patient12 %in% common_patients) %>%
  arrange(match(patient12, colnames(exp_use)))
all(colnames(exp_use) == surv_use2$patient12)
#单因素COX
cox_results <- lapply(rownames(exp_use), function(gene){
  x <- as.numeric(exp_use[gene, ])
  dat <- data.frame(
    time = surv_use2$time,
    status = surv_use2$status,
    expr = x
  )
  fit <- coxph(Surv(time, status) ~ expr, data = dat)
  s <- summary(fit)
  data.frame(
    gene = gene,
    HR = s$coefficients[1, "exp(coef)"],
    coef = s$coefficients[1, "coef"],
    pvalue = s$coefficients[1, "Pr(>|z|)"]
  )
})

cox_df <- do.call(rbind, cox_results)
cox_df <- cox_df[order(cox_df$pvalue), ]
head(cox_df)
sig_genes <- cox_df$gene[cox_df$pvalue < 0.05]
length(sig_genes)
write.csv(sig_genes,"cox_sig_genes.csv")
#NMF参数确定
nmf_mat <- exp_use[sig_genes, , drop = FALSE]
dim(nmf_mat)
set.seed(1234)
rank_range <- 2:8
nmf_est <- nmf(nmf_mat, rank = rank_range, method = "brunet", nrun = 150, seed = 1234)
pdf("NMF_rank_selection.pdf", width = 8, height = 8)
plot(nmf_est)
dev.off()
#NMF
best_k <- 3
set.seed(1234)
nmf_fit <- nmf(nmf_mat, rank = best_k, method = "brunet", nrun = 200, seed = 123456)
nmf_fit
clusters <- predict(nmf_fit)
table(clusters)
subtype_df <- data.frame(
  patient_id = names(clusters),
  subtype = paste0("C", clusters),
  stringsAsFactors = FALSE)
head(subtype_df)

feature_idx <- extractFeatures(nmf_fit, method = "max")
feature_genes <- unique(rownames(nmf_mat)[unlist(feature_idx)])

subtype_surv <- surv_use2 %>%
  mutate(patient_id = patient12) %>%
  inner_join(subtype_df, by = "patient_id")
head(subtype_surv)
table(subtype_surv$subtype)
#consensus matrix热图
cons_mat <- consensusmap(nmf_fit, annCol = data.frame(Subtype = factor(clusters)),
                         main = "Consensus map", labCol = NA, labRow = NA)
#KM
fit_km <- survfit(Surv(time, status) ~ subtype, data = subtype_surv)
p_km <- ggsurvplot(
  fit_km,
  data = subtype_surv,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"),
  ggtheme = theme_bw()
)
print(p_km)
#silhouette 分析
dist_mat <- dist(t(nmf_mat))
sil <- silhouette(as.numeric(factor(subtype_df$subtype)), dist_mat)
summary(sil)
plot(sil, border = NA, main = "Silhouette plot for NMF subtypes")
mean_sil <- summary(sil)$avg.width
mean_sil
seeds <- c(1, 11, 21, 31, 41)
cluster_list <- list()
for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  fit_tmp <- nmf(nmf_mat, rank = best_k, method = "brunet", nrun = 100, seed = seeds[i])
  cl_tmp <- predict(fit_tmp)
  cluster_list[[i]] <- cl_tmp
}
names(cluster_list) <- paste0("seed_", seeds)
ari_mat <- matrix(NA, nrow = length(cluster_list), ncol = length(cluster_list))
rownames(ari_mat) <- names(cluster_list)
colnames(ari_mat) <- names(cluster_list)
for(i in 1:length(cluster_list)){
  for(j in 1:length(cluster_list)){
    common <- intersect(names(cluster_list[[i]]), names(cluster_list[[j]]))
    ari_mat[i,j] <- adjustedRandIndex(cluster_list[[i]][common], cluster_list[[j]][common])
  }
}
ari_mat
pheatmap(ari_mat, display_numbers = TRUE, main = "ARI across different seeds")
#subsampling
set.seed(1234)
n_iter <- 50
subsample_ratio <- 0.8
subsample_results <- list()
all_patients <- colnames(nmf_mat)
for(i in 1:n_iter){
  sub_patients <- sample(all_patients, size = floor(length(all_patients) * subsample_ratio))
  sub_mat <- nmf_mat[, sub_patients, drop = FALSE]
  
  fit_sub <- nmf(sub_mat, rank = best_k, method = "brunet", nrun = 50, seed = i + 100)
  cl_sub <- predict(fit_sub)
  
  subsample_results[[i]] <- cl_sub
}
orig_cluster <- clusters
ari_sub <- c()
for(i in 1:length(subsample_results)){
  common <- intersect(names(orig_cluster), names(subsample_results[[i]]))
  ari_sub[i] <- adjustedRandIndex(orig_cluster[common], subsample_results[[i]][common])
}
summary(ari_sub)
ggplot(data.frame(ARI = ari_sub), aes(x = ARI)) +
  geom_histogram(bins = 15, fill = "#4DBBD5", color = "black") +
  theme_bw() +
  labs(title = "Subsampling robustness of NMF clustering")
#层次聚类
hc <- hclust(dist(t(nmf_mat)), method = "ward.D2")
hc_cluster <- cutree(hc, k = best_k)
table(hc_cluster)
adjustedRandIndex(clusters, hc_cluster)
plot(hc, labels = FALSE, main = "Hierarchical clustering")
rect.hclust(hc, k = best_k, border = 2:5)
#k-means
set.seed(1234)
km <- kmeans(t(nmf_mat), centers = best_k, nstart = 50)
km_cluster <- km$cluster
table(km_cluster)
adjustedRandIndex(clusters, km_cluster)
compare_df <- data.frame(
  patient_id = names(clusters),
  NMF = as.character(clusters[names(clusters)]),
  HC = as.character(hc_cluster[names(clusters)]),
  Kmeans = as.character(km_cluster[names(clusters)])
)
write.csv(compare_df, "Subtype_comparison_methods.csv", row.names = FALSE)

final_subtype <- subtype_surv %>%
  select(patient_id, sample, time, status, subtype)
write.csv(final_subtype, "KIRC_NMF_subtypes.csv", row.names = FALSE)
write.csv(cox_df, "Univariate_cox_lncRNAs_for_NMF.csv", row.names = FALSE)
write.csv(data.frame(selected_genes = sig_genes), "Selected_lncRNAs_for_NMF.csv", row.names = FALSE)
#Sankey图
library(ggalluvial)
my_cols <- c("NMF_1" = "#E64B35","NMF_2" = "#4DBBD5","NMF_3" = "#00A087","NMF_4" = "#3C5488")
plot_df <- compare_df %>%
  mutate(
    NMF = paste0("NMF_", NMF),
    HC = paste0("HC_", HC),
    Kmeans = paste0("Kmeans_", Kmeans)
  )
freq_df <- dplyr::count(plot_df, NMF, HC, Kmeans)
p2 <- ggplot(freq_df,
             aes(axis1 = NMF, axis2 = HC, axis3 = Kmeans, y = n)) +
  geom_alluvium(aes(fill = NMF), width = 0.2, alpha = 0.85) +
  geom_stratum(width = 0.2, fill = "grey95", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("NMF", "HC", "Kmeans"), expand = c(.08, .08)) +
  scale_fill_manual(values = my_cols) +
  labs(
    title = "Sankey plot of clustering concordance",
    x = NULL,
    y = "Number of samples",
    fill = "NMF subtype"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )
print(p)
#marker lncRNA热图
library(pheatmap)
library(ComplexHeatmap)
sorted_samples <- names(sort(clusters))
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
# 特征基因表达热图
fpkm <- nmf_mat[feature_genes, , drop = FALSE]
group <- clusters
group <- group[colnames(fpkm)] 
sorted_samples <- names(sort(group))
fpkm_sorted <- fpkm[, sorted_samples, drop = FALSE]
group_sorted <- group[sorted_samples]
fpkm_num <- as.matrix(fpkm_sorted)
storage.mode(fpkm_num) <- "numeric"
annotation_col <- data.frame(Group = factor(paste0("C", group_sorted)))
rownames(annotation_col) <- colnames(fpkm_num)
group_levels <- levels(annotation_col$Group)
my_palette <- c("#E64B35", "#4DBBD5", "#00A087")
ann_colors <- list(Group = setNames(my_palette[seq_along(group_levels)], group_levels))
pheatmap(fpkm_num,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         show_colnames = FALSE,
         show_rownames = TRUE,
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         annotation_legend = TRUE)
#estimate分析
KIRC_count <- read.csv("E:/RCD/data/cancer/KIRC/count/TCGA_KIRC_pcg_tum_count.csv",row.names = 1)
colnames(KIRC_count) <- gsub("\\.", "-", colnames(KIRC_count))
KIRC_count_cowmanes <- substr(colnames(KIRC_count), 1, 12)
colnames(KIRC_count) <- KIRC_count_cowmanes
KIRC_count <- KIRC_count[,colnames(nmf_mat)]
library(estimate)
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
#箱线图
rownames(scores) <- gsub("\\.", "-", rownames(scores))
identical(names(group), rownames(scores))
data <- cbind(group,scores)
library(ggplot2)
library(ggpubr)
p = ggboxplot(data,x = "group",y = "ESTIMATEScore",color = "group",add = "jitter")+#TumorPurity,StromalScore,ImmuneScore,ESTIMATEScore
  scale_color_manual(values = c("#E64B35", "#4DBBD5", "#00A087"))+
  stat_compare_means()#多组差异
p
##Cibersort分析
library(reshape2)
pcg_fpkm <- read.csv("E:/RCD/data/cancer/KIRC/fpkm/TCGA_KIRC_pcg_tum_fpkm.csv",row.names = 1)
colnames(pcg_fpkm) <- gsub("\\.", "-", colnames(pcg_fpkm))
KIRC_pcg_fpkm_colmanes <- substr(colnames(pcg_fpkm), 1, 12)
colnames(pcg_fpkm) <- KIRC_pcg_fpkm_colmanes
pcg_fpkm <- pcg_fpkm[,colnames(nmf_mat)]
library(stringr)
pcg_fpkm = 2^pcg_fpkm-1#要求输入未log处理数据
pcg_fpkm <- cbind(RowName = rownames(pcg_fpkm), pcg_fpkm)
rownames(pcg_fpkm) <- NULL
colnames(pcg_fpkm)[1] <- "Gene Symbol"
head(pcg_fpkm)
write.table(pcg_fpkm, file = "cibersort_data.txt", sep = "\t", row.names = F, col.names = TRUE, quote = FALSE)

source('Cibersort.R')
result <- CIBERSORT('LM22.txt','cibersort_data.txt', perm = 1000, QN = T)
result <- read.table("CIBERSORT-Results.txt",header = T,sep="\t",row.names = 1)
library(tidyr)
library(ggsci)
ciber <- result[,1:22]
identical(rownames(ciber),names(group))#查看行名是否一致
ciber$group <- group
ciber <- ciber %>% rownames_to_column("sample")
ciber$group <- factor(paste0("C", group[match(ciber$sample, names(group))]))
ciber_data <- ciber %>%
  pivot_longer(
    cols = -c(sample, group),
    names_to = "CIBERSORT",
    values_to = "Proportion"
  )
group_levels <- levels(ciber_data$group)
ggboxplot(ciber_data, 
          x = "CIBERSORT", 
          y = "Proportion",
          fill = "group") +
  scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087")) +
  stat_compare_means(aes(group = group),
                     method = "kruskal.test",
                     label = "p.signif",
                     symnum.args = list(
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                       symbols = c("****", "***", "**", "*", "ns"))
  ) +
  theme(text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))
##免疫细胞比例与marker基因的相关性
genes <- c("RP5-1172A22.1", "RP11-807H17.1", "RP4-764O22.1", "LINC01428", "USP30-AS1", "AC079466.1")
result <- read.table("CIBERSORT-Results.txt",header = T,sep="\t",row.names = 1)
ciber <- result[,1:22]
genes_expr <- as.data.frame(t(fpkm))
identical(rownames(ciber),rownames(genes_expr))
library(linkET)
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
#亚型比较
library(readxl)
data <- read_excel("mmc2.xlsx")
table(data$`TCGA Study`)
KIRC <- data[data$`TCGA Study` == "KIRC",]
write.csv(KIRC,"KIRC_mmc2.csv")

sample <- intersect(KIRC$`TCGA Participant Barcode`, names(group))
KIRC1 <- KIRC[KIRC$`TCGA Participant Barcode` %in% sample,]
group1 <- group[sample]
group1_df <- data.frame(
  sample = names(group1),
  x = paste0("C", as.character(group1)),
  stringsAsFactors = FALSE)
KIRC1$sample  <- KIRC1$`TCGA Participant Barcode`
df <- group1_df %>%
  inner_join(KIRC1, by = "sample") %>%
  dplyr::select(sample, my_subtype = x, immune_subtype = `Immune Subtype`)
head(df)
library(vcd)
library(mclust)
tab_mmc2  <- table(df$my_subtype, df$immune_subtype)
# 卡方 + Cramer's V
chisq.test(tab_mmc2)
assocstats(tab_mmc2)$cramer
adjustedRandIndex(df$my_subtype, df$immune_subtype)
