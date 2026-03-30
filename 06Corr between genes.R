###提取RCD相关pcg的表达谱###
rm(list = ls())
limma_pcg <- read.csv("C:/Users/Dell/Desktop/RCD-lnc/figure/Fig/Fig1/geneset.csv")
pcg_fpkm <- read.csv("F:/RCD/cancer/LUSC/fpkm/TCGA_LUSC_pcg_tum_fpkm.csv",row.names = 1)
pcg_fpkma=pcg_fpkm[limma_pcg$pcg,]
pcg_fpkm0 <-as.matrix(pcg_fpkma)
pcg_fpkm0[is.infinite(pcg_fpkm0)|is.na(pcg_fpkm0)|is.nan(pcg_fpkm0)] <- 0
pcg_fpkm0=pcg_fpkm0[which(rowSums(pcg_fpkm0) > 0),]
dim(pcg_fpkm0)
pcg_fpkm <- as.data.frame(pcg_fpkm0)
saveRDS(pcg_fpkm,"F:/RCD/cancer/LUSC/cor/LUSC_geneset_fpkm.Rdata")
###提取RCD相关lncRNA表达谱###
lnc_fpkm <- readRDS("F:/RCD/cancer/LUSC/rcdlnc/LUSC_lnc_fpkm.Rdata")
lnc <- read.csv("F:/RCD/cancer/LUSC/cor/LUSC_rcdlnc.csv")
lnc_fpkma=lnc_fpkm[lnc$lncRNA,]
lnc_fpkm0 <-as.matrix(lnc_fpkma)
lnc_fpkm0[is.infinite(lnc_fpkm0)|is.na(lnc_fpkm0)|is.nan(lnc_fpkm0)] <- 0
lnc_fpkm0=lnc_fpkm0[which(rowSums(lnc_fpkm0) > 0),]
dim(lnc_fpkm0)
lnc_fpkm1 <- as.data.frame(lnc_fpkm0)
saveRDS(lnc_fpkm1,"F:/RCD/cancer/LUSC/cor/LUSC_rcd_lnc_fpkm.Rdata")
###求相关系数###
rm(list = ls())
pcg_fpkm <- readRDS("F:/RCD/cancer/LUSC/cor/LUSC_geneset_fpkm.Rdata")
lnc_fpkm <- readRDS("F:/RCD/cancer/LUSC/cor/LUSC_rcd_lnc_fpkm.Rdata")
a1 <- as.matrix(lnc_fpkm)
a1 <- t(a1)
a2 <- as.matrix(pcg_fpkm)
a2 <- t(a2)
corresult<-data.frame(gene=character(0),lnc=character(0),
                      pair=character(0),R=numeric(0),
                      P=numeric(0))
g=1
for(i in 1:ncol(a1)){
  for (j in 1:ncol(a2)) {
    c1<-cor(as.numeric(a1[,i]),as.numeric(a2[,j]),method = "pearson")
    c2<-cor.test(as.numeric(a1[,i]),as.numeric(a2[,j]),method = "pearson")$p.value
    corresult<-rbind(corresult,data.frame(gene=colnames(a1)[i],lnc=colnames(a2)[j],
                                          pair="lnc-gene",R=c1,P=c2))
    g=g+1
    print(g)
  }
  saveRDS(corresult,file = "lnc_pcg_cor.Rdata")
}
p <- as.data.frame(corresult$P,romnames = rownames(corresult))
FDR <- p.adjust(p$`corresult$P`,method = "BH")
corresult$FDR <- FDR
corresult1 <-corresult[which(abs(corresult$R)>0.7&corresult$FDR<0.05),]
write.csv(corresult1,file = "LUSC_lnc_corresult7.csv")

#阈值选取
data <- read.csv("E:/RCD/data/cancer/UCEC/cor/UCEC_lnc_corresult4.csv",row.names = 1)
library(dplyr)
library(ggplot2)
cutoffs <- c(0.4, 0.5, 0.6, 0.7)
results_list <- list()
for (i in 1:length(cutoffs)) {
  threshold <- cutoffs[i]
  # 筛选绝对值大于当前阈值的行
  filtered_data <- data %>% filter(abs(R) > threshold)
  # 统计各项指标
  results_list[[i]] <- data.frame(
    Threshold = threshold,
    Row_Count = nrow(filtered_data),              # 保留的行数（Pair数量）
    Unique_lncRNA = n_distinct(filtered_data$gene), # 唯一的lncRNA数量（对应图中gene列）
    Unique_Gene = n_distinct(filtered_data$lnc)    # 唯一的基因数量（对应图中lnc列）
  )
}
plot_df <- bind_rows(results_list)
print(plot_df)

data <- read.csv("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/lncRES/RCD_lnc_cor_num.csv",row.names =1)
library(tidyr)
library(tibble)
plot_df <- data %>%
  # 将行名转换为名为 "CancerType" 的新列 (如果你的癌症名已经是第一列，请注释掉这行)
  rownames_to_column(var = "CancerType") %>%
  # 将所有以 X 开头的列转换成长格式
  pivot_longer(
    cols = starts_with("X"), 
    names_to = "Threshold", 
    values_to = "lncRNA_Count"
  ) %>%
  # 去掉 Threshold 列名中的 "X" 并转换为数字类型
  mutate(Threshold = as.numeric(gsub("^X", "", Threshold)))
print(head(plot_df))

p_overlap <- ggplot(plot_df, aes(x = Threshold, y = lncRNA_Count, group = CancerType, color = CancerType)) +
  
  # 绘制所有癌症的折线
  geom_line(size = 0.8, alpha = 0.7) + 
  
  # 绘制所有点
  geom_point(size = 1.5) +

  # 设置经典主题
  theme_classic() +
  
  # 细节美化
  labs(
    title = "Sensitivity Analysis: lncRNA Count vs. lncRES Threshold",
    subtitle = "Consistent trend observed across 18 cancer types",
    x = "Absolute lncRES Threshold (|sigValue| > x)",
    y = "Number of Retained Unique lncRNAs",
    color = "Cancer Type" # 图例标题
  ) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11, color = "black"),
    # 如果癌症太多，图例可以放在右侧或底部
    legend.position = "right", 
    legend.text = element_text(size = 9)
  ) +
  
  # 强制显示所有的测试阈值作为刻度
  scale_x_continuous(breaks = c(0.4, 0.5, 0.6, 0.7))

# 3. 显示图表
print(p_overlap)
