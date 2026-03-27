library(ggm)
library(tidyverse)
library(dplyr)
library(utils)
library(stringr)
library(fgsea)
library(ImmuLncRNA)

pcg_fpkm <- readRDS("C:/Users/Dell/Desktop/ImmuLnc/pcg_fpkm.RData")
lnc_fpkm <- readRDS("C:/Users/Dell/Desktop/ImmuLnc/lnc_fpkm.RData")
pathway <- readRDS("C:/Users/Dell/Desktop/ImmuLnc/GBM_immune.RData")
TumorPurity <- readRDS("C:/Users/Dell/Desktop/ImmuLnc/TumorPurity.RData")
a1<-names(TumorPurity)
a2<-gsub("-",".",a1)
names(TumorPurity) <- a2
head(TumorPurity)

#计算相关系数
par.cor <- function(mRNA_exp,lncRNA_exp,tum_pur,adjusted){    #定义函数
  fun_mtx_pcr <- function(x,y,z){
    r12=cor(t(x),t(y))
    r13=cor(t(x),z)
    r23=cor(z,t(y))
    r123=r13%*%r23
    rup=r12-r123
    rd1=sqrt(1-r13*r13)
    rd2=sqrt(1-r23*r23)
    rd=rd1%*%rd2
    rrr=rup/rd
    return(rrr)
  }
  inter_samples <- intersect(intersect(colnames(mRNA_exp),colnames(lncRNA_exp)),names(tum_pur))
  inter_tumpur <- tum_pur[inter_samples]
  if(length(inter_tumpur)<1){stop("no same samples")}
  if (adjusted){
    mRNA_exp_inter <- mRNA_exp[,inter_samples]
    lncRNA_exp_inter <- lncRNA_exp[,inter_samples]
  }else{
    mRNA_exp_inter_ori <- mRNA_exp[,inter_samples]
    lncRNA_exp_inter_ori <- lncRNA_exp[,inter_samples]
    mRow30 <- which(apply(mRNA_exp_inter_ori,1,function(v){return((sum(v==0)/length(v))>=0.3)}))
    mRemv <- mRow30
    if(length(mRemv)==0){
      mRNA_out0 <- mRNA_exp_inter_ori
    }else{
      mRNA_out0 <- mRNA_exp_inter_ori[-(mRemv),]
    }
    mRNA_exp_inter <- log2(mRNA_out0+0.001)
    lncRow50 <- which(apply(lncRNA_exp_inter_ori,1,quantile,probs=0.5)==0)
    lncRow90 <- which(apply(lncRNA_exp_inter_ori,1,quantile,probs=0.9)<=0.1)
    lncRemv <- union(lncRow50,lncRow90)
    if(length(lncRemv)==0){
      lncRNA_out0 <- lncRNA_exp_inter_ori
    }else{
      lncRNA_out0 <- lncRNA_exp_inter_ori[-(lncRemv),]
    }
    lncRNA_exp_inter <- log2(lncRNA_out0+0.001)
  }
  n=length(inter_samples)
  gn=1
  pcor <- fun_mtx_pcr(lncRNA_exp_inter,mRNA_exp_inter,inter_tumpur)
  statistic <- pcor*sqrt((n-2-gn)/(1-pcor^2))
  p.value <- 2*pnorm(-abs(statistic))
  rownames(pcor) <- rownames(lncRNA_exp_inter) ; rownames(p.value) <- rownames(lncRNA_exp_inter)
  colnames(pcor) <- rownames(mRNA_exp_inter) ; colnames(p.value) <- rownames(mRNA_exp_inter)
  par.cor.res <- list(pcor,p.value)
  names(par.cor.res) <- c("pcor.value","p.value")
  return(par.cor.res)
}

test_res <- par.cor(pcg_fpkm,lnc_fpkm,TumorPurity1,adjusted=TRUE)
sum(is.infinite(test_res$pcor.value)|is.na(test_res$pcor.value)|is.nan(test_res$pcor.value))
sum(is.infinite(test_res$p.value)|is.na(test_res$p.value)|is.nan(test_res$p.value))
test_res$pcor.value[is.na(test_res$pcor.value)|is.infinite(test_res$pcor.value)|is.nan(test_res$pcor.value)] <- 0
test_res$p.value[is.na(test_res$p.value)|is.infinite(test_res$p.value)|is.nan(test_res$p.value)] <- 0

###lncRNA与显著pcg偏相关系数
rownames(pcor2)=pcor2[,1] #第一列为行名
pcor2=pcor2[,-1]   
row_name<-rownames(pcor2)
list_interaction<-list()
func1<-function(x){
  a1<-which(x>=0.65|x<=(-0.65))
  a2<-row_name[a1]
  return(a2)
}
list_interaction<-apply(pcor2,2,func1)
data.table::rbindlist(list_interaction$PPP3CB,use.names=TRUE, fill=TRUE, idcol="Gene")

#确定相关lnc
str(pathway)
k=0.999 
options(warn = -1)

library(fgsea,warn.conflicts=F)
pValue_ori <- test_res$p.value
pcorValue_ori <- test_res$pcor.value
RS <- -log10(pValue_ori)*sign(pcorValue_ori)
fgseaRes_all <- c()
for(i in 1:nrow(RS)){
  if(sum( is.infinite(RS[i,])) != 0){
    next()
  }
  ranks <- RS[i,]
  fgseaRes <- fgsea(pathway, ranks, minSize=1, maxSize=5000, nperm=1000)
  sigValue <- c()
  for(j in 1:nrow(fgseaRes)){
    if(fgseaRes$ES[j]>0){
      sig_ij <- 1 - 2*fgseaRes$pval[j]
    }else{
      sig_ij <- 2*fgseaRes$pval[j] - 1
    }
    sigValue <- c(sigValue,sig_ij)
    
  }
  lncRNA <- rownames(RS)[i]
  fgseaRes_i <- cbind(lncRNA,fgseaRes,sigValue)
  fgseaRes_all <- rbind(fgseaRes_all,fgseaRes_i)
}
sig_ind <- which(abs(fgseaRes_all$sigValue) >= k)
sig_pairs <- fgseaRes_all[sig_ind,1:2]
sig_pairs <- as.matrix(sig_pairs)
gsea.Res <- list(sig_pairs,fgseaRes_all)
names(gsea.Res) <- c("sig_pairs","fgseaRes_all")

test_res$sig_pairs[1:2,]
test_res$fgseaRes_all[1:2,]
write.csv(sig_pairs,file = "C:/Users/Dell/Desktop/ImmuLnc/sig_pairs.csv")
fgseaRes_all1 <- select(fgseaRes_all,-leadingEdge)
write.csv(fgseaRes_all1,file = "C:/Users/Dell/Desktop/ImmuLnc/fgseaRes_all.csv")
saveRDS(gsea.Res,file = "C:/Users/Dell/Desktop/ImmuLnc/gsea.Res.RData")

data <- read.csv("E:/RCD/data/cancer/UCEC/rcdlnc/fgseaRes_all.csv",row.names = 1)
data <- data[data$padj<0.05,]
head(data)
library(dplyr)
library(ggplot2)
cutoffs <- c(0.99, 0.991, 0.992,0.993, 0.994, 0.995, 0.996, 0.997, 0.998,0.999)
count_list <- numeric(length(cutoffs))
for (i in 1:length(cutoffs)) {
  threshold <- cutoffs[i]
  # 筛选绝对值大于当前阈值的行，并统计唯一的 lncRNA 数量
  filtered_data <- data %>% filter(abs(sigValue) > threshold)
  count_list[i] <- n_distinct(filtered_data$lncRNA) 
}

# 3. 将结果整合为数据框，用于画图
plot_df <- data.frame(
  Threshold = cutoffs,
  lncRNA_Count = count_list
)

# 查看一下统计结果
print(plot_df)

data <- read.csv("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/data",row.names =1)
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
  
  # 关键：添加一条加粗的红色虚线在 0.995 处，贯穿所有折线
  geom_vline(xintercept = 0.995, linetype = "dashed", color = "red", size = 1.2) +
  
  # 在 0.995 处添加文字标注
  annotate("text", x = 0.994, y = max(plot_df$lncRNA_Count), 
           label = "Threshold = 0.995", color = "red", angle = 90, vjust = -0.5, fontface = "bold") +
  
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
  scale_x_continuous(breaks = c(0.9, 0.95, 0.98, 0.99, 0.995, 0.999))

# 3. 显示图表
print(p_overlap)

library(dplyr)
library(stringr)
standardize_lncRNA_name <- function(names_vec) {
  names_vec <- toupper(str_trim(names_vec)) # 去除空格并转大写
  # 还可以添加更多清洗规则，例如替换连字符等，但目前这两个是最常见的。
  return(names_vec)
}
# EVLncRNAs
evl_lncs <- EVLncRNAs %>%
  select(any_of(c("LncRNA name", "LncRNA.name"))) %>% 
  distinct() %>%
  setNames("lnc") %>% # 临时统一列名方便操作
  mutate(lnc = standardize_lncRNA_name(lnc)) %>%
  filter(!is.na(lnc) & lnc != "") %>%
  pull(lnc)

# Lnc2cancer
lnc2c_lncs <- Lnc2cancer %>%
  select(name) %>% # 选择 lncRNA 名称列
  distinct() %>%
  mutate(name = standardize_lncRNA_name(name)) %>%
  filter(name != "") %>%
  pull(name)

# LncRNADisease
LncRNADisease <- LncRNADisease[LncRNADisease$`ncRNA Category` == "LncRNA",]
lncrnad_lncs <- LncRNADisease %>%
  select(any_of(c("ncRNA Symbol", "ncRNA.Symbol"))) %>% 
  distinct() %>%
  setNames("lnc") %>%
  mutate(lnc = standardize_lncRNA_name(lnc)) %>%
  filter(!is.na(lnc) & lnc != "") %>%
  pull(lnc)

# 你的 RCD_lnc 列表
your_rcd_lncs <- RCD_lnc %>%
  select(lncRNA) %>%
  distinct() %>%
  mutate(lncRNA = standardize_lncRNA_name(lncRNA)) %>%
  filter(lncRNA != "") %>%
  pull(lncRNA)

result_df <- tibble(
  lncRNA = your_rcd_lncs
) %>%
  mutate(In_EVLncRNAs = lncRNA %in% evl_lncs) %>%
  mutate(In_Lnc2cancer = lncRNA %in% lnc2c_lncs) %>%
  mutate(In_LncRNADisease = lncRNA %in% lncrnad_lncs)

result_df <- result_df %>%
  mutate(Confidence_Label = case_when(
    # 逻辑 1：只要在 Lnc2cancer 中出现过，就是高置信度
    In_Lnc2cancer == TRUE ~ "High Confidence",
    
    # 逻辑 2：不在 Lnc2cancer，但在另外两个数据库中至少出现了一个
    In_Lnc2cancer == FALSE & (In_EVLncRNAs == TRUE | In_LncRNADisease == TRUE) ~ "Medium Confidence",
    
    # 逻辑 3：三个数据库都没出现
    In_EVLncRNAs == FALSE & In_Lnc2cancer == FALSE & In_LncRNADisease == FALSE ~ "Candidate",
    
    # 备用：万一有遗漏的情况（理论上不会）
    TRUE ~ "Others"
  ))

# 将结果转换为因子并排序，方便后续绘图（High -> Medium -> Candidate）
result_df$Confidence_Label <- factor(result_df$Confidence_Label, 
                                     levels = c("High Confidence", "Medium Confidence", "Candidate"))
print(table(result_df$Confidence_Label))
head(result_df)
