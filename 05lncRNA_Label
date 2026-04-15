#####1.Disease-Label
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
#####2.LncTard——Label
library(data.table)
lnctard <- read.table("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/lnctard/lnctard2.0.txt", 
                      header = TRUE, 
                      sep = "	", 
                      fill = TRUE,      # 核心参数：允许行长度不等，空位补 NA
                      quote = "",       # 核心参数：不对引号进行转义处理
                      comment.char = "", # 防止把 # 号当成注释
                      check.names = FALSE)
table(lnctard$RegulatorType)
table(lnctard$TargetType)
lnctard <- lnctard[lnctard$RegulatorType=="lncRNA",]

RCD_lnc <- read.csv("RCD-lnc.csv")
RCD_lnc_Tar <- lnctard[lnctard$Regulator %in% RCD_lnc$lncRNA,]
RCD_lnc_Tar_mRNA <- RCD_lnc_Tar[RCD_lnc_Tar$TargetType=="PCG",]

RCD_mRNA <- read.csv("E:/RCD/data/cancer/####RCD-pathway/geneset.csv",header = T)
library(dplyr)
rcd_lnc_interaction <- RCD_lnc_Tar_mRNA %>%
  # 只选择我们需要的列，减少冗余
  select(Regulator, Target) %>%
  # 连接两个表
  inner_join(RCD_mRNA, by = c("Target" = "pcg")) %>%
  # 重新排列一下列顺序，让结构更清晰
  select(Regulator, Target, pathway) %>%
  # 去重（防止源数据中有重复记录）
  distinct()
head(rcd_lnc_interaction)
table(rcd_lnc_interaction$Target)
write.csv(rcd_lnc_interaction, "LncRNA_RCD_Interaction_List.csv", row.names = FALSE)

#####3.Cis—Label
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
RCD_lnc <- read.csv("RCD-lnc.csv")
RCD_mRNA <- read.csv("E:/RCD/data/cancer/####RCD-pathway/geneset.csv",header = T)
probemap_data <- read.delim("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/CIS/gencode.v36.annotation.gtf.gene.probemap", 
                            header = TRUE, 
                            sep = "	", 
                            check.names = FALSE)
head(probemap_data)
colnames(probemap_data) <- trimws(colnames(probemap_data))
genes_anno <- probemap_data %>%
  dplyr::rename(
    seqnames = chrom, 
    gene_name = gene,
    start = chromStart,
    end = chromEnd
  ) %>%
  # 确保坐标是数值型
  mutate(
    start = as.numeric(start), 
    end = as.numeric(end)
  ) %>%
  filter(!is.na(start))
head(genes_anno)

lnc_coords <- genes_anno %>% 
  filter(gene_name %in% RCD_lnc$lncRNA) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

pcg_coords <- genes_anno %>% 
  filter(gene_name %in% RCD_mRNA$pcg) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

nearest_pairs_indices <- nearest(lnc_coords, pcg_coords, ignore.strand = TRUE)
nearest_pcg_coords <- pcg_coords[nearest_pairs_indices]
distances <- distance(lnc_coords, nearest_pcg_coords, ignore.strand = TRUE)
lnc_to_pcg_nearest_dist <- data.frame(
  lncRNA_gene_name = mcols(lnc_coords)$gene_name,
  lncRNA_id = mcols(lnc_coords)$id,
  pcg_gene_name = mcols(nearest_pcg_coords)$gene_name, # 使用 nearest_pcg_coords 的 meta-columns
  pcg_id = mcols(nearest_pcg_coords)$id,               # 使用 nearest_pcg_coords 的 meta-columns
  distance = distances
)
lnc_to_pcg_nearest_dist <- lnc_to_pcg_nearest_dist %>% filter(!is.na(pcg_gene_name))
head(lnc_to_pcg_nearest_dist)
final_res <- lnc_to_pcg_nearest_dist %>%
  mutate(type = case_when(
    is.na(distance) ~ "Different_Chr_or_Unknown", # 虽然前面 filter 掉了 NA，这里可以保留以防万一
    distance == 0 ~ "Overlap/Cis",
    distance <= 100000 ~ "Cis_Proximal(<100kb)",
    distance <= 500000 ~ "Cis_Distal(<500kb)",
    TRUE ~ "Trans"
  ))
head(final_res)
write.csv(final_res, "LncRNA_RCD_CIS_List.csv", row.names = FALSE)
library(ggplot2)

# 绘制距离分布图（只针对同一染色体上的基因，因为不同染色体距离是 NA）
ggplot(final_res %>% filter(!is.na(distance)), aes(x = distance/1000)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Distances between lncRNA and nearest PCG",
       x = "Distance (kb)",
       y = "Count") +
  geom_vline(xintercept = 100, linetype="dashed", color = "red") # 100kb 阈值线
  labs(title = "Zoomed-in View (0-500kb)", x = "Distance (kb)", y = "Count")

#整合标签
setwd("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/label")
lnctard <- read.csv("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/label/lnctard/LncRNA_RCD_Interaction_List.csv")
disease <- read.csv("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/label/Disease/lncRNA_in_database.csv",row.names = 1)
CIS <- read.csv("D:/Project/03SCI/250710_LncRCD/Computational and Structural Biotechnology Journal/Supplementary/label/CIS/LncRNA_RCD_CIS_List.csv")
lnctard <- unique(lnctard[,-3])

library(dplyr)
library(tidyr)
library(stringr)
collapse_cell <- function(x, sep = "; ") {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) NA_character_ else paste(x, collapse = sep)
}
cis2 <- CIS %>%
  rename(lncRNA = lncRNA_gene_name)
lnctard2 <- lnctard %>%
  rename(lncRNA = Regulator)
disease2 <- disease
if ("Confidence_Label" %in% colnames(disease2)) {
  disease2 <- disease2 %>% select(-Confidence_Label)
}
to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1)
  x <- tolower(trimws(as.character(x)))
  x %in% c("true", "t", "1", "yes", "y")
}
disease_presence <- disease2 %>%
  mutate(
    In_EVLncRNAs     = to_logical(In_EVLncRNAs),
    In_Lnc2cancer    = to_logical(In_Lnc2cancer),
    In_LncRNADisease = to_logical(In_LncRNADisease),
    in_disease = (In_EVLncRNAs %in% TRUE) | (In_Lnc2cancer %in% TRUE) | (In_LncRNADisease %in% TRUE)
  ) %>%
  group_by(lncRNA) %>%
  summarise(
    in_disease = any(in_disease, na.rm = TRUE),
    In_EVLncRNAs     = any(In_EVLncRNAs %in% TRUE, na.rm = TRUE),
    In_Lnc2cancer    = any(In_Lnc2cancer %in% TRUE, na.rm = TRUE),
    In_LncRNADisease = any(In_LncRNADisease %in% TRUE, na.rm = TRUE),
    .groups = "drop"
  )
cis_sum <- cis2 %>%
  group_by(lncRNA) %>%
  summarise(
    lncRNA_id     = collapse_cell(lncRNA_id),
    pcg_gene_name = collapse_cell(pcg_gene_name),
    pcg_id        = collapse_cell(pcg_id),
    distance      = collapse_cell(distance),
    type          = collapse_cell(type),
    in_cis = TRUE,
    .groups = "drop"
  )
lnctard_sum <- lnctard2 %>%
  group_by(lncRNA) %>%
  summarise(
    Target = collapse_cell(Target),
    in_lnctard = TRUE,
    .groups = "drop"
  )
all_lnc <- bind_rows(
  cis2 %>% distinct(lncRNA),
  disease2 %>% distinct(lncRNA),
  lnctard2 %>% distinct(lncRNA)
) %>%
  distinct(lncRNA)
final_table <- all_lnc %>%
  left_join(cis_sum, by = "lncRNA") %>%
  left_join(disease_presence, by = "lncRNA") %>%
  left_join(lnctard_sum, by = "lncRNA") %>%
  mutate(
    in_cis     = replace_na(in_cis, FALSE),
    in_lnctard = replace_na(in_lnctard, FALSE),
    in_disease = replace_na(in_disease, FALSE),
    score = as.integer(in_cis) + as.integer(in_disease) + as.integer(in_lnctard),
    Confidence_Label = case_when(
      score == 3 ~ "High Confidence",
      score == 2 ~ "Medium Confidence",
      score == 1 ~ "Low Confidence",
      score == 0 ~ "Candidate",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    lncRNA, score, Confidence_Label,
    in_cis, in_disease, in_lnctard,
    everything()
  )
head(final_table)
final_table_conf <- final_table %>%
  filter(Confidence_Label %in% c("High Confidence", "Medium Confidence", "Low Confidence"))
write.csv(final_table_conf, "merged_lncRNA_annotation.csv", row.names = FALSE)
