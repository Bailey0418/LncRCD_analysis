###原始数据处理###
#count数据处理
library(readr)
library(tidyverse)
library(limma)
CHOL_count <- read_tsv("F:/RCD/cancer/UCEC/count/TCGA-UCEC.htseq_counts.tsv")
probeMap <- read_tsv("F:/RCD/课题/GENE_DATA/gencode.gene.info.v22.tsv")
TCGA_gset <- CHOL_count %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "gene_id")) %>%
  select(gene_name, starts_with("TCGA") )
TCGA_gset[1:4,1:4]
TCGA_gset = as.data.frame(avereps(TCGA_gset[,-1],ID = TCGA_gset$gene_name) )
colnames(TCGA_gset) <- substring(colnames(TCGA_gset),1,15) %>% gsub("-",".",.)
TCGA_group_list <- ifelse(as.numeric(substring(colnames(TCGA_gset),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list)
pcg_info <- probeMap[which(probeMap$gene_type=="protein_coding"),]
lnc_info <- probeMap[which(probeMap$gene_type=="lincRNA"),]
#拆分保存正常与癌症count矩阵
mRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% pcg_info$gene_name,]
write.csv(mRNA_gset,"F:/RCD/cancer/UCEC/count/TCGA_UCEC_pcg_count.csv",quote = F,row.names = T)
mRNA_nor = mRNA_gset[,str_sub(colnames(mRNA_gset),14,15) > 10]
mRNA_tum = mRNA_gset[,str_sub(colnames(mRNA_gset),14,15)<=10]
write.csv(mRNA_nor,"F:/RCD/cancer/UCEC/count/TCGA_UCEC_pcg_nor_count.csv",quote = F,row.names = T)
write.csv(mRNA_tum,"F:/RCD/cancer/UCEC/count/TCGA_UCEC_pcg_tum_count.csv",quote = F,row.names = T)

lncRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% lnc_info$gene_name,]
write.csv(lncRNA_gset,"F:/RCD/cancer/BLCA/count/TCGA_BLCA_lnc_count.csv",quote = F,row.names = T)
lnc_nor = lncRNA_gset[,str_sub(colnames(lncRNA_gset),14,15) > 10]
lnc_tum = lncRNA_gset[,str_sub(colnames(lncRNA_gset),14,15) <= 10]
write.csv(lnc_nor,"F:/RCD/cancer/BLCA/count/TCGA_BLCA_lnc_nor_count.csv",quote = F,row.names = T)
write.csv(lnc_tum,"F:/RCD/cancer/BLCA/count/TCGA_BLCA_lnc_tum_count.csv",quote = F,row.names = T)

#fpkm数据处理
rm(list = ls())
CHOL_fpkm <- read_tsv("F:/RCD/cancer/UCEC/fpkm/TCGA-UCEC.htseq_fpkm.tsv")
probeMap <- read_tsv("F:/RCD/课题/GENE_DATA/gencode.gene.info.v22.tsv")
TCGA_gset <- CHOL_fpkm %>%
  inner_join(probeMap, by = c("Ensembl_ID" = "gene_id")) %>%
  select(gene_name, starts_with("TCGA") )
TCGA_gset = as.data.frame(avereps(TCGA_gset[,-1],ID = TCGA_gset$gene_name) )
colnames(TCGA_gset) <- substring(colnames(TCGA_gset),1,15) %>% gsub("-",".",.)

TCGA_group_list <- ifelse(as.numeric(substring(colnames(TCGA_gset),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list)
pcg_info <- probeMap[which(probeMap$gene_type=="protein_coding"),]
lnc_info <- probeMap[which(probeMap$gene_type=="lincRNA"),]
#拆分保存癌症与正常fpkm矩阵
mRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% pcg_info$gene_name,]
write.csv(mRNA_gset,"F:/RCD/cancer/UCEC/fpkm/TCGA_UCEC_pcg_fpkm.csv",quote = F,row.names = T)
mRNA_nor = mRNA_gset[,str_sub(colnames(mRNA_gset),14,15) > 10]
mRNA_tum = mRNA_gset[,str_sub(colnames(mRNA_gset),14,15)<=10]
write.csv(mRNA_nor,"F:/RCD/cancer/UCEC/fpkm/TCGA_UCEC_pcg_nor_fpkm.csv",quote = F,row.names = T)
write.csv(mRNA_tum,"F:/RCD/cancer/UCEC/fpkm/TCGA_UCEC_pcg_tum_fpkm.csv",quote = F,row.names = T)

lncRNA_gset <- TCGA_gset[rownames(TCGA_gset) %in% lnc_info$gene_name,]
write.csv(lncRNA_gset,"F:/RCD/cancer/BLCA/fpkm/TCGA_BLCA_lncRNA_fpkm.csv",quote = F,row.names = T)
lnc_nor = lncRNA_gset[,str_sub(colnames(lncRNA_gset),14,15) > 10]
lnc_tum = lncRNA_gset[,str_sub(colnames(lncRNA_gset),14,15)<=10]
write.csv(lnc_nor,"F:/RCD/cancer/BLCA/fpkm/TCGA_BLCA_lncRNA_nor_fpkm.csv",quote = F,row.names = T)
write.csv(lnc_tum,"F:/RCD/cancer/BLCA/fpkm/TCGA_BLCA_lncRNA_tum_fpkm.csv",quote = F,row.names = T)
#生成分析需要的pcg_fpkm矩阵
rm(list = ls())
pcg_fpkm <- read.csv("F:/RCD/cancer/UCEC/fpkm/TCGA_UCEC_pcg_tum_fpkm.csv",row.names = 1)
pcg_fpkm0 <- as.matrix(pcg_fpkm)
pcg_fpkm0[is.infinite(pcg_fpkm0)|is.na(pcg_fpkm0)|is.nan(pcg_fpkm0)] <- 0
pcg_fpkm0=pcg_fpkm0[which(rowSums(pcg_fpkm0) > 0),]
pcg_fpkm <- as.data.frame(pcg_fpkm0)
saveRDS(pcg_fpkm,file = "F:/RCD/cancer/UCEC/rcdlnc/UCEC_pcg_fpkm.RData")
