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
