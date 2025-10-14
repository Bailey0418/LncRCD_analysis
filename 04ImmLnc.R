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
