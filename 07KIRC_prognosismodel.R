#####prognosis model#####
#cox_data
setwd("")
cancer_survival <- read_tsv("KIRC_sur.txt")
cancer_survival <- cancer_survival[,c(2:4)]
names(cancer_survival) <- c("PATIENT","OS","OS.time")
cancer_survival <- cancer_survival[!duplicated(cancer_survival), ]

cancer_fpkm <- read.csv("TCGA_KIRC_lncRNA_tum_fpkm.csv",row.names = 1)
a <- read.csv("E:/RCD/data/cancer/KIRC/Apoptosis.csv")
b <- read.csv("E:/RCD/data/cancer/KIRC/Cuproptosis.csv")
c <- read.csv("E:/RCD/data/cancer/KIRC/Necroptosis.csv")
d <- read.csv("E:/RCD/data/cancer/KIRC/Pyroptosis.csv")
e <- rbind(a,b,c,d)
lnc <- unique(e$lncRNA)
rm(a,b,c,d,e)
rcd_lnc_fpkm <- cancer_fpkm[lnc,]
rcd_lnc_fpkm <- t(rcd_lnc_fpkm)
rcd_lnc_fpkm <- as.data.frame(rcd_lnc_fpkm)
library(dplyr)
rcd_lnc_fpkm$PatientID <- sapply(rownames(rcd_lnc_fpkm), function(x) {
  parts <- unlist(strsplit(x, "\\."))
  paste(parts[1:3], collapse = "-")
})

# 按病人ID聚合（取均值）
rcd_lnc_fpkm_patient <- rcd_lnc_fpkm %>%
  group_by(PatientID) %>%
  summarise(across(everything(), mean))

# 设置行名为PatientID
rcd_lnc_fpkm_patient_mat <- as.data.frame(rcd_lnc_fpkm_patient)
rownames(rcd_lnc_fpkm_patient_mat) <- rcd_lnc_fpkm_patient_mat$PatientID
rcd_lnc_fpkm_patient_mat$PatientID <- NULL

common_samples <- intersect(rownames(rcd_lnc_fpkm_patient_mat), cancer_survival$PATIENT)
length(common_samples)
expr_mat <- rcd_lnc_fpkm_patient_mat[common_samples,]
survival_info <- cancer_survival %>%
  dplyr::filter(PATIENT %in% common_samples) %>%
  dplyr::arrange(match(PATIENT, common_samples))
all.equal(rownames(expr_mat), survival_info$PATIENT)

cox_data <- cbind(survival_info,expr_mat,stringsAsFactors = TRUE)
cox_data$OS <- as.integer(cox_data$OS)
cox_data$OS.time <- as.integer(cox_data$OS.time)
write.csv(cox_data,"KIRC_cox_data.csv")

#survival analysis
library(survival)
library(survminer)
pfilter <- 0.05
cox_data <- cox_data[,-1]
uniresult <- data.frame()#单因素cox回归结果
for(i in colnames(cox_data[,3:ncol(cox_data)])){   
  unicox <- coxph(Surv(time = OS.time, event = OS) ~ cox_data[,i], data = cox_data)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  if(pvalue<pfilter){ 
    uniresult <- rbind(uniresult,
                       cbind(gene=i,
                             HR=unisum$coefficients[,2],
                             L95CI=unisum$conf.int[,3],
                             H95CI=unisum$conf.int[,4],
                             pvalue=unisum$coefficients[,5]
                       ))
  }
}   
write.csv(uniresult,file = "KIRC_单因素COX分析结果.csv",row.names = F)

unigene <- subset(cox_data,select = c("OS","OS.time",uniresult$gene))
multicox <- coxph(Surv(time = OS.time,event = OS) ~ ., data = unigene) 
multisum <- summary(multicox)
gene <- colnames(unigene)[3:ncol(unigene)]
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
multiresult <- data.frame(gene=gene,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue)
multiresult <- multiresult[multiresult$pvalue<0.05,]
write.csv(multiresult,file = "KIRC_多因素COX分析结果.csv",row.names = F)

#risk score
cancer_multi_lnc <- cox_data[, c(1:2, which(colnames(cox_data) %in% multiresult$gene))]
for(i in 1:nrow(cancer_multi_lnc)){#risk score计算和分组
  cancer_multi_lnc[i,12] <- (-1.252852)*cancer_multi_lnc[i,3]+0.9136963*cancer_multi_lnc[i,4]+ 0.4826291*cancer_multi_lnc[i,5]+
    (-0.9343253)*cancer_multi_lnc[i,6]+ 1.752132*cancer_multi_lnc[i,7]+ 0.9766936*cancer_multi_lnc[i,8]+
    0.3591962*cancer_multi_lnc[i,9]+ 0.3718698*cancer_multi_lnc[i,10]+ 0.4897705*cancer_multi_lnc[i,11]
}
names(cancer_multi_lnc)[12] <- "risk_score"
res.cut <- surv_cutpoint(cancer_multi_lnc, #数据集
                         time = "OS.time", #生存时间
                         event = "OS", #生存状态
                         variables = "risk_score" #需要计算的数据列名
)
summary(res.cut)
plot(res.cut, "risk_score", palette = "npg")
res.cat <- surv_categorize(res.cut)
head(res.cat)
write.csv(res.cat,"riskscore_group.csv")
cancer_multi_lnc$group <- res.cat$risk_score
fit <- survfit(Surv(OS.time, OS) ~ group, data = cancer_multi_lnc)
print(fit)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # 添加风险表
           risk.table.col = "strata", # 根据分层更改风险表颜色
           linetype = "strata", # 根据分层更改线型
           surv.median.line = "hv", # 同时显示垂直和水平参考线
           ggtheme = theme_bw(), # 更改ggplot2的主题
           palette = c("#E7B800", "#2E9FDF"))#定义颜色

#model ROC
library(survivalROC)
cancer_multi_lnc$OS.time=cancer_multi_lnc$OS.time/365
par(mar= c(5,5,1,1),cex.lab=1.2,cex.axis= 1.2) #先设置一下图形的边界
sROC=survivalROC(Stime=cancer_multi_lnc$OS.time, status=cancer_multi_lnc$OS, marker = cancer_multi_lnc$risk_score, predict.time =1, method="KM")#先画一个1年的图
plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red", 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=paste0("1 years"," (AUC=",sprintf("%.3f",sROC$AUC),")") #这个后面添加legend用
sROC3=survivalROC(Stime=cancer_multi_lnc$OS.time, status=cancer_multi_lnc$OS, marker = cancer_multi_lnc$risk_score, predict.time =3, method="KM")#再加一个3年的线
lines(sROC3$FP, sROC3$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="green",lwd = 2)
aucText3=paste0("3 years"," (AUC=",sprintf("%.3f",sROC3$AUC),")") #这个后面添加legend用
sROC5=survivalROC(Stime=cancer_multi_lnc$OS.time, status=cancer_multi_lnc$OS, marker = cancer_multi_lnc$risk_score, predict.time =5, method="KM")#再加一个5年的线
lines(sROC5$FP, sROC5$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="blue",lwd = 2)
aucText5=paste0("5 years"," (AUC=",sprintf("%.3f",sROC5$AUC),")") #这个后面添加legend用
legend("bottomright", c(aucText,aucText3,aucText5),lwd=2,bty="n",col=c("red","green","blue"),cex=1.2)#添加legend

library(timeROC)
tROC <-timeROC(T=cancer_multi_lnc$OS.time,delta = cancer_multi_lnc$OS,marker = cancer_multi_lnc$risk_score,
               cause = 1,times = c(1,3,5),ROC=T)
par(mar= c(5,5,1,1),cex.lab=1.5,cex.axis= 1.2) #设置图形边界
plot(tROC,time=1,col="red",title=F,lwd=2) #1年ROC
plot(tROC,time=3,col="green",add=T,title=F,lwd=2) #3年ROC
plot(tROC,time=5,col="blue",add=T,title=F,lwd=2) #5年ROC
legend(0,1, # 这里设置legend的位置坐标：横坐标0,纵坐标1。也可以像前面一样设置为"bottomright"等
       c(paste0("AUC at 1 years  ", round(tROC$AUC[1], 2)),
         paste0("AUC at 3 years  ", round(tROC$AUC[2], 2)),
         paste0("AUC at 5 years  ", round(tROC$AUC[3], 2))),
       col=c("red","green","blue"),lwd=2,cex=1.2,bty="n")
