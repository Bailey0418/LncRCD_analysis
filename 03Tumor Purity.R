#计算肿瘤纯度
rm(list = ls())
tum <- read.csv("F:/RCD/cancer/UCEC/count/TCGA_UCEC_pcg_tum_count.csv",row.names = 1)
tum = tum[which(rowSums(tum) > 0),]
dim(tum)

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
pro='UCEC'#记得改名称
scores=estimate(tum,pro)
TumorPurity = cos(0.6049872018+0.0001467884 * scores[,3])#计算肿瘤纯度
saveRDS(TumorPurity,file = "F:/RCD/cancer/UCEC/rcdlnc/UCEC_TumorPurity.RData")
