###寻找网络中的HUB节点
library(igraph)
g<-graph_from_data_frame(THCA_web,directed = F)

#1、点度中心度
g1<-degree(g)
#mode=in点入度；out=点出度；total点度中心度，三者统称绝对点中心度
g1<-data.frame(g1)
#2、接近中心度
g2<-closeness(g,mode = "all")
g2<-data.frame(g2)
#3、中间中心度
g3<-betweenness(g)
g3<-data.frame(g3)
#4、点的特征向量中心度
g4<-evcent(g,scale = F)$vector
g4<-data.frame(g4)
#page.rank特征向量中心度
g5<-page.rank(g)$vector
g5<-data.frame(g5)

tuopu<-cbind(g1,g2,g3,g4,g5)
q1<-data.frame(gene1=rownames(g1),g1=g1$g1)
q1<-q1[order(q1$g1),]
q2<-data.frame(gene2=rownames(g2),g2=g2$g2)
q2<-q2[order(q2$g2),]
q3<-data.frame(gene3=rownames(g3),g3=g3$g3)
q3<-q3[order(q3$g3),]
q4<-data.frame(gene4=rownames(g4),g4=g4$g4)
q4<-q4[order(q4$g4),]
q5<-data.frame(gene5=rownames(g5),g5=g5$g5)
q5<-q5[order(q5$g5),]

qstatistic<-data.frame(gene=q1$gene1,V=NA,Q=NA)
for(i in 1:nrow(q1)){
  r5<-(which(q5$gene5==q1$gene1[i])-1)/nrow(q1)
  r4<-(which(q4$gene4==q1$gene1[i])-1)/nrow(q1)
  r3<-(which(q3$gene3==q1$gene1[i])-1)/nrow(q1)
  r2<-(which(q2$gene2==q1$gene1[i])-1)/nrow(q1)
  r1<-(i-1)/nrow(q1)
  v0=1
  v1=(-1)^0*v0*r5
  v2=(-1)^0*v1*r4+(-1)^1*v0/2*r4^2
  v3=(-1)^0*v2*r3+(-1)^1*v1/2*r3^2+(-1)^2*v0/6*r3^3
  v4=(-1)^0*v3*r2+(-1)^1*v2/2*r2^2+(-1)^2*v1/6*r2^3+
    (-1)^3*v0/24*r2^4
  v5=(-1)^0*v4*r1+(-1)^1*v3/2*r1^2+(-1)^2*v2/6*r1^3+
    (-1)^3*v1/24*r1^4+(-1)^4*v0/120*r1^5
  Q=factorial(5)*v5
  qstatistic[i,2]<-v5
  qstatistic[i,3]<-Q
}
