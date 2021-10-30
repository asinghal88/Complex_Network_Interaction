##CN project
library(data.table)
d=fread('hyper.csv')
head(d)
library(anytime)
d$timestamp=d$Column1.4
colnames(d)=c('node_1','node_2','connection','date','time')
d$date=anytime(d$date,tz='GMT')
d$time=anytime(d$time,tz='GMT')
d$date<- format(as.POSIXct(d$date,format='%Y-%m-%d %H:%M:%S'),format='%Y-%m-%d')
d$time<- format(as.POSIXct(d$time,format='%Y-%m-%d %H:%M:%S'),format='%H:%M:%S')

library(plyr)
m=d[,c(1,2)]
m
l=count(m)
l

## 2196 unique pairs(interactions)
d1 <- data.frame(a=unlist(m, use.names = FALSE))
k=unique(d1)
#113 unique visitors in whole data
d1=d[d$date<='2009-07-01',c(1,2,4,5)]
m4=d[d$date<='2009-07-01',c(1,2)]
#test data
m0=d[d$date=='2009-07-01'& d$time>'14:00:00',c(1,2)]
d0<- data.frame(a=unlist(m0, use.names = FALSE))
k0=unique(d0)
#76 unique nodes in test data 

m01=d[d$date<='2009-07-01'& d$time<'14:00:00',c(1,2)]
d00<- data.frame(a=unlist(m01, use.names = FALSE))
k00=unique(d00)

#day1
m4$day1_6=ifelse(d1$date=='2009-06-29' & d1$time>'04:00:00'& d1$time<'06:00:00',1,0)
m4$day1_8=ifelse(d1$date =='2009-06-29'& d1$time >'06:00:00'& d1$time <'08:00:00',1,0)
m4$day1_10=ifelse(d1$date=='2009-06-29'& d1$time>'08:00:00'&d1$time<'10:00:00',1,0)
m4$day1_12=ifelse(d1$date=='2009-06-29'& d1$time>'10:00:00'&d1$time<'12:00:00',1,0)
m4$day1_14=ifelse(d1$date=='2009-06-29'& d1$time>'12:00:00'&d1$time<'14:00:00',1,0)
m4$day1_16=ifelse(d1$date=='2009-06-29'& d1$time>'14:00:00'&d1$time<'16:00:00',1,0)
m4$day1_18=ifelse(d1$date=='2009-06-29'& d1$time>'16:00:00'&d1$time<'18:00:00',1,0)
m4$day1_20=ifelse(d1$date=='2009-06-29'& d1$time>'18:00:00'&d1$time<'20:00:00',1,0)
m4$day1_22=ifelse(d1$date=='2009-06-29'& d1$time>'20:00:00'&d1$time<'22:00:00',1,0)
m4$day1_24=ifelse(d1$date=='2009-06-29'& d1$time>'22:00:00'&d1$time<'24:00:00',1,0)
#day2
m4$day2_6=ifelse(d1$date=='2009-06-30'& d1$time>'04:00:00'&d1$time<'06:00:00',1,0)
m4$day2_8=ifelse(d1$date=='2009-06-30'& d1$time>'06:00:00'&d1$time<'08:00:00',1,0)
m4$day2_10=ifelse(d1$date=='2009-06-30'& d1$time>'08:00:00'&d1$time<'10:00:00',1,0)
m4$day2_12=ifelse(d1$date=='2009-06-30'& d1$time>'10:00:00'&d1$time<'12:00:00',1,0)
m4$day2_14=ifelse(d1$date=='2009-06-30'& d1$time>'12:00:00'&d1$time<'14:00:00',1,0)
m4$day2_16=ifelse(d1$date=='2009-06-30'& d1$time>'14:00:00'&d1$time<'16:00:00',1,0)
m4$day2_18=ifelse(d1$date=='2009-06-30'& d1$time>'16:00:00'&d1$time<'18:00:00',1,0)
m4$day2_20=ifelse(d1$date=='2009-06-30'& d1$time>'18:00:00'&d1$time<'20:00:00',1,0)
m4$day2_22=ifelse(d1$date=='2009-06-30'& d1$time>'20:00:00'&d1$time<'22:00:00',1,0)
m4$day2_24=ifelse(d1$date=='2009-06-30'& d1$time>'22:00:00'&d1$time<'24:00:00',1,0)
#day3
m4$day3_6=ifelse(d1$date=='2009-07-01'& d1$time>'04:00:00'&d1$time<'06:00:00',1,0)
m4$day3_8=ifelse(d1$date=='2009-07-01'& d1$time>'06:00:00'&d1$time<'08:00:00',1,0)
m4$day3_10=ifelse(d1$date=='2009-07-01'& d1$time>'08:00:00'&d1$time<'10:00:00',1,0)
m4$day3_12=ifelse(d1$date=='2009-07-01'& d1$time>'10:00:00'&d1$time<'12:00:00',1,0)
m4$day3_14=ifelse(d1$date=='2009-07-01'& d1$time>'12:00:00'&d1$time<'14:00:00',1,0)
# m4$day3_16=ifelse(d1$date=='2009-07-01'& d1$time>'14:00:00'&d1$time<'16:00:00',1,0)
# m4$day3_18=ifelse(d1$date=='2009-07-01'& d1$time>'16:00:00'&d1$time<'18:00:00',1,0)
# m4$day3_20=ifelse(d1$date=='2009-07-01'& d1$time>'18:00:00'&d1$time<'20:00:00',1,0)
# m4$day3_22=ifelse(d1$date=='2009-07-01'& d1$time>'20:00:00'&d1$time<'22:00:00',1,0)
# m4$day3_24=ifelse(d1$date=='2009-07-01'& d1$time>'22:00:00'&d1$time<'24:00:00',1,0)




df4=aggregate(.~ node_1+node_2, m4, sum)
df5=ifelse(df4[,c(3:27)]>0,1,0)
df5=cbind(df4[,c(1:2)],df5)

df5$count=rowSums(df5[,3:27])
df6=df5[,c(1,2,28)]

write.csv(x = df6,file = 'out!.csv')
write.csv(x = df6,file = 'hyper_train.csv')




d3=d[d$date<='2009-07-01',c(1,2,4,5)]
m6=d[d$date<='2009-07-01',c(1,2)]
##test
m6$day3_16=ifelse(d3$date=='2009-07-01'& d3$time>'14:00:00'&d3$time<'16:00:00',1,0)
m6$day3_18=ifelse(d3$date=='2009-07-01'& d3$time>'16:00:00'&d3$time<'18:00:00',1,0)
m6$day3_20=ifelse(d3$date=='2009-07-01'& d3$time>'18:00:00'&d3$time<'20:00:00',1,0)
m6$day3_22=ifelse(d3$date=='2009-07-01'& d3$time>'20:00:00'&d3$time<'22:00:00',1,0)
m6$day3_24=ifelse(d3$date=='2009-07-01'& d3$time>'22:00:00'&d3$time<'24:00:00',1,0)

#Actual 
df7=aggregate(.~ node_1+node_2, m6, sum)
df8=ifelse(df7[,c(3:7)]>0,1,0)
df8=cbind(df7[,c(1:2)],df8)
df8$count=rowSums(df8[,3:7])
df9=df8[,c(1,2,8)]
#just giving 1 if interacted
connection=ifelse(df9[,3]>0,1,0)
df10=cbind(df9[,c(1:2)],connection)
write.csv(x = df10,'hyper_test.csv')
#predicted
df12=read.csv('out2!.csv')
#when threshold is .12

k5=table(df10$connection,df12$prob>=.12)
k5
TP=k5[2,2]
TN=k5[1,1]
FP=k5[1,2]
FN=k5[2,1]
prec=TP/(TP+FP)
recall=TP/(TP+FN)#TPR/Sensitivity 
FPR=FP / (TN + FP)    
accuracy=(TP+TN)/nrow(df12)
f_score=2*prec*recall/(prec+recall)

#ROC CuURVE
library(ROCR)
final=cbind(df12$prob,df10$connection)
pred <- prediction(df12$prob,df10$connection)
roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc.perf, main = "ROC Curve(Hypertext)", colorize = TRUE)
abline(a=0, b= 1)     
auc.perf = performance(pred, measure = "auc")
auc.perf@y.values

opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(roc.perf, pred))

#unique nodes which have interacted in test data and have also interacted in train data i.e 201
df11=df6[!(df6$count==0),c(1:2)]
df13=df10[!(df10$connection==0),c(1:2)]
df14=merge(df11,df13)
#so this means that out of positive examples people who have first time interacted are (321-201=120)



df15=read.csv('out!.csv')
df15=df15[,c(3,4)]
k1=count(df15)
k1$freq=NULL

plot(k1$count,k1$prob,type='l',col='blue',lwd=2,xlab='Count',ylab='Probability')

#graph
library('igraph')
graph=df12
graph$node_1=as.character(graph$node_1)
graph$node_2=as.character(graph$node_2)

gr=graph.data.frame(graph,directed = F)
plot(gr,layout=layout.fruchterman.reingold,edge.width=E(gr)$prob)
adj=get.adjacency(gr,attr='prob',sparse=F)
adj
transitivity(gr,type=c('globalundirected'))
degree(gr,v = df12$node_1,mode = c("total"))
eigen=eigen_centrality(gr,directed=F)
###



library(data.table)
require(matrixcalc)
edge_list <- fread("out!.csv", sep = "auto")
edge_list=edge_list[,c(1,2,4)]
G <- graph.data.frame(edge_list,directed=FALSE);
A1 <- as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE,attr="prob");
A2 <-  matrix.power(A1,2)
A3 <-  matrix.power(A1,3)
A <- A1+A2+A3

#taking single interaction
for(row in 1:nrow(A)) {
  for(col in 1:ncol(A)) {
    if (row>=col) {
      A[row,col]=0
    } 
  }
}
#wherever prob is greater than 1
for(row in 1:nrow(A)) {
  for(col in 1:ncol(A)) {
    if (A[row,col]>1) {
      A[row,col]=1
    } 
  }
}

#adjacency to edge list
g <- graph.adjacency(A,weighted=TRUE)
df <- get.data.frame(g)
head(df)
colnames(df) <- c('node_1','node_2','prob')

#Using the squared matrix for comparison
library(dplyr)
d0 <- merge(df, df10, by=c("node_1","node_2"))                           ## d11 is the edge list from squared matrix and d9 is the test dataset (find the intersection of the two using inner join)
d01 <- anti_join(df, d0, by = c("node_1","node_2"))                        ##remove the intersection from the edge list
d01$prob = 0
names(d01)[names(d01) == "prob"] <- "connection"
d09 <- merge(x = df10, y = d01, by = c("node_1","node_2","connection"), all = TRUE)         ##add this to the test dataset(using outer join)


b0 <- merge(df10, df, by=c("node_1","node_2"))                                     ##d11 is the edge list from squared matrix and d9 is the test dataset (find the intersection of the two using inner join)
b01 <- anti_join(df10, b0, by = c("node_1","node_2"))                                                   ##remove the intersection from the test dataset
names(b01)[names(b01) == "connection"] <- "prob"
d011 <- merge(x = df, y = b01, by = c("node_1","node_2","prob"), all = TRUE)                 ##add this to the edge list(using outer join)

#The above code is to make the # of rows same in both the edge list and the test dataset to compare.

library(ROCR)
final1=cbind(d011$prob,d09$connection)
pred1 <- prediction(d011$prob,d09$connection)
roc.perf1 = performance(pred1, measure = "tpr", x.measure = "fpr")
plot(roc.perf1, main = "ROC Curve(Hypertext)", colorize = TRUE)
abline(a=0, b= 1)     
auc.perf1 = performance(pred1, measure = "auc")
auc.perf1@y.values

opt.cut = function(perf1, pred1){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf1@x.values, perf1@y.values, pred1@cutoffs)
}
print(opt.cut(roc.perf1, pred1))
 
k5=table(d09$connection,d011$prob>=1)
k5
TP=k5[2,2]
TN=k5[1,1]
FP=k5[1,2]
FN=k5[2,1]
prec=TP/(TP+FP)
recall=TP/(TP+FN)#TPR
FPR=FP / (TN + FP)    
accuracy=(TP+TN)/nrow(d09)
f_score=2*prec*recall/(prec+recall)
