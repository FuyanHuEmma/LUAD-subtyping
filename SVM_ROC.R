                                                                                                                                                                                          .setwd('/media/hufuyan/Data/networkproject_COAD/cancer_result/cancer_network/SVM_ROC/')
#setwd('/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/ROC_test_validation/')
setwd('/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/fd_pvalue/ROC/')

library(MASS)
library(ROCR)
library(kernlab)


Path_train<-'/media/hufuyan/Data/networkproject2017/cancer/RESM/'
Path_test<-'/media/hufuyan/Data/networkproject2017/cancer/script_subtypes/validation_data/GSE40419_LC-87_RPKM_expression/'

Path_our30<-'/media/hufuyan/Data/networkproject2017/cancer/script_subtypes/'

###training data
specialgene30<-read.table(paste(Path_our30,"specialgene30.txt",sep=''),check.names=FALSE,)
specialgene30<-as.matrix(specialgene30)
dim(specialgene30)





LUAD_RSEM_best_tumor<-read.table(paste(Path_train,"LUAD_RSEM_best_tumor.txt",sep=''),check.names=FALSE,row.names=1,header=1,fill=TRUE,quote = "")
LUAD_RSEM_best_tumor<-as.matrix(LUAD_RSEM_best_tumor)
dim(LUAD_RSEM_best_tumor)


LUAD_RSEM_best_normal<-read.table(paste(Path_train,"LUAD_RSEM_best_normal.txt",sep=''),check.names=FALSE,row.names=1,header=1,fill=TRUE,quote = "")
LUAD_RSEM_best_normal<-as.matrix(LUAD_RSEM_best_normal)
dim(LUAD_RSEM_best_normal)

##get train normal data for 30
train_normal<-matrix(0,30,ncol(LUAD_RSEM_best_normal))
for (i in 1:30){
  t<-which(rownames(LUAD_RSEM_best_normal)==specialgene30[i])
  train_normal[i,]<-LUAD_RSEM_best_normal[t,]
}

train_normal<-as.matrix(train_normal)
rownames(train_normal)<-specialgene30
colnames(train_normal)<-colnames(LUAD_RSEM_best_normal)
dim(train_normal)

##get train tumor data for 30
train_tumor<-matrix(0,30,ncol(LUAD_RSEM_best_tumor))
for (i in 1:30){
  t<-which(rownames(LUAD_RSEM_best_tumor)==specialgene30[i])
  train_tumor[i,]<-LUAD_RSEM_best_tumor[t,]
}

train_tumor<-as.matrix(train_tumor)
rownames(train_tumor)<-specialgene30
colnames(train_tumor)<-colnames(LUAD_RSEM_best_tumor)
dim(train_tumor)


###train data
normal_tumor_data<-cbind(train_normal,train_tumor)
xtrain<-t(normal_tumor_data)
dim(xtrain)
ytrain<-c(rep(0,ncol(train_normal)),rep(1,ncol(train_tumor)))
#ytrain
length(ytrain)


data_TCGA<-cbind(xtrain,ytrain)
colnames(data_TCGA)<-c(colnames(xtrain),'Type')

rownames(data_TCGA)<-rownames(xtrain)
data_TCGA<-as.data.frame(data_TCGA)
dim(data_TCGA)

###testing data


GSE40419_RPKM_exp_tumor<-read.csv(paste(Path_test,"GSE40419_RPKM_exp_tumor.csv",sep=''),check.names=FALSE,row.names=1,header=1)
GSE40419_RPKM_exp_tumor<-as.matrix(GSE40419_RPKM_exp_tumor)
dim(GSE40419_RPKM_exp_tumor)


GSE40419_RPKM_exp_normal<-read.csv(paste(Path_test,"GSE40419_RPKM_exp_normal.csv",sep=''),check.names=FALSE,row.names=1,header=1)
GSE40419_RPKM_exp_normal<-as.matrix(GSE40419_RPKM_exp_normal)
dim(GSE40419_RPKM_exp_normal)

test_normal<-matrix(0,30,ncol(GSE40419_RPKM_exp_normal))
for (i in 1:30){
  t<-which(rownames(GSE40419_RPKM_exp_normal)==specialgene30[i])
  #print (specialgene30[i])
  test_normal[i,]<-GSE40419_RPKM_exp_normal[t,]
}

test_normal<-as.matrix(test_normal)
rownames(test_normal)<-specialgene30
colnames(test_normal)<-colnames(GSE40419_RPKM_exp_normal)
dim(test_normal)

##get test tumor data for 30
test_tumor<-matrix(0,30,ncol(GSE40419_RPKM_exp_tumor))
for (i in 1:30){
  t<-which(rownames(GSE40419_RPKM_exp_tumor)==specialgene30[i])
  test_tumor[i,]<-GSE40419_RPKM_exp_tumor[t,]
}

test_tumor<-as.matrix(test_tumor)
rownames(test_tumor)<-specialgene30
colnames(test_tumor)<-colnames(GSE40419_RPKM_exp_tumor)
dim(test_tumor)



xtest<-t(cbind(test_normal,test_tumor))
dim(xtest)
ytest<-c(rep(0,ncol(test_normal)),rep(1,ncol(test_tumor)))

length(ytest)






#######
dataDivide<-function(col,data,n,group_num){
  #col is a facotr type column,divide each group of the dataframe 
  #into n partitions,string type
  #data is a data.frame type in R
  #n represents the numbers which you want to divide into,default 5
  #the function return a list contain n data.frame
  #use sample(x) generate x numbers in unordered rank,then
  #divide the x numebr into n partitions
  
  #group_num=length(levels(factor(data[,col])))#
  lst1=list() #按照因子分類把原數據分成group_num份
  lst2=list() #分成等分的數據框
  lst3=list() #
  t<-which(data[,col]==group_num)
  lst1=data[t,]
 
  
  od=sample(nrow(lst1)) ##平均分成n份，並且是隨機分的，需要用到sample函數
  newdata=lst1[od,]
  len=length(od)
  cutpoint=floor(len/n)
  for(j in 1:n) {
    if(len>=cutpoint*(1+j))
    {
      lst2[[j]]=newdata[(cutpoint*(j-1)+1):(cutpoint*j),]
    }
    else
    {
      lst2[[j]]=newdata[(cutpoint*(j-1)+1):len,]
    }
  }
  return(lst2)
}
 


col<-ncol(data_TCGA)
data<-data_TCGA
n=5
group_num<-0
head(data)

##divid 58 normal into 5 parts
Train_normal_divid<-dataDivide(ncol(data_TCGA),data_TCGA,5,0)
str(Train_normal_divid)
length(Train_normal_divid)


##divid 518 tumor into 10 parts
Train_tumor_divid<-dataDivide(ncol(data_TCGA),data_TCGA,10,1)
str(Train_tumor_divid[[1]])
length(Train_tumor_divid)
dim(Train_tumor_divid[[1]])


####divid each part of 10 parts into 5 parts
Train_tumor_divid2<-list()
for (i in 1:10){
  Train_tumor_divid2[[i]]<-dataDivide(ncol(Train_tumor_divid[[i]]),Train_tumor_divid[[i]],5,1)
}
str(Train_tumor_divid2[[1]])
dim(Train_tumor_divid2[[1]][[1]])

getwd()

##svm classifier

ypredscore<-c()
true_label<-c()
for (i in 1:10){
  print (i)
  
  for (j in 1:5){
    train_data_normal<-c()
    train_data_tumor<-c()
    y_train_normal<-c()
    y_train_tumor<-c()
    for (p in 1:5) {
      if (p!=j){
        train_data_normal<-rbind(train_data_normal,Train_normal_divid[[p]][,1:30])  ## 4/5 for training
        dim(train_data_normal)
        y_train_normal<-c(y_train_normal,Train_normal_divid[[p]][,31])
        train_data_tumor<-rbind(train_data_tumor,Train_tumor_divid2[[i]][[p]][,1:30])
        dim(train_data_tumor) 
        y_train_tumor<-c(y_train_tumor,Train_tumor_divid2[[i]][[p]][,31])
      }
    }
    dim(train_data_normal)
    dim(train_data_tumor)
    x_train_normal_tumor<-rbind(train_data_normal,train_data_tumor)
    y_train_normal_tumor<-c(y_train_normal,y_train_tumor)
    print (dim(x_train_normal_tumor))
    length(y_train_normal_tumor) 
    
    sample_selected<-c(rownames(train_data_normal),rownames(train_data_tumor))
    t_label<-c()
    for (k in 1:length(sample_selected)){
      t_label[k]<-which(sample_selected[k]==rownames(xtrain))
    }
    x_test_normal_tumor<-xtrain[-t_label,]
    y_test_normal_tumor<-ytrain[-t_label]
    print (dim(x_test_normal_tumor))
    length(y_test_normal_tumor)
    svp <- ksvm(as.matrix(x_train_normal_tumor),factor(y_train_normal_tumor),kernel="vanilladot",type="C-svc",prob.model=TRUE)
    svp
    dim(xtrain)
    
    ypredscore = c(ypredscore,predict(svp,x_test_normal_tumor,type="probabilities")[,2])
    true_label<-c(true_label,y_test_normal_tumor)
  }
}
length(ypredscore) 
length(true_label)
range(ypredscore)
pred <- prediction(ypredscore,true_label==1)
length(pred)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
str(perf)

getwd()
pdf('ROC_curve_TCGA_data.pdf')
plot(perf,col='red')
#str(perf@y.values)
#as.vector(unlist(perf@y.values))
#plot(as.vector(unlist(perf@x.values)),as.vector(unlist(perf@y.values)))
#write.csv(perf,row.names=TRUE,file='our30_TCGA_perf.csv')
abline(a=0,b=1)
auc <- performance(pred,measure="auc")@y.values
auc
legend('right',legend = paste('AUC:',round(auc[[1]],6),sep=''),col = c("red"),lty= 1,lwd = 1)
dev.off()
range(ypredscore)
head(ypredscore)
pred
getwd()

####for external validation data
dim(train_normal)
dim(train_tumor)
t<-sample((ncol(train_tumor)),70)  ##random selectin 40 tumor samples
dim(train_tumor[,t])
normal_tumor_data<-cbind(train_normal,train_tumor[,t])
dim(normal_tumor_data)
xtrain2<-t(normal_tumor_data)
dim(xtrain2)
ytrain2<-c(rep(0,ncol(train_normal)),rep(1,ncol(train_tumor[,t])))
#ytrain[41]
length(ytrain2)

dim(xtest)
svp2 <- ksvm(xtrain2,factor(ytrain2),kernel="vanilladot",type="C-svc",prob.model=TRUE)
ypredscore2 <-predict(svp2,xtest,type="probabilities")
dim(xtest)
dim()
pred2 <- prediction(ypredscore2[,1],ytest==0)
range(ypredscore2[,1])
head(ypredscore2)
perf2 <- performance(pred2, measure = "tpr", x.measure = "fpr")
#str(perf2)
str(pred2)
dev.off()
getwd()
pdf('ROC_curve_Validation_data.pdf')
plot(perf2,col='red')
abline(a=0,b=1)
auc2 <- performance(pred2,measure="auc")@y.values
auc2
legend('right',legend = paste('AUC:',round(auc2[[1]],6),sep=''),col = c("red"),lty= 1,lwd = 1)
dev.off()

perf3 <- performance(pred2, measure = "prec", x.measure = "rec")
str(perf3)
plot(perf3)
plot(perf3,colorize=T,colorkey.pos='top',print.cutoffs.at=0.9,text.cex=1,text.adj=c(1.2,1.2,lwd=2))

acc.perf<- performance(pred2, measure ="acc")
plot(acc.perf)



######################
######non-probability model
dim(xtest)
svp2 <- ksvm(xtrain2,factor(ytrain2),kernel="vanilladot",type="C-svc")
ypredscore2 <-predict(svp2,xtest,type='decision')
dim(xtest)
dim()
pred2 <- prediction(ypredscore2,factor(ytest))

perf2 <- performance(pred2, measure = "tpr", x.measure = "fpr")
#str(perf2)
str(pred2)
dev.off()
getwd()
pdf('ROC_curve_Validation_data.pdf')
plot(perf,col='red')
abline(a=0,b=1)
auc2 <- performance(pred2,measure="auc")@y.values
auc2
legend('right',legend = paste('AUC:',round(auc2[[1]],5),sep=''),col = c("red"),lty= 1,lwd = 1)
dev.off()

perf3 <- performance(pred2, measure = "prec", x.measure = "rec")
str(perf3)
plot(perf3)
plot(perf3,colorize=T,colorkey.pos='top',print.cutoffs.at=0.9,text.cex=1,text.adj=c(1.2,1.2,lwd=2))

acc.perf<- performance(pred2, measure ="acc")
plot(acc.perf)





##using e1071 packages
library('e1071')
dim(xtest)
svp2 <- svm(xtrain2,factor(ytrain2),kernel="linear",type="C-classification",probability=TRUE)
ypredscore2 <-predict(svp2,xtest,probability=TRUE)
dim(xtest)


pred2 <- prediction(attr(ypredscore2,'probabilities')[,2],ytest==1)
range(attr(ypredscore2,'probabilities')[,2])
length(attr(ypredscore2,'probabilities')[,2])
perf2 <- performance(pred2, measure = "tpr", x.measure = "fpr")
#str(perf2)
str(pred2)
dev.off()
getwd()
pdf('ROC_curve_Validation_data.pdf')
plot(perf2,col='red')
abline(a=0,b=1)
auc2 <- performance(pred2,measure="auc")@y.values
auc2
legend('right',legend = paste('AUC:',round(auc2[[1]],5),sep=''),col = c("red"),lty= 1,lwd = 1)
dev.off()

perf3 <- performance(pred2, measure = "prec", x.measure = "rec")
str(perf3)
plot(perf3)
plot(perf3,colorize=T,colorkey.pos='top',print.cutoffs.at=0.9,text.cex=1,text.adj=c(1.2,1.2,lwd=2))

acc.perf<- performance(pred2, measure ="acc")
plot(acc.perf)







# mydata <- read.csv(file.choose(), header=TRUE) 
 # ROC curve and assessment of my prediction 
# 
# plot(0:1, 0:1, type="n", xlab="False positive rate", ylab="True positive 
# rate") 
# abline(0, 1, col="red") 
# 
# nsim <- 5 
# auc <- rep(NA, nsim) 
# for(i in 1:nsim) { 
#   select <- sample(nrow(mydata), round(nrow(mydata)*0.7)) 
#   data70 <- mydata[select, ]  # train 
#   data30 <- mydata[-select, ]  # test 
#   temp.glm <- glm(Death ~ Temperature, data=data70, family=binomial) 
#   pred <- prediction(data30$pred, data30$Death) 
#   perf <- performance(pred, "tpr", "fpr") 
#   plot(perf, add=TRUE) 
#   auc[i] <- attributes(performance(pred, "auc"))$y.values[[1]] # area under 
#   the ROC 
# } 
# auc 

