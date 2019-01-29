setwd('/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/PCA/')
Path_RPKM<-'/media/hufuyan/Data/networkproject2017/cancer/RESM/'

LUAD_RSEM_best_tumor<-read.table(paste(Path_RPKM,"LUAD_RSEM_best_tumor.txt",sep=''),head=T,check.names=FALSE,row.names=1)
LUAD_RSEM_best_tumor<-as.matrix(LUAD_RSEM_best_tumor)
dim(LUAD_RSEM_best_tumor)

LUAD_RSEM_best_normal<-read.table(paste(Path_RPKM,"LUAD_RSEM_best_normal.txt",sep=''),head=T,check.names=FALSE,row.names=1)
LUAD_RSEM_best_normal<-as.matrix(LUAD_RSEM_best_normal)
dim(LUAD_RSEM_best_normal)


###selected genes
Path_input<-'/media/hufuyan/Data/networkproject2017/cancer/script_subtypes/'
specialgene30<-read.table(paste(Path_input,"specialgene30.txt",sep=''),check.names=FALSE)
specialgene30<-as.matrix(specialgene30)
dim(specialgene30)

selected_gene_exp_tumor<-matrix(0,30,ncol(LUAD_RSEM_best_tumor))
for (i in 1:30){
  t<-which(specialgene30[i]==rownames(LUAD_RSEM_best_tumor))
  print (t)
  selected_gene_exp_tumor[i,]=LUAD_RSEM_best_tumor[t,]
}

selected_gene_exp_normal<-matrix(0,30,ncol(LUAD_RSEM_best_normal))
for (i in 1:30){
  t1<-which(specialgene30[i]==rownames(LUAD_RSEM_best_normal))
  print (t1)
  selected_gene_exp_normal[i,]=LUAD_RSEM_best_normal[t1,]
}



selected_gene_exp_all<-cbind(selected_gene_exp_tumor,selected_gene_exp_normal)
selected_gene_exp_all<-as.matrix(selected_gene_exp_all)
rownames(selected_gene_exp_all)<-rownames(specialgene30)
colnames(selected_gene_exp_all)<-c(rep('Tumor',ncol(LUAD_RSEM_best_tumor)),rep('Normal',ncol(LUAD_RSEM_best_normal)))
dim(selected_gene_exp_all)

###PCA
tumor_list<-c('Tumor','Normal')
colour_list<-c('red','green')
group_col<-c(rep('red',ncol(LUAD_RSEM_best_tumor)),rep('green',ncol(LUAD_RSEM_best_normal)))

library(scatterplot3d)
pdf('PCA_selected_gene30.pdf');
PCA_selected_gene30<-princomp(scale(t(selected_gene_exp_all)),cor=FALSE)
s3d<-scatterplot3d(x = PCA_selected_gene30$scores[,1], y =PCA_selected_gene30$scores[,2], z = PCA_selected_gene30$scores[,3],color=group_col,pch=20,xlab ='PC1',ylab='PC3',zlab='PC2')
legend(s3d$xyz.convert(10, 10,-4),legend = tumor_list,col = colour_list,pch=20, bty="n")
dev.off();


sink("log_PCA.txt")
summary(PCA_selected_gene30,loadings=TRUE)
sink()

###画主成分的碎石图
pdf('PCA_Ygene_screeplot.pdf');
screeplot(PCA_selected_gene30,type="lines")
dev.off();




factor(rownames(LUAD_RSEM_best_tumor)==rownames(LUAD_RSEM_best_normal))





####
###for wilcoxon selected top 30 genes
Path_wtest<-'/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/fd_pvalue/'
wtest_gene30<-read.csv(paste(Path_wtest,"wilcoxon_selected_top30_genes.csv",sep=''),check.names=FALSE,row.names=1,head=T)
wtest_gene30<-as.matrix(wtest_gene30)
dim(wtest_gene30)


wtest_selected_gene_exp_tumor<-matrix(0,30,ncol(LUAD_RSEM_best_tumor))
for (i in 1:30){
  t<-which(wtest_gene30[i]==rownames(LUAD_RSEM_best_tumor))
 # print (t)
  wtest_selected_gene_exp_tumor[i,]=LUAD_RSEM_best_tumor[t,]
}

wtest_selected_gene_exp_normal<-matrix(0,30,ncol(LUAD_RSEM_best_normal))
for (i in 1:30){
  t1<-which(wtest_gene30[i]==rownames(LUAD_RSEM_best_normal))
 # print (t1)
  wtest_selected_gene_exp_normal[i,]=LUAD_RSEM_best_normal[t1,]
}


wtest_selected_gene_exp_all<-cbind(wtest_selected_gene_exp_tumor,wtest_selected_gene_exp_normal)
wtest_selected_gene_exp_all<-as.matrix(wtest_selected_gene_exp_all)
rownames(wtest_selected_gene_exp_all)<-rownames(wtest_gene30)
colnames(wtest_selected_gene_exp_all)<-c(rep('Tumor',ncol(LUAD_RSEM_best_tumor)),rep('Normal',ncol(LUAD_RSEM_best_normal)))
dim(wtest_selected_gene_exp_all)

###PCA
tumor_list<-c('Tumor','Normal')
colour_list<-c('red','green')
group_col<-c(rep('red',ncol(LUAD_RSEM_best_tumor)),rep('green',ncol(LUAD_RSEM_best_normal)))

library(scatterplot3d)
pdf('PCA_wtest_selected_gene30.pdf');
PCA_wtest_selected_gene30<-princomp(scale(t(wtest_selected_gene_exp_all)),cor=FALSE)
s3d<-scatterplot3d(x = PCA_wtest_selected_gene30$scores[,3], y =PCA_wtest_selected_gene30$scores[,1], z = PCA_wtest_selected_gene30$scores[,2],color=group_col,pch=20,xlab ='PC3',ylab='PC1',zlab='PC2')
legend(s3d$xyz.convert(10, 10,-4),legend = tumor_list,col = colour_list,pch=20, bty="n")
dev.off();


sink("log_PCA_wtest.txt")
summary(PCA_wtest_selected_gene30,loadings=TRUE)
sink()

###画主成分的碎石图
pdf('PCA_30genes_screeplot_wtest.pdf');
screeplot(PCA_wtest_selected_gene30,type="lines")
dev.off();

