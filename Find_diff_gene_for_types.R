###
Path_RPKM<-'/media/hufuyan/Data/networkproject2017/cancer/RESM/'

LUAD_RSEM_best_tumor<-read.table(paste(Path_RPKM,"LUAD_RSEM_best_tumor.txt",sep=''),head=T,check.names=FALSE,row.names=1)
LUAD_RSEM_best_tumor<-as.matrix(LUAD_RSEM_best_tumor)
dim(LUAD_RSEM_best_tumor)

Path_cluster<-'/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/subtypes/'


SOM_cluster_result<-read.table(paste(Path_cluster,'SOM_cluster_result.txt',sep=''),check.names=FALSE,row.names=1)
SOM_cluster_result<-as.matrix(SOM_cluster_result)
dim(SOM_cluster_result)


dir_name<-paste(Path_cluster,'diff_results/',sep='')
dir.create(dir_name)

selected_100<-c()
selected_400<-c()
count_s<-0
for (i in 1:10){
  t<-which(SOM_cluster_result==i)
  cluster_exp<-LUAD_RSEM_best_tumor[,t]
  print (dim(cluster_exp))
  othercluster_exp<-LUAD_RSEM_best_tumor[,-t]
  dim(othercluster_exp)
  p_value<-c()
  p_value_diff<-c()
  count<-0
  gene_diff<-c()
  for (j in 1:nrow(LUAD_RSEM_best_tumor)){
    results=t.test(cluster_exp[j,],othercluster_exp[j,],paired =FALSE)
    p_value[j]<-results$p.value
  }
  pvalue_sorted<-p_value[order(p_value)]
  pvalue_sorted_gene<-rownames(LUAD_RSEM_best_tumor)[order(p_value)]
  #q_value<-p.adjust(p_value,method='bonferroni')
  for (t in 1:length(pvalue_sorted)){
    if (pvalue_sorted[t]<0.025){
      count<-count+1
      p_value_diff[count]<-pvalue_sorted[t]
      gene_diff[count]<-pvalue_sorted_gene[t]
    }
    if (t<=100){
      count_s<-count_s+1
      selected_100[count_s]<-pvalue_sorted_gene[t]
      selected_400[count_s]<-pvalue_sorted_gene[t]
    }
  }
  
  print (count)
  p_value_diff<-as.matrix(p_value_diff)
  rownames(p_value_diff)<-gene_diff
  
  write.csv(p_value_diff,row.names=TRUE,file=paste(dir_name,'diff_gene_cluster',i,".csv",sep=''))
  
}

print (length(selected_100))
selected_100_unique<-unique(selected_100)
print (length(selected_100_unique))
write.csv(selected_100_unique,row.names=TRUE,file=paste(dir_name,'selected_100_unique.csv',sep=''))

########################
######heatmap
##using top100 genes with the lowest p-value for each type


tumor_list<-c(1:10)
colour_list<-c('green','blue','red','yellow')

sample_list_TCGA<-colnames(LUAD_RSEM_best_tumor)
sample_list_TCGA_new<-c()

for (i in 1:length(sample_list_TCGA)){
  bb<-strsplit(sample_list_TCGA[i],'[*]')[[1]]
  cc<-strsplit(bb[1],'[-]')[[1]]
  sample_list_TCGA_new[i]<-paste(cc[1],cc[2],cc[3],cc[4],sep='-')
}

colnames(LUAD_RSEM_best_tumor)<-sample_list_TCGA_new
####
############################
sample_list_cluster<-rownames(SOM_cluster_result)
sample_list_cluster_new<-c()

for (i in 1:length(sample_list_cluster)){
  bb<-strsplit(sample_list_cluster[i],'[*]')[[1]]
  cc<-strsplit(bb[1],'[-]')[[1]]
  sample_list_cluster_new[i]<-paste(cc[1],cc[2],cc[3],cc[4],sep='-')
}

rownames(SOM_cluster_result)<-sample_list_cluster_new


############################
library(pheatmap)


draw_heatmap_subtypes<-function(input_list,exp_all,input_clustering){ 
  exp_input<-matrix(0,length(input_list),ncol(exp_all))
  for (i in 1:length(input_list)){
    #print (i)
    t1<-which(rownames(exp_all)==input_list[i]) 
    #print (t1)
    exp_input[i,]<-exp_all[t1,]
  }
  
  rownames(exp_input)<-paste("Gene", 1:400, sep = "")
  exp_input<-as.matrix(exp_input)
  colnames(exp_input)<-colnames(exp_all)
  dim(exp_input)
  
  
  group_col<-c()
  sample<-c()
  for (i in 1:ncol(exp_input)){
    t<-which(colnames(exp_input)[i]==rownames(input_clustering))
    print (t)
    group_col[i]<-colour_list[input_clustering[t]]
    sample[i]<-tumor_list[input_clustering[t]]
  }
  length(group_col)
  # colnames(exp_input)<-sample
  sample<-as.matrix(sample)
  
  rownames(sample)<-colnames(exp_input)
  colnames(sample)<-'Group'
  sample<-as.data.frame(sample)
  
  annotation_row = data.frame(
    GeneClass = factor(rep(c("Geneset1", "Geneset2", "Geneset3", "Geneset4"), c(100,100,100,100)))
  )
  rownames(annotation_row) = paste("Gene", 1:400, sep = "")
  group_col<-list(Group=c(Subtype1="green",Subtype2="blue",Subtype3="red",Subtype4="yellow"),GeneClass = c(Geneset1='forestgreen',Geneset2='#afdfe4',Geneset3='lightpink2',Geneset4='#B2DF8A'))
  

  plot_name=paste(dir_name,'subtypes_heatmap.pdf',sep='')
  pdf(plot_name);
  
  exp<-scale(t(log10(exp_input)))
  #breaksList=c(min(exp),seq(-2,2,by=1),max(exp))
  #hv<-pheatmap(log10(exp_input),cluster_cols=FALSE,cluster_rows=FALSE, scale='row',clustering_method='ward.D2',legend=TRUE,breaks=breaksList,ColSideColors=group_col,col=colorRampPalette(c("deepskyblue","black","#ffe600"))(length(breaksList)-1),show_colnames=F,labels_row='400 genes',border_color=FALSE,fontsize=10,fontsize_row=4,fontsize_col=4,annotation_col=sample,annotation_row = annotation_row,annotation_legend=TRUE,annotation_colors=group_col)
  
  bk1<-seq(-2,2,by=1)
  bk2<-c(min(exp),-1.999999)
  bk3<-c(2.001,max(exp))
  breaksList=c(bk2,bk1,bk3)
  my_palette<-c(colorRampPalette(c("deepskyblue"))(length(bk2)),colorRampPalette(c("deepskyblue","black","#ffe600"))(length(bk1)-1),colorRampPalette(c("#ffe600"))(length(bk3)))
  hv<-pheatmap(log10(exp_input),cluster_cols=FALSE,cluster_rows=FALSE, scale='row',clustering_method='ward.D2',legend=TRUE,breaks=breaksList,ColSideColors=group_col,col=my_palette,show_colnames=F,labels_row='400 genes',border_color=FALSE,fontsize=10,fontsize_row=4,fontsize_col=4,annotation_col=sample,annotation_row = annotation_row,annotation_legend=TRUE,annotation_colors=group_col)
  
  dev.off();
  return(hv)
}




dim(SOM_cluster_result)
SOM_cluster_result_ordered<-SOM_cluster_result[order(SOM_cluster_result)]
SOM_cluster_result_ordered<-as.matrix(SOM_cluster_result_ordered)
rownames(SOM_cluster_result_ordered)<-rownames(SOM_cluster_result)[order(SOM_cluster_result)]
dim(SOM_cluster_result_ordered)
length(selected_400)


LUAD_RSEM_best_tumor_new<-LUAD_RSEM_best_tumor[,order(SOM_cluster_result)]
colnames(LUAD_RSEM_best_tumor_new)<-rownames(SOM_cluster_result_ordered)
dim(LUAD_RSEM_best_tumor_new)


input_list<-selected_400
exp_all<-LUAD_RSEM_best_tumor_new
input_clustering<-SOM_cluster_result_ordered

draw_heatmap_subtypes(selected_400,LUAD_RSEM_best_tumor_new,SOM_cluster_result_ordered)

exp<-scale(t(log10(LUAD_RSEM_best_tumor_new)))
min(exp)
median(exp)
max(exp)
plot_name=paste(dir_name,'subtypes_hist.pdf',sep='')
pdf(plot_name);

hist(exp)
dev.off()
getwd()


setwd('/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/subtypes/diff_results/')
library(VennDiagram)
length(diff_gene_cluster2)
pdf(paste(dir_name,'_venn.tiff',sep=''))
venn.diagram(list(subtype1=selected_400[1:100],subtype2=selected_400[101:200],subtype3=selected_400[201:300],subtype4=selected_400[301:400]), fill=c('green','blue','red','yellow'), alpha=c(0.5,0.5,0.5,0.5), cex=0.8, cat.fontface=4, fontfamily=3, filename="all_interactions.tiff")
dev.off()