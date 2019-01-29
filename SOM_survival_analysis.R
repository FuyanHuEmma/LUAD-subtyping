
##https://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/
#som
setwd('D:/English_paper/')

Path_input<-'D:/English_paper/'



LUAD_RSEM_best_tumor<-read.table(paste(Path_input,"LUAD_RSEM_best_tumor.txt",sep=''),sep='\t',check.names=FALSE,row.names=1,header=1,fill=TRUE,quote = "")
LUAD_RSEM_best_tumor<-as.matrix(LUAD_RSEM_best_tumor)
dim(LUAD_RSEM_best_tumor)

select_dataset<-read.table('select_dataset.txt',row.names=1,header=1)
select_dataset<-as.matrix(select_dataset)
dim(select_dataset)



som_dataset<-t(select_dataset)
standardization<-function(x){x/sqrt(sum(x*x))}
sta_som_dataset<-t(apply(som_dataset,1,standardization))
library(kohonen)
sommodel<-som(sta_som_dataset,
              grid=somgrid(12,13,topo = "rectangular",neighbourhood.fct = "gaussian"),
              rlen=200,alpha=c(0.05,0.01)) 
summary(sommodel) 
sommodel$unit.classif 
sommodel$codes 
sommodel$distances 
sommodel$changes 

pdf('changes.pdf')
plot(sommodel, type="changes")
dev.off()
pdf('countplot.pdf')
plot(sommodel, type="counts", shape = "straight")#Empty units are depicted in gray
dev.off()
## show both sets of codebook vectors in the map
pdf('codeplots.pdf')
plot(sommodel, type="codes") 
dev.off()
par(mfrow = c(1,1))
pdf('qualityplot.pdf')
similarities <- plot(sommodel, type="quality", palette.name = terrain.colors)
dev.off()
pdf('mappingplot.pdf')
plot(sommodel, type="mapping",
     labels = c(rep(1,512),rep(2,58)), col = c(rep(1,512),rep(2,58)),
     main = "mapping plot")
dev.off()

  save(sommodel, file = 'som.Rdata')
  
  

pdf('neighbourdistance.pdf')
Umat <- plot(sommodel, type="dist.neighbours", main = "SOM neighbour distances")
## use hierarchical clustering to cluster the codebook vectors
som.hc <- cutree(hclust(object.distances(sommodel, "codes")), 4)


add.cluster.boundaries(sommodel, som.hc)
dev.off()

pdf('mappingplot_withcluster.pdf')
plot(sommodel, type="mapping",
     labels = c(rep(1,512),rep(2,58)), col = c(rep(1,512),rep(2,58)),
     main = "mapping plot")
som.hc <- cutree(hclust(object.distances(sommodel, "codes")), 4)
add.cluster.boundaries(sommodel, som.hc)
dev.off()

write.table(cutree(hclust(object.distances(sommodel, "codes")), 4),file='hcluster.txt')



######
#mydata <- sommodel$codes[[1]]
#wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var)) 
#for (i in 2:15) {
#  wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
#}
#plot(wss)


##
class_target<-rep(0,570)
for(i in 1:570){
  for(j in 1:156){
    if(sommodel$unit.classif[i]==j) class_target[i]=som.hc[j]
  }
}
#class_target
write.table(class_target,file = "class_target.txt")

lungsubtype<-class_target[1:512] #·Î°©Ñù±¾
lungsubtype<-as.matrix(lungsubtype)
rownames(lungsubtype)<-colnames(LUAD_RSEM_best_tumor)



write.csv(lungsubtype,row.names=T,file = "cluster_result.csv")


####
LUAD_clin_merged<-read.table(paste(Path_input,"LUAD_clin_merged.txt",sep=''),sep='\t',check.names=FALSE,row.names=1,fill=TRUE,quote = "")
LUAD_clin_merged<-as.matrix(LUAD_clin_merged)
dim(LUAD_clin_merged)


####
event<-c()
time_to_event<-c()
sample_event<-c()
count<-0
t1<-which(rownames(LUAD_clin_merged)=='patient.vital_status')
t2<-which(rownames(LUAD_clin_merged)=='patient.days_to_last_followup')
t3<-which(rownames(LUAD_clin_merged)=='patient.days_to_death')

s1<-which(rownames(LUAD_clin_merged)=='patient.bcr_patient_barcode')

for (i in 1:ncol(LUAD_clin_merged)){
#print (i)
if (LUAD_clin_merged[t1,i]=="alive" && !is.na(LUAD_clin_merged[t2,i])){
count<-count+1
event[count]<-0

time_to_event[count]<-as.numeric(LUAD_clin_merged[t2,i])/30
sample_event[count]<-toupper(LUAD_clin_merged[s1,i])
}
if (LUAD_clin_merged[t1,i]=="dead" && !is.na(LUAD_clin_merged[t3,i])){
count<-count+1
event[count]<-1
#print (time_to_event)
time_to_event[count]<-as.numeric(LUAD_clin_merged[t3,i])/30

sample_event[count]<-toupper(LUAD_clin_merged[s1,i])

}
}

print (count)
t<-which(event==0)
mean(time_to_event[t])/12

####
cluster_result<-read.csv(paste(Path_cluster,"cluster_result.csv",sep=''),head=T,row.names=1,check.names=FALSE)
cluster_result<-as.matrix( cluster_result)
dim(cluster_result)
sample_name<-rownames(cluster_result)


sample_list_TCGA_new<-c()  

for (i in 1:length(sample_name)){
bb<-strsplit(sample_name[i],'[-]')[[1]]
sample_list_TCGA_new[i]<-paste(bb[1],bb[2],bb[3],sep='-')
}

sample_list_clinical<-sample_event   
length(sample_list_clinical)


common_sample<-intersect(sample_list_TCGA_new,sample_list_clinical) 
print (paste('LUAD: the no.of common samples for clinical data and TCGA is ',length(common_sample),sep=''))

sample_label_TCGA<-c()
sample_label_clinical<-c()
for (i in 1:length(common_sample)){
t1=which(sample_list_TCGA_new==common_sample[i])
t2=which(sample_list_clinical==common_sample[i])
sample_label_TCGA[i]<-t1[1]
sample_label_clinical[i]<-t2[1]
}


event<-event[sample_label_clinical]
time_to_event<-time_to_event[sample_label_clinical]
length(event)
length(time_to_event)

library(survival)

y<-Surv(time_to_event,event)

#########################################
####divide samples based on the clustring results using hubgene exp
sample_divide<-c()
for (i in 1:nrow(cluster_result)){
if (cluster_result[i]==1){
sample_divide[i]<-1
} else if (cluster_result[i]==2){
sample_divide[i]<-2
} else if (cluster_result[i]==3){
sample_divide[i]<-3
} else {
sample_divide[i]<-4
}
}

###


sample_divide_clin<-c()
for (i in 1:length(common_sample)){
t1=which(sample_list_TCGA_new==common_sample[i])
#print (t1)
sample_divide_clin[i]<-sample_divide[t1]
}

diff0<-survdiff(y~sample_divide_clin)
pvalue<- 1 - pchisq(diff0$chisq, length(diff0$n) - 1)
pvalue

pdf('survival_plot_4.pdf')

fit<-survfit(y~sample_divide_clin)
plot(fit,xlab='Time in months',ylab='Survival probality',
col=c('green','blue','red','yellow'),
lty = c('solid', 'solid','solid', 'solid'),lwd=c(2,2,2,2))

legend('topright', c('subtype 1 (n=164)','subtype 2 (n=206)','subtype 3 (n=113)','subtype 4 (n=9)'),
col=c('green','blue','red','yellow'),
lty = c('solid', 'solid','solid', 'solid'),lwd=c(2,2,2,2), bty="n")
legend('bottomleft',paste('p-value=',0.0047,sep=''), bty="n")

dev.off()
###################
###################
###################

t4<-c()
count<-0
for (i in 1:nrow(cluster_result)){
if (cluster_result[i]==1){
sample_divide[i]<-1
} else if (cluster_result[i]==2){
sample_divide[i]<-2
} else if (cluster_result[i]==3){
sample_divide[i]<-3
} else {
count<-count+1
sample_divide[i]<-4
t4[count]<-i
}
}

sample_name<-rownames(cluster_result)[-t4]


sample_list_TCGA_new<-c() 

for (i in 1:length(sample_name)){
bb<-strsplit(sample_name[i],'[-]')[[1]]
sample_list_TCGA_new[i]<-paste(bb[1],bb[2],bb[3],sep='-')
}

sample_list_clinical<-sample_event   
length(sample_list_clinical)


common_sample2<-intersect(sample_list_TCGA_new,sample_list_clinical) 
print (paste('LUAD: the no.of common samples for clinical data and TCGA is ',length(common_sample2),sep=''))



sample_label_TCGA<-c()
sample_label_clinical<-c()
for (i in 1:length(common_sample2)){
t1=which(sample_list_TCGA_new==common_sample2[i])
t2=which(sample_list_clinical==common_sample2[i])
sample_label_TCGA[i]<-t1[1]
sample_label_clinical[i]<-t2[1]
}


event<-event[sample_label_clinical]
time_to_event<-time_to_event[sample_label_clinical]
length(event)
length(time_to_event)

library(survival)

y<-Surv(time_to_event,event)


sample_divide_clin<-c()
sample_divide2<-sample_divide[-t4]

for (i in 1:length(common_sample2)){
t1=which(sample_list_TCGA_new==common_sample2[i])
#print (t1)
sample_divide_clin[i]<-sample_divide2[t1]
}

diff0<-survdiff(y~sample_divide_clin)
pvalue<- 1 - pchisq(diff0$chisq, length(diff0$n) - 1)
pvalue

pdf('survival_plot_3.pdf')

fit<-survfit(y~sample_divide_clin)
plot(fit,xlab='Time in months',ylab='Survival probality',
col=c('green','blue','red','yellow'),
lty = c('solid', 'solid','solid', 'solid'),lwd=c(2,2,2,2))

legend('topright', c('subtype1(n=164)','subtype2(n=206)','subtype3(n=113)'),
col=c('green','blue','red'),
lty = c('solid', 'solid','solid', 'solid'),lwd=c(2,2,2,2))

legend('bottomleft',paste('pvalue=',0.0019,sep=''))

dev.off()

#############
library(survminer)

plot_name= paste(Path_clinical_out,type,'_survival_plot.pdf',sep='')
pdf(plot_name)
survp1<-ggsurvplot(fit,data=all_inf,pval = TRUE, break.time.by=12,gtheme = theme_minimal(),xlab='Time in months',legend.labs =c('Subtyping I','Subtyping II','Subtyping III','Subtyping IV','Subtyping V','Subtyping VI','Subtyping VII','Subtyping VIII','9'),palette= c('green','blue','red','yellow','hotpink','purple','tan3','salmon','tan3'))
#survp1<-ggsurvplot(fit,data=all_inf,pval = TRUE, break.time.by=12,gtheme = theme_minimal(),xlab='Time in months',legend.labs =c('Subtyping I','Subtyping II','Subtyping III','Subtyping IV','Subtyping V','Subtyping VI','Subtyping VII','Subtyping VIII','9','10'),palette= c('green','blue','red','yellow','hotpink','purple','tan3','salmon','tan3','salmon'))
print(survp1$plot, newpage = FALSE)

dev.off();

