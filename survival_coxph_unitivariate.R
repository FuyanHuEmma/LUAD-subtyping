path_1<-'/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/subtypes/'
Path_clinical_out<-paste(path_1,'/clinical_analysis/',sep='')
dir.create(Path_clinical_out)

setwd('/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/subtypes/diff_results/')

###

Path_input<-'/media/hufuyan/Data/networkproject2017/cancer/clinical_data/'

LUAD_clin_merged<-read.table(paste(Path_input,"LUAD_clin_merged.txt",sep=''),sep='\t',check.names=FALSE,row.names=1,fill=TRUE,quote = "")
LUAD_clin_merged<-as.matrix(LUAD_clin_merged)
dim(LUAD_clin_merged)


##compare samples
Path_cluster<-'/media/hufuyan/Data/networkproject2017/cancer_result/LUAD_subtypes/subtypes/'


SOM_cluster_result<-read.table(paste(Path_cluster,'SOM_cluster_result.txt',sep=''),check.names=FALSE,row.names=1)
SOM_cluster_result<-as.matrix(SOM_cluster_result)
dim(SOM_cluster_result)
###


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


t<-which(event==0)
mean(time_to_event[t])/12
length(event)
#########


#######
fixed_effects<-c("patient.gender","patient.stage_event.pathologic_stage","patient.tobacco_smoking_history" ,"patient.age_at_initial_pathologic_diagnosis")


##filter out the one without full alive/death information

t_gender<-which(rownames(LUAD_clin_merged)=='patient.gender')
t_stage<-which(rownames(LUAD_clin_merged)=='patient.stage_event.pathologic_stage')
t_smoking<-which(rownames(LUAD_clin_merged)=='patient.tobacco_smoking_history')
#t_retrospective<-which(rownames(LUAD_clin_merged)=='patient.tissue_retrospective_collection_indicator')
#t_race<-which(rownames(LUAD_clin_merged)=='patient.race_list.race')
#t_initial<-which(rownames(LUAD_clin_merged)=='patient.year_of_initial_pathologic_diagnosis')
t_age<-which(rownames(LUAD_clin_merged)=='patient.age_at_initial_pathologic_diagnosis')

feature_all<-c(t_gender,t_stage,t_smoking,t_age)
name_feature<-c('gender','stage','smoking_indicator','age')

t_status<-which(rownames(LUAD_clin_merged)=='patient.vital_status')
t_followup<-which(rownames(LUAD_clin_merged)=='patient.days_to_last_followup')
t_days<-which(rownames(LUAD_clin_merged)=='patient.days_to_death')

t_patient<-which(rownames(LUAD_clin_merged)=='patient.bcr_patient_barcode')


library(survival)

for (k in 1:length(feature_all)){
  k<-2
  feature<-feature_all[k]
  deleted_col<-c()
  deleted_col_feature<-c()
  count<-0
  count_feature<-0
  for (i in 1:ncol(LUAD_clin_merged)){
    t<-is.na(LUAD_clin_merged[feature,i])
    
    a<-which(t=='TRUE')
    if (length(a)>0){
      print (1)
      count<-count+1
      count_feature<-count_feature+1
      deleted_col[count]<-i
      deleted_col_feature[count_feature]<-i
    }
    if (LUAD_clin_merged[t_status,i]=="alive"){ 
      if (is.na(LUAD_clin_merged[t_followup,i])){
        count<-count+1
        deleted_col[count]<-i
        print (2)
      }
    }
    if (LUAD_clin_merged[t_status,i]=="dead"){
      if (is.na(LUAD_clin_merged[t_days,i])){
        count<-count+1
        print (3)
        print (LUAD_clin_merged[t_patient,i])
        deleted_col[count]<-i
      }
    }
  }
  deleted_col<-unique(deleted_col)
  print (paste(length(deleted_col),'sampels need to be delected due to the NA'))
  print (paste(length(deleted_col_feature),'sampels need to be delected due to the NA of feature'))
  
  all_features_label<-c(t_patient,t_status,t_followup,t_days,feature)
  tumor_clinical_data_all_features<-LUAD_clin_merged[all_features_label,-deleted_col]
  dim(tumor_clinical_data_all_features)
  
  t_patient_new<-which(rownames(tumor_clinical_data_all_features)=='patient.bcr_patient_barcode')
  sample_list_TCGA<-rownames(SOM_cluster_result)
  length(sample_list_TCGA)
  
  sample_list_TCGA_new<-c()
  
  for (i in 1:length(sample_list_TCGA)){
    bb<-strsplit(sample_list_TCGA[i],'[*]')[[1]]
    cc<-strsplit(bb[1],'[-]')[[1]]
    sample_list_TCGA_new[i]<-paste(cc[1],cc[2],cc[3],sep='-')
  }
  
  sample_list_clinical<-toupper(tumor_clinical_data_all_features[t_patient_new,])
  length(sample_list_clinical)
  
  
  common_sample<-intersect(sample_list_TCGA_new,sample_list_clinical)
  print (paste('LUAD',': the no.of common samples for clinical data and TCGA exp is ',length(common_sample),sep=''))
  
  sample_label_TCGA<-c()
  sample_label_clinical<-c()
  for (i in 1:length(common_sample)){
    t1=which(sample_list_TCGA_new==common_sample[i])
    t2=which(sample_list_clinical==common_sample[i])
    sample_label_TCGA[i]<-t1[1]
    sample_label_clinical[i]<-t2[1]
  }
  
  ####clinical data
  tumor_clinical_data_used<-tumor_clinical_data_all_features[,sample_label_clinical]
  dim(tumor_clinical_data_used)
  tumor_clinical_data_used<-data.frame(tumor_clinical_data_used, stringsAsFactors=FALSE)
  

  ######
  tumor_clinical_data_used_new<-tumor_clinical_data_used
  
  if (k=1){
    ##gender
    ###[1] "FEMALE" "MALE"  
    t6<-which(rownames(tumor_clinical_data_used)=="patient.gender")
    gender_sex0<-tumor_clinical_data_used[t6,]
    gender_sex<-c()
    for (i in 1:length(gender_sex0)){
      if (gender_sex0[i]=="female" ){
        gender_sex[i]<-0
      } else {
        gender_sex[i]<-1
      }
    }
    tumor_clinical_data_used_new[t6,]<-gender_sex
  }
  if (k=2){
    tt<-which(rownames(tumor_clinical_data_used)=="patient.stage_event.pathologic_stage")
    tumor_stage0<-as.character(toupper(tumor_clinical_data_used[tt,]))
    
    tumor_stage<-c()
    for (i in 1:length(tumor_stage0)){
      tumor_stage0[i]<-gsub('A','',tumor_stage0[i])
      tumor_stage0[i]<-gsub('B','',tumor_stage0[i])
      tumor_stage0[i]<-gsub('C','',tumor_stage0[i])
      tumor_stage0[i]<-gsub('a','',tumor_stage0[i])
      tumor_stage0[i]<-gsub('b','',tumor_stage0[i])
      tumor_stage0[i]<-gsub('c','',tumor_stage0[i])
      tumor_stage0[i]<-gsub('c','',tumor_stage0[i])
      
      if (tumor_stage0[i]=="STGE 0"){
        tumor_stage[i]<-0
      }
      if (tumor_stage0[i]=="STGE I"){
        tumor_stage[i]<-1
      }
      if (tumor_stage0[i]=="STGE II"){
        tumor_stage[i]<-2
      }
      if (tumor_stage0[i]=="STGE III"){
        tumor_stage[i]<-3
      }
      if (tumor_stage0[i]=="STGE IV"){
        tumor_stage[i]<-4
      }
      if (tumor_stage0[i]=="STGE X"){
        tumor_stage[i]<-5
      }
    }
    tumor_clinical_data_used_new[tt,]<-tumor_stage
  }
  if (k=3){
    t5<-which(rownames(tumor_clinical_data_used)=="patient.tobacco_smoking_history")
    smoking_history0<-tumor_clinical_data_used[t5,]
    smoking_history<-c()
    for (i in 1:length(smoking_history0)){
      if (smoking_history0[i]==1 ){  ##==1 means 'Lifelong Non-smoker'
        smoking_history[i]<-0
      } else {
        smoking_history[i]<-1
      }
    }
    tumor_clinical_data_used_new[t5,]<-smoking_history
    
  }
  
  event<-c()
  time_to_event<-c()
  for (i in 1:ncol(tumor_clinical_data_used_new)){
    if (tumor_clinical_data_used_new[2,i]=="alive"){
      event[i]<-0
      time_to_event[i]<-as.numeric(as.character(tumor_clinical_data_used_new[3,i]))/30
    }
    if (tumor_clinical_data_used_new[2,i]=="dead"){
      event[i]<-1
      time_to_event[i]<-as.numeric(as.character(tumor_clinical_data_used_new[4,i]))/30
    }
  }
  tumor_clinical_data_used_new[2,]<-event
  
  
  
  tumor_clinical_data_used_new<-t(tumor_clinical_data_used_new)
  dim(tumor_clinical_data_used_new)

  sink(paste(name_feature[k],'univariate.txt',sep='_'))
  if (k==5){
    res.cox2<-coxph(Surv(time_to_event,event)~tumor_clinical_data_used_new[,5])
  } else{
    res.cox2<-coxph(Surv(time_to_event,event)~as.numeric(as.character(tumor_clinical_data_used_new[,5])))
  }
  summary(res.cox2)
  sink()
}

table(tumor_clinical_data_used_new[,5])

median(as.numeric(as.character(tumor_clinical_data_used_new[,5])))
range(as.numeric(as.character(tumor_clinical_data_used_new[,5])))

sink('cluster_univariate.txt')
res.cox2<-coxph(Surv(time_to_event,event)~factor(SOM_cluster_result[sample_label_TCGA]))
summary(res.cox2)
sink()
length(SOM_cluster_result[sample_label_TCGA])
survdiff(Surv(time_to_event,event)~factor(SOM_cluster_result[sample_label_TCGA]),rho = 0)
length(event)



t<-which(SOM_cluster_result[sample_label_TCGA]==3)
table(tumor_clinical_data_used_new[t,5])

table(tumor_clinical_data_used_new[,5])
tumor_clinical_data_used_new[t,2]

colnames(tumor_clinical_data_used_new)
tumor_clinical_data_used_new[t,]