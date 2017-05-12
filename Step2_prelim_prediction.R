######################################################
# This file is a crude 'preview' run of a basic
# prediction model.
#
# The actual analysis is done in Python,
# in folder `analysis_code/`
######################################################

library(caTools)
library(e1071)
library(randomForest)
library(caret)
######################################################
# Remove Sparse features
# Use temporary data frame to work with
######################################################

get_minimum_colSum_cols<-function(df,min){
  result<-as.vector(0)
  for(col_num in 1:length(df)){
    if(sum(df[,col_num])>=min & sum(df[,col_num])<500){
      result <- c(result, col_num)
    }
    if(col_num %% 2000 == 0){
      print(paste0(col_num / length(df)," %"))
    }
  }
  return(result)
}
col_to_keep<-get_minimum_colSum_cols(all_tcrd_reduced,20)
temp<-all_tcrd_reduced[,col_to_keep]

######################################################
# exclude certain features from tables for Experimentation
# ######################################################
# temp<-temp[,!grepl("expression",colnames(temp))]

temp$target_id<-unlist(all_tcrd_reduced_target_id)
######################################################
# Merge our gold standards list
######################################################

temp<-merge(temp,all_gold[,c('id','tdl','idgfam','y')],by.x='target_id',by.y='id')

######################################################
# Exclude the target_id column for training our data
######################################################

temp<-temp[,c(-1)]


######################################################
# Split the sample at even proportions using our 'y'
# which indicates target's association with a withdrawn drug
######################################################
spl<-sample.split(temp$y,SplitRatio = 0.8)
all_tcrd_reduced_train<-temp[spl,]
all_tcrd_reduced_test<-temp[!spl,]


# Make sure y is a factor
all_tcrd_reduced_train$y<-as.factor(all_tcrd_reduced_train$y)
model<-randomForest(y~.-tdl-idgfam,data=all_tcrd_reduced_train)
pred <- predict(model,all_tcrd_reduced_test[,c(-which(colnames(all_tcrd_reduced_test)=="y"))])
table(pred, all_tcrd_reduced_test$y)

View(varImp(model))

