setwd('/Users/ja/Documents/Computational Methods/Project/Model')

# Loading packages
library(caTools)
library(e1071)
library(randomForest)
library(caret)
library(glmnet)
library(pROC)
library(randomForest)

#############################
## Loading + Prepping Data
#############################
load("~/Documents/Computational Methods/Project/Model/all_tcrd_05122017.Rda")
all_gold <- read.csv('~/Documents/Computational Methods/Project/GoldStandards2.csv')
# Check for NA in dataframe
# any(is.na(all_tcrd))
# Check column names 
# check <- as.data.frame(colnames(all_tcrd))

cols_to_subtract<-c(which(colnames(all_tcrd)=="tdl"),
                    which(colnames(all_tcrd)=="idgfam"),
                    which(colnames(all_tcrd)=="predict"))
all_tcrd<-data.matrix(all_tcrd[,-cols_to_subtract])

# Colsums - removed features with <20 non-zeros values (out of 1129 samples)
temp<-all_tcrd[,colSums(all_tcrd) > 20]

# Merge the gold standard columns using IDs
temp<-merge(temp,all_gold[,c('id','y','tdl','idgfam')],by.x='target_id',by.y='id')

# Exclude the target_id column for training
id_index <- c(which(colnames(temp)=="target_id"))
temp<-temp[,c(-id_index)]
temp <- subset(temp,select=-c(tdl,idgfam))

write.csv(temp,"NewXy.csv")

