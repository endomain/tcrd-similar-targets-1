Rprof ( tf <- "log.log",  memory.profiling = TRUE )
all_tcrd<-bind_rows(cum_results,combined_IC_GPCR)
Rprof ( NULL ) ; print ( summaryRprof ( tf )  )

all_tcrd[is.na(all_tcrd)]<-0

# combine all the gold
gold_GPCR<-read.csv('toxic_targets/Gold_GPCR.csv')
gold_Kinase<-read.csv('toxic_targets/gold_Kinase.csv')
gold_IC<-read.csv('toxic_targets/gold_IC.csv')
all_gold<-bind_rows(gold_IC,gold_GPCR)
all_gold<-bind_rows(all_gold,gold_Kinase[,c(-1)])

all_tcrd_reduced<-all_tcrd[,-c(which(colnames(all_tcrd)=="tdl_class"))]
# Select for training
all_tcrd_reduced<-all_tcrd_reduced %>% filter(target_id %in% all_gold$id)

#############################################
# DO NOT EDIT ABOVE TAKES FOREVER!
#############################################

########################
# Edit below
#######################

# get rid of rare features
temp<-all_tcrd_reduced[,colSums(all_tcrd_reduced) > 10]
# merge our gold standards
temp<-merge(temp,all_gold[,c('id','y')],by.x='target_id',by.y='id')
# Exclude the target_id column for training
temp<-temp[,c(-1)]

# Split the sample at even distribution using our y gold standard label
library(caTools)
spl<-sample.split(temp$y,SplitRatio = 0.7)
all_tcrd_reduced_train<-temp[spl,]
all_tcrd_reduced_test<-temp[!spl,]

library(e1071)
library(randomForest)
library(caret)

# Make sure y is a factor
all_tcrd_reduced_train$y<-as.factor(all_tcrd_reduced_train$y)

# model <- naiveBayes(y~.,data=combined_IC_GPCR_results_reduced_train)
model<-randomForest(y~.,data=all_tcrd_reduced_train)
pred <- predict(model,all_tcrd_reduced_test[,c(-which(colnames(all_tcrd_reduced_test)=="y"))])
table(pred, all_tcrd_reduced_test$y)
