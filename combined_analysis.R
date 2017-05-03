combined_IC_GPCR<-bind_rows(IC_results,gcpr_results)
combined_IC_GPCR[is.na(combined_IC_GPCR)]<-0

# I am scared to lose the original combined_IC_GPCR
# get rid of tdl_class and target_id
combined_IC_GPCR_reduced<-combined_IC_GPCR[,-c(which(colnames(combined_IC_GPCR)=="tdl_class"))]
all_gold<-rbind(gold_GPCR,gold_IC)
library(dplyr)
combined_IC_GPCR_reduced<-combined_IC_GPCR_reduced %>% filter(target_id %in% all_gold$id)
#############################
## Edit from here for new training
#############################

temp<-combined_IC_GPCR_reduced[,colSums(combined_IC_GPCR_reduced) > 10]
# Merge the gold standard 'y' column using IDs
temp<-merge(temp,all_gold[,c('id','y')],by.x='target_id',by.y='id')
# Exclude the target_id column for training
temp<-temp[,c(-1)]

# Split the sample at even distribution using our y gold standard label
library(caTools)
spl<-sample.split(temp$y,SplitRatio = 0.8)
combined_IC_GPCR_results_reduced_train<-temp[spl,]
combined_IC_GPCR_results_reduced_test<-temp[!spl,]

# Make sure we got rid of class
library(e1071)
library(randomForest)
library(caret)

combined_IC_GPCR_results_reduced_train$y<-as.factor(combined_IC_GPCR_results_reduced_train$y)
# model <- naiveBayes(y~.,data=combined_IC_GPCR_results_reduced_train)
model<-randomForest(y~.,data=combined_IC_GPCR_results_reduced_train)
pred <- predict(model,combined_IC_GPCR_results_reduced_test[,c(-which(colnames(combined_IC_GPCR_results_reduced_test)=="y"))])
table(pred, combined_IC_GPCR_results_reduced_test$y)

## For randomForest
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)

###########################
# GLMNET
###########################
y <- combined_IC_GPCR_results_reduced_train$y
label_col<-which(colnames(combined_IC_GPCR_results_reduced_train)=="y")
x        <- as.matrix(data.frame(combined_IC_GPCR_results_reduced_train[,c(-label_col)],
                                 y))
fit.lasso<-glmnet(x,
                  y,
              alpha=1,
              family='multinomial')
# # Plot it what the hell is that!
# plot(model,xvar="lambda")
fit.lasso.cv<-cv.glmnet(x,
                        y,
                     alpha=1,
                     nfolds=5, 
                     family='multinomial',
                     type.measure = "class")
plot(fit.lasso, xvar="lambda")
fit.lasso$lambda.min

glm_test_data<-as.matrix(data.frame(combined_IC_GPCR_results_reduced_test[,c(-which(colnames(combined_IC_GPCR_results_reduced_test)=="y"))]))
pred<-predict(cv.glmmod,newx=glm_test_data,x="lambda.min",type="class")
table(pred, combined_IC_GPCR_results_reduced_test$y)
