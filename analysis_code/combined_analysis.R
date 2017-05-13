setwd('/Users/ja/Documents/Computational Methods/Project/Model')

# Loading packages
library(caTools)
library(e1071)
library(randomForest)
library(caret)
library(glmnet)
library(pROC)
library(randomForest)

# Reading in data from ch2 feature selection
temp <- read.csv('/Users/ja/Documents/Computational Methods/Project/Model/NewXy_chi2.csv')

# Split the gold standards sample
set.seed(1)
spl<-sample.split(temp$y,SplitRatio = 0.8)
train<-temp[spl,]
test<-temp[!spl,]

# Need to double check if there's extra rows (not all targets are gold standards targets)
#------------------------------------------------------------------------------------------------
# Set training x and y
train$y <- as.factor(train$y)
y <- train$y

label_col<-which(colnames(train)=="y")
x <- as.matrix(train[,c(-label_col)])

# set testing x and y
testx <- test[,c(-label_col)]
testx <- as.matrix(testx)

test$y <- as.factor(test$y)
testy <- test$y

###########################
# Naive Bayes
###########################
set.seed(1)
nb.fit <- naiveBayes(x,y)

# Get prediction by probability
nb.pred <- predict(nb.fit,testx,type = 'raw')
nb.pred <- nb.pred[,2]
nb.roc <-roc(testy,nb.pred)

# Get prediciton by 0-1 
#nb.pred_01 <- predict(nb.fit,testx,type = 'response')
#confusionMatrix(nb.pred,testy)

###########################
# GLMNET
###########################
# CV - Lasso
fit.cv<-cv.glmnet(x,y,alpha=1,nfolds=5,family='multinomial',type.measure = "class")
  #to plot lambda by fold
  pdf(file="CV_Lasso_misclassification_lambda.pdf")
  plot(fit.cv)
  dev.off()

# Get none-zero coefficients 
coef.cv <- coef(fit.cv)$`1`
coef.cv <- as.matrix(coef.cv)
coef.cv <- as.data.frame(coef.cv[coef.cv != 0,])
colnames(coef.cv) <- 'coefficient'
write.csv(coef.cv,"L1LogisticRegression_Coefficients3.csv")

# prediction by 0-1 response 
ypred.cv <- predict(fit.cv,testx,s="lambda.min",type="class")
confusionMatrix(ypred.cv,testy)

# prediciton by probability
ypred.cv_act <- predict(fit.cv,testx,s="lambda.min",type="response")
  # selected the predicted probability for 1 (rather than 0)
ypred.cv_act <- ypred.cv_act[,2,]

# Get ROC
# 1) get AUC based on 0-1 classification
#auc <- roc(testy,ypred.cv_act)
#testy <- as.integer(testy)

# 2) get AUC based on probability estimates
cv.lasso.roc <-roc(testy,ypred.cv_act)
# to get AUC value
auc(cv.lasso.roc)

# 3) plot individual ROC curve
#pdf(file="ROC_cvLasso.pdf")
#plot.roc(cv.lasso.roc,xlim=c(1,0),main="ROC curve of Logistic Regression (CV-Lasso)",col="blue")
#legend("topleft",legend=c("AUC=0.7713"))
#dev.off()

###########################
# Random Forest
###########################
rf_fit <- randomForest(x=x,y=y,importance = TRUE,proximity = TRUE)
# prediction by 0-1 response or probability
rf_pred <- predict(rf_fit,testx,type = 'response')
rf_predp <- predict(rf_fit,testx,type = 'prob')
  # selected the predicted probability for 1 (rather than 0)
rf_predp <- rf_predp[,2]

# Get ROC
# 1) Get absolute misclassification based on 0-1 loss
confusionMatrix(rf_pred,testy)
# 2) Get prediction accuracy based on probability
rf.roc <-roc(testy,rf_predp)
rf.auc <- auc(rf.roc)

# 3) Plot individual ROC curve
pdf(file="ROC_rf.pdf")
plot.roc(rf.roc,xlim=c(1,0),main="ROC curve of Random Forest",col="blue")
legend("topleft",legend=c("AUC=0.7367"))
dev.off()

# Important features for random forest
important_rf <- as.data.frame(rf_fit$importance[,c(2,3,4)])
colnames(important_rf) <- c("MeanDecreaseAccuracyClass","MeanDecreaseAccuracy","MeanDecreaseGini")
important_rf <- important_rf[order(important_rf$MeanDecreaseAccuracyClass,decreasing = TRUE),]
write.csv(important_rf,'ImportantFeaturesByRF.csv')

#----------------------------------------------------------------------
#save variables in R data file
#save(x,y,testx,testy,file='Variables.Rdata')

#----------------------------------------------------------------------
# Plot combined ROC curve
# Need to load svm roc results

pdf(file="mergedROC2.pdf")
plot.roc(cv.lasso.roc,xlim=c(1,0),main="ROC curves comparison",col=c(1))
lines(nb.roc$sensitivities,nb.roc$specificities,col=c(4))
lines(rf.roc$sensitivities,rf.roc$specificities,col=c(6))
lines(svml.roc$sensitivities,svml.roc$specificities,col=c(2))
#lines(svmk.roc$sensitivities,svmk.roc$specificities,col=c(3))
legend("bottomright",col=c(1,4,6,2), lty=1, legend=c("LogisticRegression-L1","Naive Bayes","Random Forest","SVM-linear"))
dev.off()


#nb.roc
#cv.lasso.roc
#rf.auc
#svml.roc


