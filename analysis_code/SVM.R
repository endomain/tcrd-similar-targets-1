#SVM
setwd("/Users/ja/Documents/Computational Methods/Project/Model")
# Use the same x, y, trainx, trainy, testx, testy for the "combined_analysis_JA.R" file

#------------------------------------------------------------------------------------
# SVM - linear 

# Tuning parameters
tc <- tune.control(cross=10)
cost_array <- c(seq(0.001,0.01,0.0005))
tune.out <- tune(svm,train.x = x,train.y = y, probability=TRUE, kernel="linear",ranges=list(cost=cost_array),tunecontrol = tc)
lperform2 <- tune.out$performance

# Ploting linear SVM misclassifiation
pdf('LinearSVM_misclass.pdf')
plot(y=lperform2$error,x=lperform2$cost,xlab='soft-margin parameter',ylab='misclassification rate',main = 'Misclassification rate for linear SVM based on soft-margin parameter')
dev.off()

# Select best performing model for prediction
bestl <- tune.out$best.model
pred_y_l <- predict(bestl,testx,probability = TRUE) 
prob <- attr(pred_y_l,"probabilities")[,1]
lmisclass <- mean(pred_y_l!=testy)
svml.roc <- roc(testy,prob)

