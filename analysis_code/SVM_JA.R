#SVM
setwd("/Users/ja/Documents/Computational Methods/Project/Model")
#load("~/Documents/Computational Methods/Project/Model/GPCR_IC.Rdata")

#------------------------------------------------------------------------------------
# SVM - linear 

# Tuning parameters
tc <- tune.control(cross=10)
cost_array <- c(seq(0.001,0.05,0.001))
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

#------------------------------------------------------------------------------------
# SVM - kernel 

# Tuning parameters
g_value = seq(0.01,0.035,0.005)
c_value = seq(2,12,1)
tune.out_r <- tune(svm,train.x = x,train.y = y, probability = TRUE, kernel="radial",ranges=list(gamma=g_value,cost=c_value))
hm <- tune.out_r$performance[,1:3]

# Ploting Kernel SVM misclassifiation
pdf('GaussianSVM_misclass.pdf')
p <- ggplot(hm,aes(cost,gamma)) + geom_tile(aes(fill=error),color="white") + scale_fill_gradient(low="white",high="steelblue")
p + ggtitle("SVM Misclassification based on kernel bandwidth and soft-margin parameter") + xlab("soft margin parameter(cost)") + ylab("kernel bandwidith (gamma)")
dev.off()

# Select best performing model for prediction
rperform <- tune.out_r$performances
bestk <- tune.out_r$best.model
pred_y_k <- predict(bestk,testx,probability = TRUE)
prob_k <- attr(pred_y_k,"probabilities")[,1]
kmisclass <- mean(pred_y_k!=testy)
svmk.roc <- roc(testy,prob_k)

#auc(svml.roc)
#auc(svmk.roc)
