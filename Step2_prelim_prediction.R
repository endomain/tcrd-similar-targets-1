library(caTools)
library(e1071)
library(randomForest)
library(caret)
######################################################
# Remove Sparse features
# Use temporary data frame to work with
######################################################
temp<-all_tcrd_reduced[,colSums(all_tcrd_reduced) > 40]
######################################################
# Merge our 
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
which(colnames(all_tcrd_reduced_train)=="tdl")
which(colnames(all_tcrd_reduced_train)=="idgfam")

model<-randomForest(y~.-tdl-idgfam,data=all_tcrd_reduced_train)
pred <- predict(model,all_tcrd_reduced_test[,c(-which(colnames(all_tcrd_reduced_test)=="y"))])
table(pred, all_tcrd_reduced_test$y)

View(varImp(model))


######################################################
bah<-all_tcrd
idl_idgfam<-read.csv("all_tcrd_tdl_idgfam.csv")

colnames(idl_idgfam)[1]<-"target_id"
# bah<-cbind(bah,data.frame(idl_idgfam[match(all_tcrd$target_id, idl_idgfam$target_id),])$tdl)
# colnames(bah)[88598]<-'tdl'
bah<-cbind(bah,data.frame(idl_idgfam[match(all_tcrd$target_id, idl_idgfam$target_id),])$idgfam)
colnames(bah)[88598]<-'idgfam'
bah$predict<-predict(model,bah)

table(bah$predict,bah$tdl)
table(bah$predict,bah$idgfam)
# table(bah$predict,bah$chembl_activity_CHEMBL42_Ki)
# table(bah$predict,bah$drug_activity_Cytosol)
# table(bah$tdl,bah$drug_activity_Cytosol)
# table(bah$predict,bah$tdl)

library(heatmap3)
rownames(bah)<-bah[,which(colnames(bah)=="target_id")]

cols_to_subtract<-c(which(colnames(bah)=="target_id"),
                    which(colnames(bah)=="tdl"),
                    which(colnames(bah)=="idgfam"),
                    which(colnames(bah)=="predict"))

heatmap_mat<-data.matrix(bah[,-cols_to_subtract])
heatmap_mat<-heatmap_mat[,colSums(heatmap_mat) > 40]

library(ComplexHeatmap)
library(colorspace)
dev.off()
idgfam_annotation<-data.frame()
idgfam_annotation<-data.frame(idgfam=bah[bah$target_id == as.numeric(rownames(heatmap_mat))]$idgfam,
                              tdl=bah[bah$target_id == as.numeric(rownames(heatmap_mat))]$tdl,
                              predict=bah[bah$target_id == as.numeric(rownames(heatmap_mat))]$predict)

pdf("heatmap.pdf",width=40,height=50)

# ha1 = HeatmapAnnotation(df = df,
#                         col = list(type = c("a" = "red", "b" = "blue"),
#                                    age = colorRamp2(c(0, 20), c("white", "red"))))
# ha2 = HeatmapAnnotation(df = data.frame(age = sample(1:20, 10)),
#                         col = list(age = colorRamp2(c(0, 2), c("IC", "GPCR","Kinase"))))



h1<-Heatmap(heatmap_mat, km=12,
        name="h1",
        # top_annotation = ha1, 
        # bottom_annotation = ha2,
        col = colorRamp2(c(0, 1.2), c("white", "blue")),
        use_raster=FALSE)

# ha1 = rowAnnotation(df = idgfam_annotation$idgfam, 
#                    col = list(idgfam = c("IC" = "red", 
#                                        "GPCR" = "blue", 
#                                        "Kinase"="green")))
ha2 = rowAnnotation(df = idgfam_annotation, 
                    col = list(predict = c("1"="blue",
                                           "0"="white")))
# 
# ha2 = rowAnnotation(df = idgfam_annotation, 
#                    col = list(tdl = c("Tclin" = "green", 
#                                          "Tbio" = "blue", 
#                                          "Tchem"="red",
#                                       "Tdark"="orange")))


h1+ha2

dev.off()

library(qgraph)
library(parcor)
cor<-cor(heatmap_mat)
qgraph(cor,graph="cor",layout="spring",minimum=0.2,label.cex=5)
