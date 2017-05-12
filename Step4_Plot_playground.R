library(ComplexHeatmap)
library(colorspace)
library(circlize)
library(heatmap3)

######################################################
data_viz_df<-all_tcrd
idl_idgfam<-read.csv("all_tcrd_tdl_idgfam.csv")

colnames(idl_idgfam)[1]<-"target_id"
data_viz_df<-cbind(data_viz_df,data.frame(idl_idgfam[match(all_tcrd$target_id, idl_idgfam$target_id),])$idgfam)
colnames(data_viz_df)[161039]<-'idgfam'
data_viz_df$predict<-predict(model,data_viz_df)

######################################################

table(data_viz_df$predict,data_viz_df$tdl)
table(data_viz_df$predict,data_viz_df$idgfam)

# cols_to_subtract<-c(which(colnames(data_viz_df)=="target_id"),
#                     which(colnames(data_viz_df)=="tdl"),
#                     which(colnames(data_viz_df)=="idgfam"),
#                     which(colnames(data_viz_df)=="predict"))
# 
# data_viz_df<-data_viz_df[,-cols_to_subtract]
data_viz_df<-data_viz_df[,get_colSums(data_viz_df,20)]
heatmap_mat<-data.matrix(data_viz_df)


idgfam_annotation<-data.frame()
idgfam_annotation<-data.frame(idgfam=data_viz_df[data_viz_df$target_id == as.numeric(rownames(heatmap_mat))]$idgfam,
                              tdl=data_viz_df[data_viz_df$target_id == as.numeric(rownames(heatmap_mat))]$tdl,
                              predict=data_viz_df[data_viz_df$target_id == as.numeric(rownames(heatmap_mat))]$predict)

pdf("heatmap.pdf",width=40,height=70)

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
# #                                        "Kinase"="green")))
ha2 = rowAnnotation(df = idgfam_annotation,
                    col = list(predict = c("1"="blue",
                                           "-1"="red")))
# 
# ha2 = rowAnnotation(df = idgfam_annotation, 
#                    col = list(tdl = c("Tclin" = "green", 
#                                          "Tbio" = "blue", 
#                                          "Tchem"="red",
#                                       "Tdark"="orange")))


h1

dev.off()
######################################################
# Features Correlation matrix
# as a network cluster
######################################################
library(qgraph)
library(parcor)
cor<-cor(heatmap_mat)

pdf("correlation_network.pdf",width=40,height=50)
qgraph(cor,graph="cor",
       layout="spring",
       minimum=0.5,
       label.cex=5)
dev.off()