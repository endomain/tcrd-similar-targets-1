library(RMySQL)
library(dplyr)

mydb = dbConnect(MySQL(), user='root', password='jhoon11', dbname='pharos', host='localhost')

kinase_ids<-dbSendQuery(mydb,'SELECT id from target where idgfam="Kinase";')
kinase_ids<-fetch(kinase_ids, n=-1)
kinase_ids<-kinase_ids$id

result<-data.frame()

dbSendQuery(mydb,
            'CREATE TEMPORARY TABLE kinases (id int) as (SELECT id from target where idgfam="Kinase");')

tdl_labels_query<-dbSendQuery(mydb,
                            'SELECT tdl from target where id IN (select id from kinases);')

tdl_labels = fetch(tdl_labels_query, n=-1)

query_all<-function(kinase_ids){
  for(id in kinase_ids){
    query<-dbSendQuery(mydb,sprintf(
                'SELECT 
    	CONCAT(
                "chembl_activity",
                "_",
                cmpd_chemblid,
                "_",
                act_type
    ) AS results
                FROM chembl_activity as results
                WHERE target_id=%s
                
                UNION
                
                SELECT DISTINCT(CONCAT("drug_activity","_",REPLACE(go_term," ","_"))) as results
                FROM compartment 
                WHERE protein_id=%s
                
                UNION
                
                SELECT 
                CONCAT(
                "expression_",
                Replace(tissue," ","_"),
                "_",
                Replace(qual_value," ","_")
                ) AS expression
                FROM expression 
                WHERE protein_id=%s 
                # AND etype="HPM Protein"
                GROUP BY tissue	
                
                UNION
                
                SELECT (CONCAT("generif_pubmed","_",pubmed_ids))
                FROM generif 
                WHERE protein_id=%s
                
                UNION
                
                
                SELECT REPLACE(CONCAT("gene_attribute","_",REPLACE(gene_attribute.name,CONCAT("/",gene_attribute_type.name),""))," ","_") as results
                FROM gene_attribute
                JOIN gene_attribute_type ON gene_attribute_type.id=gene_attribute.gat_id 
                WHERE 
                gat_id=103 AND 
                protein_id=%s
                
                UNION
                
                SELECT CONCAT("panther_class_id_",pcid) AS results 
                FROM panther_class 
                WHERE ID IN (SELECT panther_class_id from p2pc WHERE protein_id="%s")
                
                UNION
                
                SELECT (CONCAT("pathway","_",pwtype,"_",id_in_source)) as results
                FROM pathway 
                WHERE protein_id=%s AND pwtype="Reactome"
                
                UNION
                
                SELECT (CONCAT("phenotype","_",REPLACE(term_id,":","_"))) as results
                FROM phenotype 
                WHERE protein_id=%s AND ptype="JAX/MGI Human Ortholog Phenotype"
                
                UNION
                
                SELECT CONCAT("ppi","_",(
                CASE 
                WHEN protein1_id = %s THEN protein2_id 
                ELSE protein1_id END)) AS ppi_protein_id
                FROM ppi WHERE protein1_id=%s OR protein2_id=%s
                
                UNION
                
                
                SELECT CONCAT("target2disease","_",REPLACE(name," ","_")) as target2disease_id
                FROM target2disease
                WHERE target_id=%s;
                ',id,id,id,id,id,id,id,id,id,id,id,id))
    
    data = fetch(query, n=-1)
    all_variables<-append(data$result,names(result))
    print(length(all_variables))
    print(id)
    row<-data.frame(as.list(setNames(c(ifelse(all_variables %in% data$results,1,0)),
                                     all_variables)))
    result<-bind_rows(result,row)
    
  }
  return(result)
}

final_result1_30<-query_all(kinase_ids[1:30])
final_result31_60<-query_all(kinase_ids[31:60])
final_result61_90<-query_all(kinase_ids[61:90])
final_result91_120<-query_all(kinase_ids[91:120])
final_result121_150<-query_all(kinase_ids[121:150])
final_result151_180<-query_all(kinase_ids[151:180])
final_result181_210<-query_all(kinase_ids[181:210])
final_result211_240<-query_all(kinase_ids[211:240])
final_result241_260<-query_all(kinase_ids[241:260])
final_result261_280<-query_all(kinase_ids[261:280])
final_result281_310<-query_all(kinase_ids[281:310])
final_result311_340<-query_all(kinase_ids[311:340])
final_result341_370<-query_all(kinase_ids[341:370])
final_result381_410<-query_all(kinase_ids[381:410])
final_result401_430<-query_all(kinase_ids[401:430])
final_result431_460<-query_all(kinase_ids[431:460])
final_result461_490<-query_all(kinase_ids[461:490])
final_result491_520<-query_all(kinase_ids[491:520])
final_result521_540<-query_all(kinase_ids[521:540])
final_result541_578<-query_all(kinase_ids[541:578])
cum_results<-bind_rows(final_result1_30,final_result31_60)
cum_results<-bind_rows(cum_results,final_result61_90)
cum_results<-bind_rows(cum_results,final_result91_120)
cum_results<-bind_rows(cum_results,final_result121_150)
cum_results<-bind_rows(cum_results,final_result151_180)
cum_results<-bind_rows(cum_results,final_result181_210)
cum_results<-bind_rows(cum_results,final_result211_250)
cum_results<-bind_rows(cum_results,final_result251_300)
cum_results<-bind_rows(cum_results,final_result301_350)
cum_results<-bind_rows(cum_results,final_result351_400)
cum_results<-bind_rows(cum_results,final_result401_430)
cum_results<-bind_rows(cum_results,final_result431_460)
cum_results<-bind_rows(cum_results,final_result461_490)
cum_results<-bind_rows(cum_results,final_result491_520)
cum_results<-bind_rows(cum_results,final_result521_540)
cum_results<-bind_rows(cum_results,final_result541_578)

# daisy_dist<-daisy(clean_data,metric="euclidean")
# kinase_hclust<-hclust(daisy_dist,method="complete")
library(flexclust)
library(dendextend)
library(RColorBrewer)
library(ggplot2)

kinase_10_clusters_dend<-as.dendrogram(kinase_hclust)
kinase_10_clusters = cutree(kinase_hclust, k=12, order_clusters_as_data = FALSE)

# So that we can color the branches of 
palette <- brewer.pal(n = 12, name = "Paired")

kinase_10_clusters_dend %>%
  set("branches_k_color",
      value=palette[kinase_10_clusters])%>%
  set("labels_col", 
      palette[kinase_10_clusters]) %>%
  plot(main = "Dendrogram for Hierarchical Clustering for bcwisc\nwith true labels")

palette <- brewer.pal(n = 2, name = "PuOr")

colored_bars(colors = palette[factor(tdl_labels$tdl!="Tclin")],
             dend = kinase_10_clusters_dend,
             rowLabels = results$target_id[which(tdl_labels$tdl=="Tclin")])

results$target_id[which(tdl_labels$tdl=="Tclin")]
###########################################################
# Naive Bayes classifier filtered on the clusters
# 
# for TDL class
###########################################################
kinase_10_cluster_labels = cutree(kinase_hclust, k=10, order_clusters_as_data = TRUE)
results_m<-as.matrix(cum_results)
results<-cum_results
# Append 10 clusters
results$kdl_class<-tdl_labels$tdl
results$kinase_10_cluster_labels<-kinase_10_cluster_labels
  
library(e1071)

# For cluster 5, run naive based of TDL class
cluster_5<-subset(results,kinase_10_cluster_labels==5)
cluster_5_target_ids<-cluster_5$target_id
cluster_5_kinase_10_cluster_labels<-cluster_5$kinase_10_cluster_labels
cluster_5_kdl_class<-cluster_5$kdl_class

which(colnames(cluster_5)=="kinase_10_cluster_labels")
which(colnames(cluster_5)=="kdl_class")
which(colnames(cluster_5)=="target_id")

cluster_5<-cluster_5[,c(1:84581,84583)]
cluster_5<-as.matrix(cluster_5)

model <- naiveBayes(cluster_5_kdl_class,cluster_5)
pred <- predict(model,cluster_5)
table(pred, cluster_5$kdl_class)

# For cluster 8, run naive based of TDL class
cluster_8<-subset(results,kinase_10_cluster_labels==8)
cluster_8_target_ids<-cluster_8$target_id
cluster_8_kinase_10_cluster_labels<-cluster_8$kinase_10_cluster_labels
cluster_8_kdl_class<-cluster_8$kdl_class

which(colnames(cluster_8)=="kinase_10_cluster_labels")
which(colnames(cluster_8)=="kdl_class")
which(colnames(cluster_8)=="target_id")

cluster_8<-cluster_8[,c(1:84581,84583)]
# cluster_8<-as.matrix(cluster_8)



model <- naiveBayes(kdl_class~.,data=cluster_8)
pred <- predict(model,cluster_8)

table(pred, cluster_8_kdl_class)

