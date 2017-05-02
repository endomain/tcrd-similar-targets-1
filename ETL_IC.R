library(RMySQL)
library(dplyr)

mydb = dbConnect(MySQL(), user='root', password='jhoon11', dbname='pharos', host='localhost')

IC_ids<-dbSendQuery(mydb,'SELECT id from target where idgfam="IC";')
IC_ids<-fetch(IC_ids, n=-1)
IC_ids<-IC_ids$id

result<-data.frame()

dbSendQuery(mydb,
            'CREATE TEMPORARY TABLE ion_channel (id int) as (SELECT id from target where idgfam="IC");')

IC_tdl_labels_query<-dbSendQuery(mydb,
                                   'SELECT tdl from target where id IN (select id from ion_channel);')

IC_tdl_labels = fetch(IC_tdl_labels_query, n=-1)

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

IC_1_30<-query_all(IC_ids[1:30])
IC_31_60<-query_all(IC_ids[31:60])
IC_61_90<-query_all(IC_ids[61:90])
IC_91_120<-query_all(IC_ids[91:120])
IC_121_150<-query_all(IC_ids[121:150])
IC_151_180<-query_all(IC_ids[151:180])
IC_181_210<-query_all(IC_ids[181:210])
IC_211_240<-query_all(IC_ids[211:240])
IC_241_270<-query_all(IC_ids[241:270])
IC_271_300<-query_all(IC_ids[271:300])
IC_301_330<-query_all(IC_ids[301:330])
IC_331_342<-query_all(IC_ids[331:342])

IC_results<-bind_rows(IC_1_30,IC_31_60)
IC_results<-bind_rows(IC_results,IC_61_90)
IC_results<-bind_rows(IC_results,IC_91_120)
IC_results<-bind_rows(IC_results,IC_121_150)
IC_results<-bind_rows(IC_results,IC_151_180)
IC_results<-bind_rows(IC_results,IC_181_210)
IC_results<-bind_rows(IC_results,IC_211_240)
IC_results<-bind_rows(IC_results,IC_241_270)
IC_results<-bind_rows(IC_results,IC_271_300)
IC_results<-bind_rows(IC_results,IC_301_330)
IC_results<-bind_rows(IC_results,IC_331_342)

IC_results$target_id<-unlist(IC_ids)
IC_results$tdl_class<-unlist(IC_tdl_labels)
IC_results[is.na(IC_results)]<-0
write.csv(IC_results,"IC_results.csv")
# Figure out which column numbers are 'tdf_class' and 'target_id'
which(colnames(IC_results)=="tdl_class")
which(colnames(IC_results)=="target_id")
## Reduce number of features
IC_results_reduced<-IC_results[,c(-37033,-37034)]
which(colnames(IC_results_reduced)=="tdl_class")
which(colnames(IC_results_reduced)=="target_id")
IC_results_reduced<-IC_results_reduced[,colSums(IC_results_reduced) > 50]

library(e1071)

## Trying to predict tdl_class an a proof of concept
IC_results_reduced$tdl_class<-unlist(IC_tdl_labels)
IC_results_reduced$tdl_class<-as.factor(IC_results_reduced$tdl_class)
model <- naiveBayes(tdl_class~.,data=IC_results_reduced)
which(colnames(IC_results_reduced)=="tdl_class")
pred <- predict(model,IC_results_reduced[,c(-232)])
table(pred, IC_results_reduced$tdl_class)

### Now for real
gold_IC<-read.csv('toxic_targets/gold_IC.csv')
IC_results_reduced<-IC_results[,c(-37033)]
IC_results_reduced<-IC_results_reduced %>% filter(target_id %in% gold_IC$id)
temp<-IC_results_reduced[,colSums(IC_results_reduced) > 10]
temp<-merge(temp,gold_IC[,c('id','y')],by.x='target_id',by.y='id')
temp<-temp[,c(-1)]
library(caTools)
spl<-sample.split(temp$y,SplitRatio = 0.8)
IC_results_reduced_train<-temp[spl,]
IC_results_reduced_test<-temp[!spl,]

# Make sure we got rid of class

IC_results_reduced_train$y<-as.factor(IC_results_reduced_train$y)
model <- naiveBayes(y~.,data=IC_results_reduced_train)
# model<-glmnet(IC_results_reduced_train,y)
pred <- predict(model,IC_results_reduced_test[,c(-which(colnames(IC_results_reduced_test)=="y"))])
table(pred, IC_results_reduced_test$y)

