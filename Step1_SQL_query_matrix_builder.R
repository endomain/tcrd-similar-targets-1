library(RMySQL)
library(dplyr)

mydb = dbConnect(MySQL(), user='root', password='jhoon11', dbname='pharos', host='localhost')

# Get all kinase ids
dbSendQuery(mydb,
            'CREATE TEMPORARY TABLE kinases (id int) as (SELECT id from target where idgfam="Kinase");')
kinase_ids<-dbSendQuery(mydb,'SELECT id from target where idgfam="Kinase";')
kinase_ids<-fetch(kinase_ids, n=-1)
kinase_ids<-kinase_ids$id
kinase_tdl_labels_query<-dbSendQuery(mydb,
                                     'SELECT tdl from target where id IN (select id from kinases);')

kinase_tdl_labels = fetch(kinase_tdl_labels_query, n=-1)

# Get all IC ids and labels

dbSendQuery(mydb,
            'CREATE TEMPORARY TABLE ion_channel (id int) as (SELECT id from target where idgfam="IC");')

IC_ids<-dbSendQuery(mydb,'SELECT id from target where idgfam="IC";')
IC_ids<-fetch(IC_ids, n=-1)
IC_ids<-IC_ids$id
IC_tdl_labels_query<-dbSendQuery(mydb,
                                 'SELECT tdl from target where id IN (select id from ion_channel);')
IC_tdl_labels = fetch(IC_tdl_labels_query, n=-1)

# Get all GPCR ids
dbSendQuery(mydb,
            'CREATE TEMPORARY TABLE gcpr (id int) as (SELECT id from target where idgfam="GPCR");')
GPCR_ids<-dbSendQuery(mydb,'SELECT id from target where idgfam="GPCR";')
GPCR_ids<-fetch(GPCR_ids, n=-1)
GPCR_ids<-GPCR_ids$id
GPCR_tdl_labels_query<-dbSendQuery(mydb,
                                   'SELECT tdl from target where id IN (select id from gcpr);')

GPCR_tdl_labels = fetch(GPCR_tdl_labels_query, n=-1)


result<-data.frame()
## This function queries our database
query_all<-function(kinase_ids){
  for(id in kinase_ids){
    query<-dbSendQuery(mydb,sprintf(
      'SELECT CONCAT("chembl_activity","_",cmpd_chemblid,"_",act_type) AS results
      FROM chembl_activity as results
      WHERE target_id=%s
      
      UNION
      
      SELECT DISTINCT(CONCAT("drug_activity","_",REPLACE(go_term," ","_"))) AS results
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
      AND etype="Consensus"
      GROUP BY tissue	
      
      UNION
      
      #SELECT (CONCAT("generif_pubmed","_",pubmed_ids))
      #FROM generif 
      #WHERE protein_id=%s
      
      #UNION
      
      #SELECT REPLACE(CONCAT("gene_attribute","_",REPLACE(gene_attribute.name,CONCAT("/",gene_attribute_type.name),""))," ","_") as results
      #FROM gene_attribute
      #JOIN gene_attribute_type ON gene_attribute_type.id=gene_attribute.gat_id 
      #WHERE 
      #gat_id=103 AND 
      #protein_id=%s
      
      #UNION
      
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
    all_variables<-make.names(unique(append(make.names(data$result),names(result))))
    # print(all_variables)
    print(length(all_variables))
    print(id)
    row<-data.frame(as.list(setNames(c(ifelse(all_variables %in% data$result,1,0)),
                                     all_variables)))
    # print(row)
    result<-bind_rows(result,row)
    
  }
  return(result)
}
debug<-query_all(kinase_ids[1])

debug2<-query_all(kinase_ids[2])

names(debug)

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
final_result341_380<-query_all(kinase_ids[341:380])
final_result381_410<-query_all(kinase_ids[381:410])
final_result411_430<-query_all(kinase_ids[411:430])
final_result431_460<-query_all(kinase_ids[431:460])
final_result461_490<-query_all(kinase_ids[461:490])
final_result491_520<-query_all(kinase_ids[491:520])
final_result521_540<-query_all(kinase_ids[521:540])
final_result541_578<-query_all(kinase_ids[541:578])
kinase_results<-bind_rows(final_result1_30,final_result31_60)
kinase_results<-bind_rows(kinase_results,final_result61_90)
kinase_results<-bind_rows(kinase_results,final_result91_120)
kinase_results<-bind_rows(kinase_results,final_result121_150)
kinase_results<-bind_rows(kinase_results,final_result151_180)
kinase_results<-bind_rows(kinase_results,final_result181_210)
kinase_results<-bind_rows(kinase_results,final_result211_240)
kinase_results<-bind_rows(kinase_results,final_result241_260)
kinase_results<-bind_rows(kinase_results,final_result261_280)
kinase_results<-bind_rows(kinase_results,final_result281_310)
kinase_results<-bind_rows(kinase_results,final_result311_340)
kinase_results<-bind_rows(kinase_results,final_result341_380)
kinase_results<-bind_rows(kinase_results,final_result381_410)
kinase_results<-bind_rows(kinase_results,final_result411_430)
kinase_results<-bind_rows(kinase_results,final_result431_460)
kinase_results<-bind_rows(kinase_results,final_result461_490)
kinase_results<-bind_rows(kinase_results,final_result491_520)
kinase_results<-bind_rows(kinase_results,final_result521_540)
kinase_results<-bind_rows(kinase_results,final_result541_578)


gcpr_1_30<-query_all(GPCR_ids[1:30])
gcpr_31_60<-query_all(GPCR_ids[31:60])
gcpr_61_90<-query_all(GPCR_ids[61:90])
gcpr_91_120<-query_all(GPCR_ids[91:120])
gcpr_121_150<-query_all(GPCR_ids[121:150])
gcpr_151_180<-query_all(GPCR_ids[151:180])
gcpr_181_210<-query_all(GPCR_ids[181:210])
gcpr_211_240<-query_all(GPCR_ids[211:240])
gcpr_241_270<-query_all(GPCR_ids[241:270])
gcpr_271_300<-query_all(GPCR_ids[271:300])
gcpr_301_330<-query_all(GPCR_ids[301:330])
gcpr_331_360<-query_all(GPCR_ids[331:360])
gcpr_361_390<-query_all(GPCR_ids[361:390])
gcpr_391_406<-query_all(GPCR_ids[391:406])

gcpr_results<-bind_rows(gcpr_31_60,gcpr_61_90)
gcpr_results<-bind_rows(gcpr_results,gcpr_121_150)
gcpr_results<-bind_rows(gcpr_results,gcpr_151_180)
gcpr_results<-bind_rows(gcpr_results,gcpr_181_210)
gcpr_results<-bind_rows(gcpr_results,gcpr_211_240)
gcpr_results<-bind_rows(gcpr_results,gcpr_241_270)
gcpr_results<-bind_rows(gcpr_results,gcpr_301_330)
gcpr_results<-bind_rows(gcpr_results,gcpr_331_360)
gcpr_results<-bind_rows(gcpr_results,gcpr_361_390)
gcpr_results<-bind_rows(gcpr_results,gcpr_391_406)
# Order of adding chunks with large number of features make this go a lot faster
gcpr_results<-bind_rows(gcpr_results,gcpr_91_120)
gcpr_results<-bind_rows(gcpr_results,gcpr_271_300)
gcpr_results<-bind_rows(gcpr_results,gcpr_1_30)



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

######################################################
# Add appropriate Labels for each class
######################################################

kinase_results$target_id<-unlist(kinase_ids)
kinase_results$tdl<-unlist(kinase_tdl_labels)

gcpr_results$target_id<-unlist(GPCR_ids)
gcpr_results$tdl<-unlist(GPCR_tdl_labels)

IC_results$target_id<-unlist(IC_ids)
IC_results$tdl<-unlist(IC_tdl_labels)

######################################################
# Combine all three into one giant data frame
######################################################
Rprof ( tf <- "log.log",  memory.profiling = TRUE )
all_tcrd<-bind_rows(kinase_results,gcpr_results)
all_tcrd<-bind_rows(all_tcrd,IC_results)
Rprof ( NULL ) ; print ( summaryRprof ( tf )  )
######################################################
# Zero all NA values. 
# Because we are checking for 'presence' of an association
# and not obtaining any quantitative values, we felt that
# 0ing the features not listed was a safe assumption.
######################################################

all_tcrd[is.na(all_tcrd)]<-0

######################################################
# Combine all the gold standards, extracted from
# WITHDARWN database.
######################################################
gold_GPCR<-read.csv('toxic_targets/Gold_GPCR.csv')
gold_Kinase<-read.csv('toxic_targets/gold_Kinase.csv')
gold_IC<-read.csv('toxic_targets/gold_IC.csv')
all_gold<-bind_rows(gold_IC,gold_GPCR)
# Kinase csv has an extra row
all_gold<-bind_rows(all_gold,gold_Kinase[,c(-1)])


######################################################
# 
######################################################
all_tcrd_reduced<-all_tcrd[,-c(which(colnames(all_tcrd)=="tdl"))]
# Select for training
all_tcrd_reduced<-all_tcrd_reduced %>% filter(target_id %in% all_gold$id)