library(RMySQL)
library(dplyr)

mydb = dbConnect(MySQL(), user='root', password='SET_PASSWORD', dbname='pharos', host='localhost')

all_gold<-read.csv('toxic_targets/GoldStandards.csv')


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

# Get all GPCR ids and labels
dbSendQuery(mydb,
            'CREATE TEMPORARY TABLE gcpr (id int) as (SELECT id from target where idgfam="GPCR");')
GPCR_ids<-dbSendQuery(mydb,'SELECT id from target where idgfam="GPCR";')
GPCR_ids<-fetch(GPCR_ids, n=-1)
GPCR_ids<-GPCR_ids$id
GPCR_tdl_labels_query<-dbSendQuery(mydb,
                                   'SELECT tdl from target where id IN (select id from gcpr);')

GPCR_tdl_labels = fetch(GPCR_tdl_labels_query, n=-1)

# Get all idfgam=NULL ids and TDL labels
dbSendQuery(mydb,
            'CREATE TEMPORARY TABLE null_idg (id int) as (SELECT id from target where idgfam IS NULL);')
null_ids<-dbSendQuery(mydb,'SELECT id from target where idgfam IS NULL;')
null_ids<-fetch(null_ids, n=-1)
null_ids<-null_ids$id
# Get only gold standard IDS
null_ids<-null_ids[null_ids %in% all_gold$id]
null_tdl_labels_query<-dbSendQuery(mydb,
                                   'SELECT id,tdl from target where id IN (select id from null_idg);')

null_tdl_labels = fetch(null_tdl_labels_query, n=-1)
null_tdl_labels<-null_tdl_labels %>% filter(id %in% null_ids)

# Get all idfgam=NR and TDL labels
dbSendQuery(mydb,
            'CREATE TEMPORARY TABLE nr (id int) as (SELECT id from target where idgfam ="NR");')
nr_ids<-dbSendQuery(mydb,'SELECT id from target where idgfam ="NR";')
nr_ids<-fetch(nr_ids, n=-1)
nr_ids<-nr_ids$id
# Get only gold standard IDS
nr_ids<-nr_ids[nr_ids %in% all_gold$id]
nr_tdl_labels_query<-dbSendQuery(mydb,
                                   'SELECT id,tdl from target where id IN (select id from nr);')

nr_tdl_labels = fetch(nr_tdl_labels_query, n=-1)
nr_tdl_labels<-nr_tdl_labels %>% filter(id %in% nr_ids)


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
    # Each row is added after transforming the combined feature list of existing
    # data frame and checks if new features are listed in the new dataset.
    row<-data.frame(as.list(setNames(c(ifelse(all_variables %in% data$result,1,0)),
                                     all_variables)))
    # print(row)
    result<-bind_rows(result,row)
    
  }
  return(result)
}


########################################################################
# Ugly coding, but had to to this to get around R memory limits and speed 
# issues related to querying and expanding each feature.
########################################################################
# Query and Bind KINASES
########################################################################

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

########################################################################
# Query and Bind GPC receptors
########################################################################
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

########################################################################
# Query and row bind Ion channels
########################################################################
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

########################################################################
# Query and row bind IDG family = NULL
########################################################################
null_1_30<-query_all(null_ids[1:30])
null_31_60<-query_all(null_ids[31:60])
null_61_90<-query_all(null_ids[61:90])
null_91_120<-query_all(null_ids[91:120])
null_121_150<-query_all(null_ids[121:150])
null_151_180<-query_all(null_ids[151:180])
null_181_210<-query_all(null_ids[181:210])
null_211_240<-query_all(null_ids[211:240])
null_241_270<-query_all(null_ids[241:270])
null_271_300<-query_all(null_ids[271:300])
null_301_330<-query_all(null_ids[301:330])
null_331_360<-query_all(null_ids[331:360])
null_361_390<-query_all(null_ids[361:390])
null_391_420<-query_all(null_ids[391:420])
null_421_450<-query_all(null_ids[421:450])
null_451_480<-query_all(null_ids[451:480])
null_481_510<-query_all(null_ids[481:510])
null_511_540<-query_all(null_ids[511:540])
null_541_570<-query_all(null_ids[541:570])
null_571_600<-query_all(null_ids[571:600])
null_601_630<-query_all(null_ids[601:630])
null_631_660<-query_all(null_ids[631:660])
null_661_690<-query_all(null_ids[661:690])
null_691_720<-query_all(null_ids[691:720])
null_721_750<-query_all(null_ids[721:750])
null_751_763<-query_all(null_ids[751:763])

null_results<-bind_rows(null_1_30,null_31_60)
null_results<-bind_rows(null_results,null_61_90)
null_results<-bind_rows(null_results,null_91_120)
null_results<-bind_rows(null_results,null_121_150)
null_results<-bind_rows(null_results,null_151_180)
null_results<-bind_rows(null_results,null_181_210)
null_results<-bind_rows(null_results,null_211_240)
null_results<-bind_rows(null_results,null_241_270)
null_results<-bind_rows(null_results,null_271_300)
null_results<-bind_rows(null_results,null_301_330)
null_results<-bind_rows(null_results,null_331_360)
null_results<-bind_rows(null_results,null_361_390)
null_results<-bind_rows(null_results,null_391_420)
null_results<-bind_rows(null_results,null_421_450)
null_results<-bind_rows(null_results,null_451_480)
null_results<-bind_rows(null_results,null_481_510)
null_results<-bind_rows(null_results,null_511_540)
null_results<-bind_rows(null_results,null_541_570)
null_results<-bind_rows(null_results,null_571_600)
null_results<-bind_rows(null_results,null_601_630)
null_results<-bind_rows(null_results,null_631_660)
null_results<-bind_rows(null_results,null_661_690)
null_results<-bind_rows(null_results,null_691_720)
null_results<-bind_rows(null_results,null_721_750)
null_results<-bind_rows(null_results,null_751_763)

########################################################################
# Query and row bind IDG family = NR
########################################################################

nr_1_28<-query_all(nr_ids[1:28])
nr_results<-nr_1_28

########################################################################
# Add appropriate target_ids and labels.
# *** In the right order of the queries - especially for GPCR
########################################################################

null_results$target_id<-unlist(null_ids)
null_results$tdl<-unlist(null_tdl_labels$tdl)

nr_results$target_id<-unlist(nr_ids)
nr_results$tdl<-unlist(nr_tdl_labels$tdl)

kinase_results$target_id<-unlist(kinase_ids)
kinase_results$tdl<-unlist(kinase_tdl_labels)

########################################################################
# caveat for above optimization for bind_rows for GPCR is that
# is that it messed up the original target_id label order.
########################################################################
gcpr_results$target_id<-unlist(c(GPCR_ids[31:60],
                                 GPCR_ids[61:90],
                                 GPCR_ids[121:150],
                                 GPCR_ids[151:180],
                                 GPCR_ids[181:210],
                                 GPCR_ids[211:240],
                                 GPCR_ids[241:270],
                                 GPCR_ids[301:330],
                                 GPCR_ids[331:360],
                                 GPCR_ids[361:390],
                                 GPCR_ids[391:406],
                                 GPCR_ids[91:120],
                                 GPCR_ids[271:300],
                                 GPCR_ids[1:30]
))
gcpr_results$tdl<-unlist(c(GPCR_tdl_labels$tdl[31:60],
                           GPCR_tdl_labels$tdl[61:90],
                           GPCR_tdl_labels$tdl[121:150],
                           GPCR_tdl_labels$tdl[151:180],
                           GPCR_tdl_labels$tdl[181:210],
                           GPCR_tdl_labels$tdl[211:240],
                           GPCR_tdl_labels$tdl[241:270],
                           GPCR_tdl_labels$tdl[301:330],
                           GPCR_tdl_labels$tdl[331:360],
                           GPCR_tdl_labels$tdl[361:390],
                           GPCR_tdl_labels$tdl[391:406],
                           GPCR_tdl_labels$tdl[91:120],
                           GPCR_tdl_labels$tdl[271:300],
                           GPCR_tdl_labels$tdl[1:30]
))
nr_results$target_id<-nr_ids
nr_results$tdl<-unlist(nr_tdl_labels$tdl)

IC_results$target_id<-unlist(IC_ids)
IC_results$tdl<-unlist(IC_tdl_labels)

######################################################
# Combine all 5 classes/non-classes into one giant data frame:
#
# `all_tcrd`
#
######################################################

Rprof ( tf <- "log.log",  memory.profiling = TRUE )
all_tcrd<-bind_rows(kinase_results,gcpr_results)
all_tcrd<-bind_rows(all_tcrd,IC_results)
all_tcrd<-bind_rows(all_tcrd,nr_results)
all_tcrd<-bind_rows(all_tcrd,null_results)
Rprof ( NULL ) ; print ( summaryRprof ( tf )  )


######################################################
# Zero all NA values. 
# Because we are checking for 'presence' of an association
# and not obtaining any quantitative values, we felt that
# 0ing the features not listed was a safe assumption.
######################################################
Rprof ( tf <- "log.log",  memory.profiling = TRUE )


######################################################
# R memory limit forces us to fragment this process
# as well.
######################################################
all_tcrd_1_300<-all_tcrd[1:300,]
all_tcrd_301_600<-all_tcrd[301:600,]
all_tcrd_601_900<-all_tcrd[601:900,]
all_tcrd_901_1200<-all_tcrd[901:1200,]
all_tcrd_1201_1500<-all_tcrd[1201:1500,]
all_tcrd_1501_1800<-all_tcrd[1501:1800,]
all_tcrd_1801_2000<-all_tcrd[1801:2000,]
all_tcrd_2001_2117<-all_tcrd[2001:2117,]

zeroer<-function(df){
  # Temporary store categorical variable like TDL
  
  temp_tdls<-unlist(df$tdl)
  
  mat<-data.matrix(df)
  mat[is.na(mat)]<-0
  df<-data.frame(mat)
  # Reattach target_ids and tdl
  
  df$tdl<-temp_tdls
  return(df)
}

all_tcrd_1_300<-zeroer(all_tcrd_1_300)
all_tcrd_301_600<-zeroer(all_tcrd_301_600)
all_tcrd_601_900<-zeroer(all_tcrd_601_900)
all_tcrd_901_1200<-zeroer(all_tcrd_901_1200)
all_tcrd_1201_1500<-zeroer(all_tcrd_1201_1500)
all_tcrd_1501_1800<-zeroer(all_tcrd_1501_1800)
all_tcrd_1801_2000<-zeroer(all_tcrd_1801_2000)
all_tcrd_2001_2117<-zeroer(all_tcrd_2001_2117)

all_tdl_class<-all_tcrd$tdl

all_tcrd_zeroed<-rbind(all_tcrd_1_300,all_tcrd_301_600)
all_tcrd_zeroed<-rbind(all_tcrd_zeroed,all_tcrd_601_900)
all_tcrd_zeroed<-rbind(all_tcrd_zeroed,all_tcrd_901_1200)
all_tcrd_zeroed<-rbind(all_tcrd_zeroed,all_tcrd_1201_1500)
all_tcrd_zeroed<-rbind(all_tcrd_zeroed,all_tcrd_1501_1800)
all_tcrd_zeroed<-rbind(all_tcrd_zeroed,all_tcrd_1801_2000)
all_tcrd_zeroed<-rbind(all_tcrd_zeroed,all_tcrd_2001_2117)

all_tcrd<-all_tcrd_zeroed
rm(all_tcrd_zeroed)

Rprof ( NULL ) ; print ( summaryRprof ( tf )  )

######################################################
# Combine all the gold standards, extracted from
# WITHDARWN database.
######################################################
all_gold<-read.csv('toxic_targets/GoldStandards.csv')


######################################################
# The caveat to chunking queries based on feature set
# Size is that the order of GPCR IDs get messed up
######################################################

all_tcrd$target_id<-c(kinase_ids,unlist(c(GPCR_ids[31:60],
                                          GPCR_ids[61:90],
                                          GPCR_ids[121:150],
                                          GPCR_ids[151:180],
                                          GPCR_ids[181:210],
                                          GPCR_ids[211:240],
                                          GPCR_ids[241:270],
                                          GPCR_ids[301:330],
                                          GPCR_ids[331:360],
                                          GPCR_ids[361:390],
                                          GPCR_ids[391:406],
                                          GPCR_ids[91:120],
                                          GPCR_ids[271:300],
                                          GPCR_ids[1:30]
)),
IC_ids,
nr_ids,
null_ids)

all_tcrd_reduced<-all_tcrd[,-c(which(colnames(all_tcrd)=="tdl"))]

# Select for training sets, based on our manually minined gold standard list
# This is our "training/testing" set

all_tcrd_reduced<-all_tcrd_reduced[all_tcrd_reduced$target_id %in% all_gold$id,]
all_tcrd_reduced_target_id<-unlist(all_tcrd_reduced$target_id)
