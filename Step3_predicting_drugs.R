library(caTools)
library(e1071)
library(randomForest)
library(caret)
library(RMySQL)
library(dplyr)
data_viz_df<-all_tcrd
data_viz_df$idgfam<-0
temp_predict_all<-predict(model,data_viz_df)

target_adverse<-data.frame(target_id=all_tcrd$target_id,
             adverse=unlist(temp_predict_all))

# Get positive from drug_activity in PHAROS
mydb = dbConnect(MySQL(), user='root', password='jhoon11', dbname='pharos', host='localhost')
qry<-paste("SELECT drug,count(drug) AS n from drug_activity 
           WHERE target_id IN (",paste(target_adverse[target_adverse$adverse==1,]$target_id, collapse=","),") GROUP BY drug;"  )
drug_pred_counts<-dbSendQuery(mydb,qry)
drug_pred_counts<-fetch(drug_pred_counts, n=-1)
drug_pred_counts<-drug_pred_counts%>%
  arrange(desc(n))

# Get negative from drug_activity in PHAROS
qry2<-paste("SELECT drug,count(drug) AS n_minus from drug_activity 
           WHERE target_id IN (",paste(target_adverse[target_adverse$adverse==0,]$target_id, collapse=","),") GROUP BY drug;"  )
n_minus<-dbSendQuery(mydb,qry2)
n_minus<-fetch(n_minus, n=-1)
n_minus<-n_minus%>%
  arrange(desc(n_minus))

# merge two queries
drug_pred_counts<-merge(drug_pred_counts,n_minus,by='drug')
drug_pred_counts<-drug_pred_counts %>% 
  mutate(diff=(n-n_minus))
# 
# # Apply it on out CT AE drug counts
# ct_ae_drug_counts<-merge(ct_ae_drug_counts, drug_pred_counts, by = "drug", all.x = TRUE)
# ct_ae_drug_counts %>% 
#   na.omit() %>%
#   arrange(desc(diff))


##################################
# Get positive from chembl in PHAROS
mydb = dbConnect(MySQL(), user='root', password='jhoon11', dbname='pharos', host='localhost')
qry<-paste("SELECT distinct(drug) FROM drug_activity 
           WHERE target_id IN (",paste(target_adverse[target_adverse$adverse==1,]$target_id, collapse=","),") GROUP BY drug"  )
qry<-paste("SELECT distinct(target_id) FROM chembl_activity 
           WHERE cmpd_name_in_ref IN (",qry,")")
qry<-paste("SELECT drug,count(drug) AS n from drug_activity 
           WHERE target_id IN (",qry,") group by drug;")


drug_pred_counts<-dbSendQuery(mydb,qry)
drug_pred_counts<-fetch(drug_pred_counts, n=-1)
drug_pred_counts<-drug_pred_counts%>%
  arrange(desc(n))


# Get negative from chembl in PHAROS
qry2<-paste("SELECT distinct(drug) FROM drug_activity 
           WHERE target_id IN (",paste(target_adverse[target_adverse$adverse==0,]$target_id, collapse=","),") GROUP BY drug"  )
qry2<-paste("SELECT distinct(target_id) FROM chembl_activity 
           WHERE cmpd_name_in_ref IN (",qry2,")")
qry2<-paste("SELECT drug,count(drug) AS n_minus from drug_activity 
           WHERE target_id IN (",qry2,") group by drug;")


n_minus<-dbSendQuery(mydb,qry2)
n_minus<-fetch(n_minus, n=-1)
n_minus<-n_minus%>%
  arrange(desc(n_minus))

# merge two queries
drug_pred_counts<-merge(drug_pred_counts,n_minus,by='drug',all=TRUE)
drug_pred_counts[is.na(drug_pred_counts)]<-0

drug_pred_counts<-drug_pred_counts %>% 
  mutate(score=(n-n_minus))
# 
# # Apply it on our CT AE drug counts
# ct_ae_drug_counts<-merge(ct_ae_drug_counts, drug_pred_counts, by = "drug", all.x = TRUE)
# ct_ae_drug_counts %>% 
#   na.omit() %>%
#   arrange(desc(diff))