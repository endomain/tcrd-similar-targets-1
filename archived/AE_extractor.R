############################################################
# First looking at all adverse events in clinical trials
############################################################

ct_ae<-read.csv('drug_data/SingleSeriousAEdrugs_orderbyAEProb+StudySize.csv')

includer<-function(df,items){
  results_df<-data.frame()
  df[grepl(kinase_grep,df$Treatment,ignore.case=TRUE)|grepl(kinase_grep,df$Treatment_description,ignore.case=TRUE),]
}

kinase_grep<-"cediranib|enzastaurin|necitumumab|alectinib|cobimetinib|trastuzumab|emtansine|palbociclib|lenvatinib|linaclotide|ceritinib|masitinib|cabozantinib|ponatinib|trametinib|sunitinib|dabrafenib|erlotinib|dasatinib|fasudil|ibrutinib|idelalisib|nesiritide|palifermin|nintedanib|regorafenib|ingenol mebutate|alvocidib|ramucirumab|becaplermin|ripasudil|pazopanib|gefitinib|sorafenib|imatinib|nilotinib|lapatinib|cetuximab|panitumumab|pertuzumab|vandetanib|vemurafenib|crizotinib|ruxolitinib|axitinib|afatinib|bosutinib|tofacitinib|trastuzumab|ridaforolimus|mecasermin|osimertinib|quercetin|ruboxistaurin"
test<-includer(ct_ae,kinase_grep)
############################################################
# 2) clinical trials with AEs as a cause of halting trials
############################################################
phase12_ct_ae<-read.csv('drug_data/TerminatedPhase12_due to AE or efficacy.csv')
phase12_ct_ae<-includer(phase12_ct_ae,kinase_grep)

#
# given a our filtered file, looks for 
# drug names in columns 'Treatment' and 'Treatment_description'
# 
drug_ct_ae_counter<-function(df,kinase_names){
  result_df<-data.frame()
  for(drug in kinase_names){
    result_df<-rbind(result_df,data.frame(drug=drug,
               count=sum(grepl(drug,phase12_ct_ae$Treatment,ignore.case = TRUE)|grepl(drug,phase12_ct_ae$Treatment_description,ignore.case = TRUE)))
    )
  }
  return(result_df)
}

ct_ae_drug_counts<-drug_ct_ae_counter(phase12_ct_ae,kinase_names) %>% filter(count>0)

#####################
# DrugBank
#####################

library(xml2)
drug_bank<-read_xml("drug_data/drug_bank.xml")
