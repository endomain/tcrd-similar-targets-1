# tcrd-similar-targets

*Jing Ai* and *Jung Hoon Son*

__BINF G4002 Computation Methods__

__Final Project__

**Utilizing TCRD database (used by Pharos) to predict drug targets associated with adverse drug events. **

------------

## Compiled Data Set:

#### Main Data File: `all_tcrd.Rda`

R workspace file, containing `all_tcrd` data frame for our analysis. This data frame can be generated from `Step1_SQL_query_matrix_builder.R` file on any computer with a local instance of Pharos/TCRD database installed on MySQL. SQL data dump is available from:  (http://juniper.health.unm.edu/tcrd/download/tcrd_v4.4.2.sql.gz)

Individual IDG family based data sets (GPCR, Ion channels, Kinases, NRs, null-classes) have been compiled and is made available in the folder `individual_class`.

### 

------------

## Code Description:

------

### Step1: SQL Query, Vector represention of targets, and Working Matrix

##### `Step1_SQL_query_matrix_builder.R`

This code queries each druggable target + additional targets identified via WITHDRAWN database from a local instance of TCRD MySQL database. It expands the feature set with each new target vector representation added. 

##### `toxic_targets/GoldStandards.csv`
This file contains curated information from WITHDRAWN database, implicating drug targets associated with withdrawn drugs around the world (rows where y=0). Together with the protein targets classified as Tclin in TCRD that are not associated with any withdrawn drugs (rows where y=1), they serves as our gold standards for denoting which targets have been implicated with drug adverse events and which have not.

------

### Step2: Feature selection 
##### `analysis_code/DataMappingColsums.R`
The file contains code for feature reduction. Since a large number of features contained zeros, we first removed the columns with less than 20 non-zero values (out of 1129 samples). 

##### `analysis_code/FeatureSelection.ipynb`
This file contains code for feature selection codes (in python). 
We applied chi-square feature selection and filered out features with chi-square p-values of <0.05

`Step2_prelim_prediction.R` is used for preliminary `randomForest` model in R, for quickly testing our test sets.

------

### Step 3: Predictive modeling of drug targets
##### `analysis_code/combined_analysis_JA.R`
The files contains the modeling codes for L1-Logistic Regression (5-fold CV), Naive Bayes, Random Forrest based on cross validation.  
##### `analysis_code/SVM_JA.R`
The file contains the modeling code for Support Vector Machine (Linear and Kernel). 

------

### Step 4: Future work (biologically validating the predicted targets)

##### `future_work/predicting_drugs.R`
This code can be used to generate an adverse drug scoring system.  We query the entire TCRD database for drugs associated with targets (among the 2117 we predicted) predicted by our algorithm. We obtain a score using 

> Score = `# of Predicted Adverse Targets` - `# of Predicted non-Adverse targets`. 

In which drugs with higest score we expect to have highest likelihood of being implicated with an adverse drug effect. Current

##### `future_work/Plot_playground.R`
This is an experimental file, used to generate exploratory heatmaps and network visualizations. 



