## tcrd-similar-targets

*Jing Ai* and *Jung Hoon Son*

**Utilizing TCRD database (used by Pharos) to predict drug targets associated with adverse drug events. **

## Compiled Data Set:

`

------------

## Code Description:

------------

## Step1: SQL Query, Vector represention of targets, and Working Matrix

#### `Step1_SQL_query_matrix_builder.R`

This code queries each druggable target + additional targets identified via WITHDRAWN database from a local instance of TCRD MySQL database. It expands the feature set with each new target vector representation added. 

------------

## Step2: Feature selection 
#### `analysis_code/DataMappingColsums.R`
The file contains code for feature reduction. Since a large number of features contained zeros, we first removed the columns with less than 20 non-zero values (out of 1129 samples). 

#### `analysis_code/FeatureSelection.ipynb`
This file contains code for feature selection codes (in python). 
We applied chi-square feature selection and filered out features with chi-square p-values of <0.05

`Step2_prelim_prediction.R` is used for preliminary `randomForest` model in R, for quickly testing our test sets.

------------

## Step 3: Predictive modeling of drug targets
### `analysis_code/combined_analysis_JA.R`
The files contains the modeling codes for L1-Logistic Regression (5-fold CV), Naive Bayes, Random Forrest based on cross validation.  
### `analysis_code/SVM_JA.R`
The file contains the modeling code for Support Vector Machine (Linear and Kernel). 

------------

## Step 4: Future work (biologically validating the predicted targets)

#### `future_work/predicting_drugs.R`
This code can be used to generate an adverse drug scoring system.  We query the entire TCRD database for drugs associated with targets (among the 2117 we predicted) predicted by our algorithm. We obtain a score using `# of Predicted Adverse Targets` - `# of Predicted non-Adverse targets`. 

Ones with higest score we expect to have highest likelihood of being implicated with an adverse drug effect. 

#### `future_work/Plot_playground.R`
This is an experimental file, used to generate exploratory heatmaps and network visualizations. 
