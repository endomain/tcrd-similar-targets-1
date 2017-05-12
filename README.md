## tcrd-similar-targets

*Jing Ai* and *Jung Hoon Son*

**Utilizing TCRD database (used by Pharos) to predict drug targets associated with adverse drug events. **

## Code Description:

------------

## Step1_SQL_query_matrix_builder.R 

This code queries each druggable target + additional targets identified via WITHDRAWN database from a local instance of TCRD MySQL database. It expands the feature set with each new target vector representation added. 

------------

## Step2: Feature selection 
### DataMappingColsums.R
The file contains code for feature reduction. Since a large number of features contained zeros, we first removed the columns with less than 20 non-zero values (out of 1129 samples). 

### FeatureSelection.ipynb
This file contains code for feature selection codes (in python). 
We applied chi-square feature selection and filered out features with chi-square p-values of <0.05

------------

## Step3 Predictive modeling of drug targets
### Combined_analysis.R
The files contains the modeling codes for L1-Logistic Regression (5-fold CV), Naive Bayes, Random Forrest based on cross validation.  
### SVM.R
The file contains the modeling code for Support Vector Machine (Linear and Kernel). 

------------

## Step4_Future work (biologicall validating the predicted targets)

### predicting_drugs.R
This code can be used to generate an adverse drug scoring system.  We query the entire TCRD database for drugs associated with targets (among the 2117 we predicted) predicted by our algorithm. We obtain a score using `# of Predicted Adverse Targets` - `# of Predicted non-Adverse targets`. 

Ones with higest score we expect to have highest likelihood of being implicated with an adverse drug effect. 

### Plot_playground.R
This is an experimental file, used to generate exploratory heatmaps and network visualizations. 
