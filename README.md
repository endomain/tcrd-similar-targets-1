## tcrd-similar-targets

*Jing Ai* and *Jung Hoon Son*

**Utilizing TCRD database (used by Pharos) to predict drug targets associated with poor safety profiles of drugs. **

## Code Description:

------------

## Step1_SQL_query_matrix_builder.R 

This code queries each druggable target + additional targets identified via WITHDRAWN database from a local instance of TCRD MySQL database. 

It expands the feature set with each new target vector representation added. It likely needs optimization, but for our purposes we ran it rarely. 

------------

## Step2_prelim_prediction.R

This file runs a very rudimentary training and prediction with `randomForest`. 
Simple feature reduction with column sums greater than a designated number is used, drastically reducing our feature set.

------------

## Step3_predicting_drugs.R

We use this to generate an adverse drug scoring system. 

We query the entire TCRD database for drugs associated with targets (among the 2117 we predicted) predicted by our algorithm. We obtain a score using `# of Predicted Adverse Targets` - `# of Predicted non-Adverse targets`. 

Ones with higest score we expect to have highest likelihood of being implicated with an adverse drug effect. 

------------

## Step4_Plot_playground.R

This is an experimental file, used to generate exploratory heatmaps and network visualizations. 