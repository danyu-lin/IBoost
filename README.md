# **IBoost**

# **IBoost: An Integrative Boosting Approach for Predicting Survival Time With Multiple Genomics Platforms**

IBoost is an R function for predicting survival time using multiple types of potentially high-dimensional genomics and clinical data. It implements the I-Boost method of Wong *et al.* (2017) to estimate a sparse model for the survival time.

## **USAGE**

**IBoost** (**X, Y, data.type, method**=”permute”, **iter.max**=2000, **v=**ifelse(method**==**“permute”,0.1,1), **m.stop**=5, **alpha.series**=c(0.05,seq(0.1,1,by=0.1)), **n.fold**=5, **seed**=12345678)

## **ARGUMENTS**

| Argument         | Default                         | Description                                                                                                                                                                                                                                      |
|:-----------------|:--------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **X**            | {NONE}                          | The matrix of all predictors. Each row represents a subject (with total of *n* rows) and each column represents a feature (with total of *d* columns). Each column of **X** should be standardized.                                              |
| **Y**            | {NONE}                          | The survival time, represented by a Surv object.                                                                                                                                                                                                 |
| **data.type**    | {NONE}                          | The list of indices representing the types of the predictors. Each element of the list is a vector of integers between 1 and *d* that corresponds to the column numbers of a type of predictors in **X**. The indices should be non-overlapping. |
| **method**       | permute                         | The version of I-Boost to be used; “permute” for I-Boost-Permutation, and “CV” for I-Boost-CV.                                                                                                                                                   |
| **iter.max**     | 2000                            | The maximum number of iterations.                                                                                                                                                                                                                |
| **v**            | ifelse(method==”permute”,0.1,1) | The penalizing factor for the current estimate at each iteration.                                                                                                                                                                                |
| **m.stop**       | 5                               | The stopping criterion; if the parameters are not updated consecutively for **m.stop** iterations, then the algorithm terminates.                                                                                                                |
| **alpha.series** | c(0.05,seq(0.1,1,by=0.1))       | The set of values to be considered for the tuning parameter α in cross-validation; only applicable when **method**=”CV”.                                                                                                                         |
| **n.fold**       | 5                               | The number of folds in cross-validation for the selection of tuning parameters; only applicable when **method**=”CV”.                                                                                                                            |
| **seed**         | 12345678                        | The initial random seed for partitioning the data set for cross-validation or generating permutation data sets.                                                                                                                                  |

## **VALUE**

The function returns a list of two elements, **beta** and **iter.no**. **beta** is a vector of estimated regression parameters of the selected variables. **iter.no** is the number of iterations used in the estimation.

## **DATA AND EXAMPLE**

### **Data file**

The files “TCGA_8cancer_rmmis_1.csv” to “TCGA_8cancer_rmmis_6.csv” contain the pan-cancer data set from The Cancer Genome Atlas (TCGA) studied by Wong *et al.* (2017). (They should be one big data file but are split due to size constrains.) Each column of the combined data file gives data on the overall survival and different types of predictors of a patient. The survival (or censoring) time in months and the event indicator are given in the rows labeled “OS_Months_5yr” and “OS_Event_5yr”, respectively. The label of each predictor begins with the abbreviation of the data type, followed by an underscore “\_” and the name of the predictor. The abbreviations of data types are as follows: Gene represents raw gene expression; Module represents gene module; Clinical represents clinical variable; CNV represents copy number variant; Mutation represents somatic mutation; miRNA represents micro-RNA expression; and Protein represents protein expression. Gender is coded as female = 0 and male = 1. Pathologic stage T is dichotomized into T1 (0) and T2–T4 (1). Pathologic stage N is dichotomized into N0 (0) and N1–N3 (1).

The file “DataSplit.csv” indicates whether patients belong to the training or testing set in each of the 30 splits considered in Wong *et al.* (2017). The variable “splitx” indicates, by the values 0 versus 1, whether the patient belongs to the training or testing set in the *x*th split, respectively.

### **Example – Evaluation of LASSO, Elastic Net, or I-Boost Using TCGA Data**

The R script “main_program.R” performs the training and testing procedure on the TCGA data using LASSO, elastic net, or I-Boost, as described in Wong *et al.* (2017). The user should specify, at the beginning of the script, the analysis method, the subset of patients to be analyzed, the set of predictors, and the training and testing data split by setting the values of the following variables:

| Variable    | Description                                                                                                                                      |
|:------------|:-------------------------------------------------------------------------------------------------------------------------------------------------|
| **method**  | The analysis method; “permute” = I-Boost-Permutation; “CV” = I-Boost-CV; “lasso” = LASSO; “en” = elastic net.                                    |
| **type**    | The subset of patients to be analyzed; 1=LUAD patients; 2=KIRC patients; 3=all pan-cancer patients.                                              |
| **predset** | The index for the set of predictors, which ranges from 1 to 95; see the file “config.csv” for the set of predictors that each number represents. |
| **s**       | The split of training and testing data sets to be adopted. For s=1, 2, …, 30, the *s*th split in the file “DataSplit.csv” is used.               |

To run the script, the files “boosting_functions.R”, “DataSplit.csv”, and “TCGA_8cancer_rmmis_1.csv” to “TCGA_8cancer_rmmis_6.csv” should be present in the home directory. The packages “Hmisc”, “glmnet”, “R.utils”, “survival”, and “data.table” should be installed.

The R script generates three output files, with names “summary-”, “model-”, or “predict-” followed by labels of the set of predictors used, subset of patients analyzed, analysis method used, and the training and testing data split adopted. For example, the analysis results generated using I-Boost-CV on the LUAD data set with clinical and protein expression data as predictors and the 1^st^ training and testing data split are presented in the files “summary-26.Clinical.Protein-LUAD-1-CV.csv”, “model-26.Clinical.Protein-LUAD-1-CV.csv”, and “predict-26.Clinical.Protein-LUAD-1-CV.csv”.

The “summary-” file presents:

| Row name           | Description                                                                                                                                                                                                                                                                                                 |
|:-------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Model              | The set of predictors used.                                                                                                                                                                                                                                                                                 |
| alpha/Iteration No | The value of the tuning parameter *α* for method = “lasso” or “en”, or the number of iterations used for method = “CV” or “permute”.                                                                                                                                                                        |
| lambda             | The value of the tuning parameter *λ* for method = “lasso” or “en”; the value is NA for I-Boost.                                                                                                                                                                                                            |
| Total Size         | Total number of variables selected in the final model.                                                                                                                                                                                                                                                      |
| GeneExp Size       | Number of gene expression variables selected in the final model; is NA if gene expression data are not used as predictors.                                                                                                                                                                                  |
| Module Size        | Number of gene modules selected in the final model; is NA if gene modules are not used as predictors.                                                                                                                                                                                                       |
| Clinical Size      | Number of clinical variables selected in the final model; is NA if clinical variables are not used as predictors.                                                                                                                                                                                           |
| CNV Size           | Number of copy number variables selected in the final model; is NA if copy number data are not used as predictors.                                                                                                                                                                                          |
| Mutation Size      | Number of somatic mutations selected in the final model; is NA if somatic mutations are not used as predictors.                                                                                                                                                                                             |
| miRNA Size         | Number of miRNA expression variables selected in the final model; is NA if miRNA expression data are not used as predictors.                                                                                                                                                                                |
| Protein Size       | Number of protein expression variables selected in the final model; is NA if protein expression data are not used as predictors.                                                                                                                                                                            |
| Train N            | Number of patients in the training set.                                                                                                                                                                                                                                                                     |
| Train Log Rank P   | Patients are classified into risk groups of high risk versus low risk according to whether the estimated risk score is higher than or lower than the median risk score in the training set. This entry presents the p-value of the log-rank test between the two groups among patients in the training set. |
| Train HR           | The estimated log-hazard ratio of high risk versus low risk (defined in the above entry) under a Cox proportional hazards model among patients in the training set.                                                                                                                                         |
| Train COX P        | The p-value for the risk group under the Cox model specified in the above entry.                                                                                                                                                                                                                            |
| Train C-Index      | The estimated C-index for the risk score among patients in the training set.                                                                                                                                                                                                                                |
| Test N             | Number of patients in the testing set.                                                                                                                                                                                                                                                                      |
| Test Log Rank P    | Patients are classified into risk groups of high risk versus low risk according to whether the estimated risk score is higher than or lower than the median risk score a the testing set. This entry presents the p-value of the log-rank test between the two groups among patients in the testing set.    |
| Test HR            | The estimated log-hazard ratio of high risk versus low risk (defined in the above entry) under a Cox proportional hazards model among patients in the testing set.                                                                                                                                          |
| Test COX P         | The p-value for the risk group under the Cox model specified in the above entry.                                                                                                                                                                                                                            |
| Test C-Index       | The estimated C-index for the risk score among patients in the testing set.                                                                                                                                                                                                                                 |

The “model-” file presents the selected variables and their estimated regression parameters. Here is an example of a “model-” file:

|                       |                 |
|-----------------------|-----------------|
| **Feature**           | **Coefficient** |
| Protein_Smad4         | -0.8613         |
| Protein_Shc_pY317     | -0.5928         |
| Protein_Chk2_pT68     | -0.1763         |
| Clinical_age          | 0.011           |
| Protein_Cyclin_B1     | 0.2522          |
| Clinical_pathologic_N | 0.8538          |

The “predict-” file presents the risk scores. Here is an example of the top few lines of a “predict-” file:

|              |                 |          |
|--------------|-----------------|----------|
| **ID**       | **Testing Set** | **Risk** |
| TCGA.05.4249 | 0               | -1.119   |
| TCGA.05.4250 | 0               | 1.5743   |
| TCGA.05.4384 | 1               | -0.2918  |
| TCGA.05.4396 | 0               | 0.6906   |
| TCGA.05.4397 | 0               | 1.8356   |
| TCGA.05.4398 | 1               | 0.6884   |
| TCGA.05.4405 | 1               | -0.309   |
| TCGA.05.4417 | 0               | -0.5835  |
| TCGA.05.4418 | 0               | 1.247    |

The variable ID denotes the TCGA barcode. The variable Testing Set, by values 1 versus 0, denotes whether the patient belongs to the testing set or the training set, respectively. The variable risk denotes the estimated risk score.

## **REFERENCE**

Wong, K. Y., Fan, C., Tanioka, M., Parker, J. S., Nobel, A. B., Zeng, D., Lin, D. Y., and Perou, C. M. (2017) An Integrative Boosting Approach for Predicting Survival Time With Multiple Genomics Platforms. *Submitted*.

## **DOWNLOAD**

#### **R script containing the function “IBoost” \[updated November 20 2017\]**

executable (zip archive) **»** [boosting_functions.R](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/boosting_functions.zip)

#### **R script for evaluation of variable selection methods using TCGA \[updated November 20 2017\]**

executable (zip archive) **»** [main_program.R](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/main_program.zip)

#### **Data files \[updated November 20 2017\]**

zip archive **»** [TCGA_8cancer_rmmis_1.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/TCGA_8cancer_rmmis_1.zip)

zip archive **»** [TCGA_8cancer_rmmis_2.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/TCGA_8cancer_rmmis_2.zip)

zip archive **»** [TCGA_8cancer_rmmis_3.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/TCGA_8cancer_rmmis_3.zip)

zip archive **»** [TCGA_8cancer_rmmis_4.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/TCGA_8cancer_rmmis_4.zip)

zip archive **»** [TCGA_8cancer_rmmis_5.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/TCGA_8cancer_rmmis_5.zip)

zip archive **»** [TCGA_8cancer_rmmis_6.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/TCGA_8cancer_rmmis_6.zip)

zip archive **»** [index_and_split.zip](http://dlin.web.unc.edu/wp-content/uploads/sites/1568/2017/11/index_and_split.zip)

## **VERSION HISTORY**

| Version | Date      | Description             |
|:--------|:----------|:------------------------|
| 1.0     | Nov. 2017 | First version released. |
