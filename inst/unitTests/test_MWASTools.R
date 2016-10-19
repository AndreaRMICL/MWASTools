######################### TEST MWASTools ###################################
library(RUnit)

## Load data
data(metabolic_data)
data(target_metabolites)
data(clinical_data)

## Define variables
BMI = clinical_data[,4]
diabetes = clinical_data[,3]
Age_Gender = clinical_data[,1:2]
ppm = as.numeric(colnames(metabolic_data))
metabolites = colnames(target_metabolites)
QC_data = metabolic_data[507:516,]

#### Test MWAS_stats ####

## Test for association between diabetes and target_metabolites
T2D_model1 = MWAS_stats (target_metabolites, diabetes,
                         assoc_method = "logistic", mt_method = "BH",
                         output = "pvalues", metabo_ids = metabolites)

T2D_model2 = MWAS_stats (target_metabolites, diabetes,
                         assoc_method = "logistic", mt_method = "BY",
                         output = "pvalues", metabo_ids = metabolites)

T2D_model3 = MWAS_stats (target_metabolites, diabetes,
                         assoc_method = "logistic", mt_method = "bonferroni",
                         output = "pvalues", metabo_ids = metabolites)

T2D_model4 = MWAS_stats (target_metabolites, diabetes,
                         assoc_method = "logistic", mt_method = "none",
                         output = "pvalues", metabo_ids = metabolites)

## Test for association between diabetes and target_metabolites (age & gender adjusted)
T2D_model5 = MWAS_stats (target_metabolites, diabetes, CF_matrix = Age_Gender,
                         assoc_method = "logistic", mt_method = "BH",
                         output = "pvalues", metabo_ids = metabolites)

## Test for association between BMI and target_metabolites
BMI_model1 = MWAS_stats (target_metabolites, BMI,
                        assoc_method = "spearman", mt_method = "BH",
                        output = "pvalues", metabo_ids = metabolites)

## Test for association between BMI and target_metabolites (age & gender adjusted)
BMI_model2 = MWAS_stats (target_metabolites, BMI, CF_matrix = Age_Gender,
                        assoc_method = "spearman", mt_method = "BH",
                        output = "pvalues", metabo_ids = metabolites)

## Test for association between BMI and target_metabolites (output = models)
BMI_model3 = MWAS_stats (target_metabolites, BMI, CF_matrix = Age_Gender,
                         assoc_method = "spearman", mt_method = "BH",
                         output = "models", metabo_ids = metabolites)

## Attempt to do logistic regression with a non-binary predictive variable
BMI_model4 = try(MWAS_stats (target_metabolites, BMI, assoc_method = "logistic"), silent = TRUE)

## Attempto to apply MWAS_stats on variables with inconsistent dimension
BMI_model5 = try(MWAS_stats (target_metabolites, head(BMI), assoc_method = "spearman"), silent = TRUE)

## Attempto to apply MWAS_stats with invalid assoc_method
BMI_model6 = try(MWAS_stats (target_metabolites, BMI, assoc_method = "X"), silent = TRUE)

## Attempto to apply MWAS_stats with invalid assoc_method
BMI_model7 = try(MWAS_stats (target_metabolites, BMI, assoc_method = "spearman",
                             mt_method = "X"), silent = TRUE)


## Tests
test.MWAS_stats <- function() {

    checkTrue(is.matrix(T2D_model1))
    checkTrue(is.matrix(T2D_model2))
    checkTrue(is.matrix(T2D_model3))
    checkTrue(is.matrix(T2D_model4))
    checkTrue(is.matrix(T2D_model5))
    checkTrue(is.matrix(BMI_model1))
    checkTrue(is.matrix(BMI_model2))
    checkTrue(is.list(BMI_model3))
    checkException(BMI_model4)
    checkException(BMI_model5)
    checkException(BMI_model6)
    checkException(BMI_model7)
}
test.MWAS_stats()

################################################################################

#### Test QC_CV ####
metabo_CV1 =  QC_CV (QC_data, metabo_ids = ppm)
metabo_CV2 =  QC_CV (QC_data)

## Tests
test.metabo_CV <- function() {
  checkTrue(is.vector(metabo_CV1))
  checkTrue(is.vector(metabo_CV2))
}
test.metabo_CV()

################################################################################

#### Test QC_PCA ####
PCA_model =  QC_PCA (target_metabolites)

## Tests
test.PCA_model <- function() {
  checkTrue(is.list(PCA_model))
}
test.PCA_model()

################################################################################

#### Test STOCSY_NMR ####

## Attempt to run STOCSY with a ppm_query outside the ppm range
STOCSY = try(STOCSY_NMR(metabolic_data, ppm, ppm_query = 10.01), silent = TRUE)

## Tests
test.STOCSY <- function() {
  checkException(STOCSY)
}
test.STOCSY()

################################################################################

#### Test MWAS_filter ####

## Test for association between diabetes and target_metabolites
T2D_model1 = MWAS_stats (target_metabolites, diabetes,
                         assoc_method = "logistic", mt_method = "BH",
                         output = "pvalues", metabo_ids = metabolites)

## Filter T2D_model1 by pvalue
pvalue_filter1 <- MWAS_filter(T2D_model1, type = "pvalue", alpha_th = 0.001)

## Attempt to filter T2D_model1 with too stringent criteria
pvalue_filter2 <- try(MWAS_filter(T2D_model1, type = "pvalue", alpha_th = 0), silent = TRUE)

## Tests
test.MWAS_filter <- function() {
  checkTrue(is.matrix(pvalue_filter1))
  checkException(pvalue_filter2)
}
test.MWAS_filter()

################################################################################

#### Test CV_filter ####

## Calculate CVs
CV_metabo <-  QC_CV (QC_data)

## Filter metabolic_data by CV
metabo_CVfiltered <- CV_filter(metabolic_data, CV_metabo, CV_th = 0.30)
metabo_CVfiltered2 <- CV_filter(metabolic_data, CV_metabo, CV_th = 0.15)

## Tests
test.CV_filter <- function() {
  checkTrue(is.matrix(metabo_CVfiltered) & is.matrix(metabo_CVfiltered2))
}
test.CV_filter()

