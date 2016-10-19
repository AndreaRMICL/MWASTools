\name{MWAS_bootstrapping}
\alias{MWAS_bootstrapping}
\title{MWAS bootstrap resampling}
\description{
This function generates bootstrap replicates (non-parametric resampling) of a model
testing for association between a given metabolite and a disease phenotype,
and calculates the confidence intervals of model coefficients.
}
\usage{
MWAS_bootstrapping (metabolite, disease, CF_matrix = NULL, assoc_method,
                    iterations = 10000)
}
\arguments{
  \item{metabolite}{numeric vector corresponding to the levels (e.g. concentration,
  peak intensity) of the metabolite of interest.}
  \item{disease}{numeric vector corresponding to the response to be modeled.
  For logistic regression, disease must be a two-value (e.g. 0,1) numeric vector.
  NA values are allowed.}
  \item{CF_matrix}{optional numeric matrix corresponding to the covariates to be
  included in the model ((e.g. age or gender)). The columns of the matrix must correspond
  to the covariates and the rows to the samples. NA values are allowed.}
  \item{assoc_method}{character constant indicating the association method that
  will be used. Possible values for assoc_method are: "pearson" (pearson correlation),
  "spearman" (spearman correlation), "kendall" (kendall correlation), "linear"
  (linear regression) or "logistic" (logistic regression).}
  \item{iterations}{numeric value indicating the number of bootstrap replicates}
}
\value{
A list with 3 elements, each list element reporting the following information:
i) object of class "boot"; ii) summary of the previous object; iii) 95-confidence
interval of the metabolite model coefficient. For more details, check the function
"boot()" from the "boot" package.
}
\references{
Davison AC, Hinkley, DV. (1997). Bootstrap Methods and Their Application. Cambridge
University Press.
}
\examples{
## Load data
data(target_metabolites)
data(clinical_data)

## Define variables
OHbutyrate <- target_metabolites[, 2]
diabetes <- clinical_data[, 3]

## Bootstrapping
MWAS_bootstrapping (OHbutyrate, diabetes, assoc_method = "logistic", iterations = 1000)
}

