\name{CV_filter}
\alias{CV_filter}
\title{Filter metabolic data by CV}
\description{
This function allows filtering a matrix of metabolic variables based on the coefficient of variation (CV) of each variable across the quality control (QC) samples. See also function "QC_CV()".
}
\usage{
CV_filter(metabo_matrix, CV_metabo, CV_th = 0.30)
}
\arguments{
 \item{metabo_matrix}{numeric matrix containing the metabolic data (e.g. NMR peak
  intensities or metabolite concentrations). The columns of the matrix must correspond
  to the metabolic variables and the rows to the samples.}
  \item{CV_metabo}{numeric vector containing the CVs of the metabolic variables.
  See function "QC_CV()".}
  \item{CV_th}{numeric value indicating the CV threshold. Only features with CV below
  CV_th will be retained in the matrix.}
}
\value{
A numeric matrix corresponding to the CV-filtered metabo_matrix.
}
\references{
Dumas ME, et al. (2006). Assessment of analytical reproducibility of 1H NMR spectroscopy
based metabonomics for large-scale epidemiological research: the INTERMAP Study.
Analytical Chemistry, 78, 2199-1208.
}
\examples{
## Load data
data(metabolic_data)

## Select subset of QC samples
QC_data <- metabolic_data[507:516,]

## Calculate CVs
CV_metabo <-  QC_CV(QC_data)

## Filter metabolic_data by CV
metabo_CVfiltered <- CV_filter(metabolic_data, CV_metabo, CV_th = 0.30)
metabo_CVfiltered2 <- CV_filter(metabolic_data, CV_metabo, CV_th = 0.15)
}
