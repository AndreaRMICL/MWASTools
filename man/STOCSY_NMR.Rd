\name{STOCSY_NMR}
\alias{STOCSY_NMR}
\title{Statistical Total Correlation Spectroscopy - Academic use only}
\description{
This function calculates STOCSY between an NMR signal of interest and all the NMR
variables, representing a useful tool for NMR molecular identification and assignment.
The results are represented in a pseudo-NMR spectrum displaying the covariance (height) and the Pearson correlation coefficient (color) of all spectral variables with the variable of interest (driver signal).
}
\usage{
STOCSY_NMR (metabo_matrix, ppm, ppm_query, alpha_th = 0.05, xlab = "ppm",
            ylab = "covariance", size_lab = 12, size_axis = 12, xlim = NULL, ylim = NULL,
            xbreaks = waiver(), xnames = waiver(), ynames = waiver(), ybreaks = waiver())
}
\arguments{
  \item{metabo_matrix}{numeric matrix containing the metabolic data (e.g. NMR peak
  intensities or metabolite concentrations). The columns of the matrix must correspond
  to the metabolic variables and the rows to the samples.}
  \item{ppm}{numeric vector defining the chemical shift scale.}
  \item{ppm_query}{numeric value (at least 2 decimals) corresponding to the driver ppm.}
  \item{alpha_th}{numeric value indicating the significance threshold. NMR variables
  with BH-adjusted p-value above this threshold will be neglected.}
  \item{xlab}{character vector specifying a title for the x-axis.}
  \item{ylab}{character vector specifying a title for the y-axis. }
  \item{size_lab}{numeric value indicating the font size of x- and y- axis titles.}
  \item{size_axis}{numeric value indicating the font size of x- and y- axis labels. }
  \item{xlim}{numeric vector containing the minimum and maximum values of the x-axis.
  Notice that ppm is displayed in reverse scale (e.g. xlim = c(2, 1)).}
  \item{ylim}{numeric vector containing the minimum and maximum values of the y-axis.}
  \item{xbreaks}{numeric vector indicating the positions of the breaks of the x-axis.}
  \item{xnames}{character vector (same length as xbreaks) containing the labels
  of each break of the x-axis.}
  \item{ybreaks}{numeric vector indicating the positions of the breaks of the y-axis.}
  \item{ynames}{character vector (same length as ybreaks) containing the labels
  of each break of the y-axis.}
}
\value{
A plot displaying the Pearson correlation coefficient (color) and covariance (height)
between all spectral variables and the driver signal.
}
\references{
Cloarec O, et al.(2005). Statistical total correlation spectroscopy: An exploratory approach for latent biomarker identification from metabolic 1H NMR data sets. Analytical Chemistry, 77, 1282–1289.

Nicholson JK, et al. (2005). Method for the identification of molecules and biomarkers using chemical, biochemical and biological data. US 20070043518 A1
}
\examples{
## Load data
data(metabolic_data)

## Define ppm variable and ppm_query
ppm <- as.numeric(colnames(metabolic_data))
ppm_query <- 1.04 # valine signal

## STOCSY
STOCSY_NMR (metabolic_data, ppm, ppm_query)
STOCSY_NMR (metabolic_data, ppm, ppm_query, alpha_th = 0, xlim = c(1.06, 1))
}

