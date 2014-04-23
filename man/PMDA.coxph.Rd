\name{PMDA.coxph}
\alias{PMDA.coxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple Imputation using Poor Man's Data Augmentation scheme (PMDA) for Cox proportional hazard regression}
\description{for internal uses}
\usage{
PMDA.coxph(formula, data, imax = 25, k = 10, th0 = 1e-03)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{formula}{a formula}
  \item{data}{data used in this function}
  \item{imax}{maximum number of iteration}
  \item{k}{number of dataset to generate for each iteration}
  \item{th0}{note used currently}
  }

\details{
This function is called by MIICD.coxph with option PMDA (Poor Man's Data Augmentation)
}
\value{
\item{beta}{estimation of the coefficient(s)}
\item{sr_sigma}{estimation of the standard errors}
\item{n_iter}{number of iteration used}
\item{sigma1}{mean of standard errors not augmented for the between imputations composant of the variance}
}

\references{
PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1, p. 199-203.

Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. 

Schenker, N. and Welsh, A. (1988). Asymptotic results for multiple imputation. The Annals of Statistics pages 1550-1566.

Tanner, M. A. and Wong, W. H. (1987). An application of imputation to an estimation problem in grouped lifetime analysis. Technometrics 29, 23-32.

Wei, G. C., & Tanner, M. A. (1991). Applications of multiple imputation to the analysis of censored regression data. Biometrics, 47(4), 1297-1309.
}
\author{
Marc Delord <mdelord@gmail.com>
}
%\note{}
\seealso{PMDA.crreg}
%\examples{}
\keyword{Interval Censored Data}
\keyword{Multiple Imputation}
\keyword{Poor Man's Data Augmentation}
\keyword{Cox Proportional Hazards Regression}
