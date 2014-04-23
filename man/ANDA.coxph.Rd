\name{ANDA.coxph}
\alias{ANDA.coxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple Imputation using Asymptotic Normal Data Augmentation scheme (ANDA) for Cox proportional hazards model}
\description{For internal uses}
\usage{
ANDA.coxph(formula, data , imax = 25 , k = 10 , th0 = 1e-03 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula}
  \item{data}{a data frame}
  \item{imax}{maximum number of iteration}
  \item{k}{number of dataset generated at each iteration}
  \item{th0}{note used currently}
}
\details{
This function is called by MIICD.coxph with option ANDA (Asymptotic Normal Data Augmentation)
}
\value{
\item{beta}{estimation of the coefficient(s)}
\item{sr_sigma}{estimation of the standard errors}
\item{n_iter}{number of iteration used}
\item{sigma1}{mean of standard errors not augmented for the between imputations composant of the variance}
}
\references{
Pan, W. (2000). A multiple imputation approach to cox regression with interval-censored data. Biometrics 56, 199-203.

Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys.

Schenker, N. and Welsh, A. (1988). Asymptotic results for multiple imputation. The Annals of Statistics pages 1550-1566.

Tanner, M. A. and Wong, W. H. (1987). An application of imputation to an estimation problem in grouped lifetime analysis. Technometrics 29, 23-32.
}
\author{
Marc Delord (mdelord@gmail.com)
}
%\note{}
\seealso{
\code{\link{coxph}}
}
\keyword{Interval Censored Data}
\keyword{Multiple Imputation}
\keyword{Asymptotic Normal Data Augmentation}
\keyword{Cox Proportional Hazards Regression}