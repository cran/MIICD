\name{ANDA.crreg}
\alias{ANDA.crreg}
\title{Multiple Imputation using Asymptotic Normal Data Augmentation scheme (ANDA) for Competing risks regression model}
\description{For internal use}
\usage{
ANDA.crreg(formula, data, imax = 25, k = 10, th0 = 1e-03,
status, trans, cens, keep, model = "FG")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula}
  \item{data}{a data frame}
  \item{imax}{maximum number of iteration}
  \item{k}{number of dataset generated at each iteration}
  \item{th0}{note used currently}
  \item{status}{name of the column were the status of observation is indicated}
  \item{trans}{what is the transition of interest in the status column}
  \item{cens}{code used for right censored data in the status column}
  \item{keep}{must be the same as the covariate used in formula}
  \item{model}{which model shall be used currently FG state for Fine and Gray proportional hazards regression model}
}
\details{
This function is called by MIICD.crreg with option ANDA (Asymptotic Normal Data Augmentation)
}
\value{
%beta=beta_iter,sr_sigma=sigma_iter,n_iter=i,conv=th_iter,sigma1=sigma1_iter
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
\code{\link{PMDA.crreg}}
\code{\link{crr}}
}
\keyword{Interval Censored Data}
\keyword{Multiple Imputation}
\keyword{Asymptotic Normal Data Augmentation}
\keyword{Competing Risks Regression}
