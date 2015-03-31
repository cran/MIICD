\name{MIICD}
\alias{MIICD}
\docType{package}
\title{Data Augmentation and Multiple Imputation for Interval Censored Data}
\description{This is version 2.0 of the MIICD package: bugs have been corrected and 2 sets of functions have been added to the previous MIICD functions. 
Proposed functions apply data augmentation and/or the multiple imputation approach of Pan (2000) to handle interval censored data with the Cox proportional hazards regression and interval censored competing risks data with the Fine & Gray proportional subdistribution hazards regression. 3 main sets of functions are proposed: \itemize{
\item \code{MI.surv} and \code{MI.ci}: these functions use the multiple imputation approach to compute the survival function and the cumulative incidence function with confidence intervals, for interval censored data and interval censored competing risks data respectively.  \code{print} and \code{plot} methods are available.
\item \code{DA.surv} and \code{DA.ci}: these functions use data augmentation schemes with multiple imputation at each iteration to compute the survival function and the cumulative incidence function with confidence intervals, for interval censored data and interval censored competing risks data respectively. \code{print} and \code{plot} methods are available.
\item \code{MIICD.coxph} and \code{MIICD.crreg}: both functions use data augmentation schemes with multiple imputation at each iteration to handle interval censored data with the Cox propotional hazards regression, and interval censored competing risks data with the Fine and Gray subdistribution hazard regression.
\itemize{
\item \code{MIICD.coxph} compute the baseline survival function and Cox regression parameters. 
\item \code{MIICD.rcreg} compute the baseline cumulative incidence function and Fine & Gray regression parameters.
}
For both functions, the \emph{Poor Man's Data augmentation} (\code{PMDA})  and the \emph{Asymptotic Normal Data Augmentation} (\code{ANDA}) shemes are implemented (by using the \code{methode} parameter) (Pan 2000; Wei and Tanner 1991; Tanner and Wong 1987a). \code{print} and \code{plot} methods are available.
}

MIICD functions call iteratively the \code{coxph} (library \code{survival}) or the \code{FGR} (library \code{riskRegression}) function on imputed datasets and derive multiple estimates from imputed data sets. Finally multiple estimates are combined following multiple imputation rules (Rubin (1987); Schenker and Welsh (1988); Tanner and Wong (1987b)) to update parameter estimates. When ANDA is chosen, \code{mvrnorm} (package \code{MASS}) is called.
}
\details{
\tabular{ll}{
Package: \tab MIICD\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-02-11\cr
License: \tab GPL-3 \cr
}

}
\author{
Marc Delord <mdelord@gmail.com>,
Maintainer: Marc Delord <mdelord@gmail.com>
}
\references{
PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1, p. 199-203.

Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. 

Schenker, N. and Welsh, A. (1988). Asymptotic results for multiple imputation. The Annals of Statistics pages 1550-1566.

Tanner, M. A. and Wong, W. H. (1987). An application of imputation to an estimation problem in grouped lifetime analysis. Technometrics 29, 23-32.

Wei, G. C., & Tanner, M. A. (1991). Applications of multiple imputation to the analysis of censored regression data. Biometrics, 47(4), 1297-1309.
}
\keyword{Competing Risks Interval Censored Data}
\keyword{Interval Censored Data}
\keyword{Multiple Imputation}
\keyword{Cox proportional hazards Regression}
\keyword{Fine and Gray proportional subdistribution hazards Regression}
\keyword{Poor Man's Data Augmentation}
\keyword{Asymptotic Normal Data Augmentation}

\seealso{
\code{\link{Surv}}
\code{\link{survfit}}
\code{\link{coxph}}
\code{\link{mvrnorm}}
}
\examples{
res <- MIICD.coxph(formula = ~ treatment , data = bcos,
k = 10, m = 10, method = 'ANDA', verbose = FALSE)
res
plot(res)
plot(res , type = 'coef' , coef = 1)
plot(res , type = 'sigma' , coef = 1)
}
