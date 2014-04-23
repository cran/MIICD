\name{MIICD}
\alias{MIICD}
\docType{package}
\title{
Multiple Imputation for proportional hazard regression with interval censored data.
}
\description{
This set of functions applies the multiple imputation approach of Pan (2000) to handle interval censored data using the Cox proportional hazards model, possibly with competing outcomes using the Fine and Gray Propotional Hazards model (or Cox Proportional Hazard Model). It calls iteratively the \code{coxph} (library \code{survival}) or the \code{crr} (library \code{cmprsk}) function on imputed datasets and derives multiple estimates from imputed data sets. Finally it combines multiple estimates following multiple some imputation rules (Rubin (1987); Schenker and Welsh (1988); Tanner and Wong (1987b)) to update parameter estimation. The process stops once the desired number of iterations have been reached. Both the Poor Man's Data Augmentation (PMDA) and the Asymptotic Normal Data Augmentation (ANDA) scheme are implemented (Pan 2000; Wei and Tanner 1991; Tanner and Wong 1987a). When ANDA is chosen, the function calls \code{mvrnorm} (package \code{MASS}).
}
\details{
\tabular{ll}{
Package: \tab MIICD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-02-17\cr
License: \tab GPL-2\cr
}
Multiple Imputation for proportional hazard regression with interval censored data.
}
\author{
Marc Delord

Maintainer: mdelord@gmail.com

}
\references{
PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1, p. 199-203.

Rubin, D. B. (1987). Multiple imputation for nonresponse in surveys. 

Schenker, N. and Welsh, A. (1988). Asymptotic results for multiple imputation. The Annals of Statistics pages 1550-1566.

Tanner, M. A. and Wong, W. H. (1987). An application of imputation to an estimation problem in grouped lifetime analysis. Technometrics 29, 23-32.

Wei, G. C., & Tanner, M. A. (1991). Applications of multiple imputation to the analysis of censored regression data. Biometrics, 47(4), 1297-1309.
}
\keyword{Interval Censored Data}
\keyword{Multiple Imputation}
\keyword{Competing Risks Regression}
\keyword{Poor Man's Data Augmentation}
\keyword{Asymptotic Normal Data Augmentation}


\seealso{
\code{\link{coxph}}
\code{\link{crr}}
\code{\link{FGR}}
\code{\link{mvrnorm}}
\code{\link{cuminc}}
}
\examples{
data(cml)
head(cml)
res <- MIICD.crreg( formula= ~  treatment + SOKAL , data = cml , imax = 10 ,
k = 10 , status = 'status' , cens = 'cens' , trans = 'MMR' ,
method = 'ANDA' , model = 'FG')

res$df

}
