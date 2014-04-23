\name{MIICD.crreg}
\alias{MIICD.crreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple Imputation for competing risks regression}
\description{This function uses the multiple imputation approach of Pan (2000) to handle interval censored data with competing outcomes using the Fine and Gray Propotional Hazards model (or Cox Proportional Hazard Model). It calls iteratively the \code{crr} function on imputed datasets and derives multiple estimates from imputed data sets. Finally it combines multiple estimates following multiple imputation rules (Rubin 1987; Schenker and Welsh 1988; Tanner and Wong 1987b) to update parameter estimation. The process stops once the desired number of iteration have been reached. Both the Poor Man's Data Augmentation (PMDA) as well as the Asymptotic Normal Data Augmentation (ANDA) scheme are implemented (Pan 2000; Wei and Tanner 1991; Tanner and Wong 1987a). 
When ANDA is chosen, the function calls \code{mvrnorm} from package \code{MASS}.
}
\usage{
MIICD.crreg(formula, data, imax = 25, k = 10, th0 = 1e-3, status, trans,
cens, keep, method = c("PMDA", "ANDA"), model = c("FG", "Cox"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula}
  \item{data}{a data frame to look for interval censored data and covariates}
  \item{imax}{maximum number of iteration}
  \item{k}{number of dataset generated at each iteration}
  \item{th0}{note used currently}
  \item{status}{collumn name were the status of observation is indicated}
  \item{trans}{what is the transition of interest in the status comlumn}
  \item{cens}{code used for censored data in the status column}
  \item{keep}{must be the same as the covariate used in formula}
  \item{method}{which imputation scheme shall be used PMDA or ANDA}
  \item{model}{which model shall be used currently FG state for Fine and Gray proportional hazards regression model}
}
\details{
The data frame MUST have one column named 'left' one column 'right'. Interval censored data are observation for which 'left' < 'right'. Inf in the 'right' column stands for right sensored observations. Competing event must have exact failure time, that is 'left' = 'right'.
}
\value{
\item{Mean beta}{estimation of coefficients computed using the last iterations}
\item{Mean sigma}{estimation of sigma computed using the last iterations}
\item{betas}{list of coefficient estimates for each iteration}
\item{call}{function call}
\item{df}{results returned in a data frame}
\item{niter}{number of iteration}
\item{conv}{successive mean of coefficient estimates over iterations}
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
\seealso{
\code{\link{MIICD.coxph}}
\code{\link{crr}}
\code{\link{mvrnorm}}
}
\examples{
data(cml)
head(cml)
#res<-MIICD.crreg( formula= ~  treatment + SOKAL  , data = cml , imax = 10 ,
#k = 10 , th0 = 1e-1 ,  status = 'status' , cens = 'cens' , trans = 'MMR' ,
#method = 'ANDA' , model = 'FG' , keep=c('treatment','SOKAL')
#res$df

}
\keyword{Interval Censored Data}
\keyword{Multiple Imputation}
\keyword{Asymptotic Normal Data Augmentation}
\keyword{Competing Risks Regression}