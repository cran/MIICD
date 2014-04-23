\name{MIICD.coxph}
\alias{MIICD.coxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple Imputation for Cox proportional hazards model}
\description{
This function uses the multiple imputation approach of Pan (2000) to handle interval censored data with the Cox Propotional Hazards model. It calls iteratively the \code{coxph} function on imputed datasets and derives multiple estimates from imputed data sets. Finally it combines multiple estimates following multiple some imputation rules (Rubin 1987; Schenker and Welsh 1988; Tanner and Wong 1987b) to update parameter estimation. The process stops once the desired number of iterations have been reached. Both the Poor Man's Data Augmentation (PMDA) as well as the Asymptotic Normal Data Augmentation (ANDA) scheme are implemented (Pan 2000; Wei and Tanner 1991; Tanner and Wong 1987a). 
}
\usage{
MIICD.coxph(formula, data, imax = 25, k = 10, th0 = 1e-03, method = c("PMDA", "ANDA"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{formula}{a formula}
  \item{data}{a data frame to look for interval censored data and covariates}
  \item{imax}{maximum number of iteration}
  \item{k}{number of dataset generated at each iteration}
  \item{th0}{note used currently}
  \item{method}{choose between the Poor Man's Data Augmentation (PMDA) or Asymptotic Normal Data Augmentation (ANDA) scheme}
}

\details{
The data frame MUST have one column named 'left' one column 'right' and one column named 'status'. Interval censored data are observation for which 'left' < 'right'. Inf in the 'right' column stands for right sensored observations.
}
\value{
\item{Mean beta}{estimation of coefficients computed using the last iterations}
\item{Mean sigma}{estimation of sigma computed using the last iterations}
\item{betas}{list of coefficient estimates for each iteration}
\item{call}{function call}
\item{df}{results returned in a data frame}
\item{niter}{number of iteration}
\item{conv}{successive mean of coefficient estimates over iterations}
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
\code{\link{MIICD.crreg}}
\code{\link{coxph}}
}
\examples{
data(bcos)
head(bcos)
#res <- MIICD.coxph( formula = ~ treatment , data = bcos ,
#imax = 10 , k = 10 , th0 = 0 , method = 'ANDA' )
#res$df
}
  % Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Interval Censored Data}
\keyword{Multiple Imputation}
\keyword{Asymptotic Normal Data Augmentation}
\keyword{Poor Man's Data Augmentation}
\keyword{Cox Proportional Hazards Regression}
