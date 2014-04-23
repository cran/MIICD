\name{bcos}
\alias{bcos}
\docType{data}
\title{Breast Cosmesis Study data}
\description{bcos datas}
\usage{data(bcos)}
\format{
  A data frame with 94 observations on the following 3 variables.
  \describe{
    \item{\code{left}}{a numeric vector}
    \item{\code{right}}{a numeric vector}
    \item{\code{treatment}}{a factor with levels \code{Rad} \code{RadChem}}
  }
}
%\details{}
\source{
Finkelstein, D.M., and Wolfe, R.A. (1985). A semiparametric model for regression analysis of interval-censored failure time data. Biometrics 41: 731-740.
}
\references{
FINKELSTEIN, Dianne M. A proportional hazards model for interval-censored failure time data. Biometrics, 1986, p. 845-854.

}
\examples{
data(bcos)
head(bcos)
}
\keyword{Interval Censored Data}
\keyword{datasets}
