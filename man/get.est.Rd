\name{get.est}
\alias{get.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Cox and Breslow estimates.
}
\description{
For internal uses.
}
\usage{
get.est( x , data , t1 , formula )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{time to events.}
  \item{data}{original dataset.}
  \item{t1}{unique single times at which the survival is to be estimatied.}
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
}
%\details{}
\value{
\item{surv}{An estimate of the survival function}
\item{sigma}{An estimate of the variance covariance matrix}
\item{beta}{An estimate of the coeficients}
}
\references{
PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1, p. 199-203.
}
\author{Marc Delord <mdelord@gmail.com>}
%\note{}
\seealso{
\code{\link{coxph}}
}
%\examples{}