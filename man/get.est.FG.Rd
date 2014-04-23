\name{get.est.FG}
\alias{get.est.FG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute estimates.
}
\description{
This function uses the formula interface for Fine-Gray regression competing risk (\code{FGR} from the \code{riskRegression} library) models to compute coeficients and variance estimates. It uses also the \code{cuminc} function from the \code{cmprsk} library to estimate the baseline cause specific hazard.

For internal use.
}
\usage{
get.est.FG( x , data , t1 , status , trans , cens , keep , formula )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{time to events.}
  \item{data}{original dataset.}
  \item{t1}{unique single times at which the survival is to be estimatied.}
   \item{status}{columns indicating the status of observations}
  \item{trans}{transission of interest}
  \item{cens}{censure code}
  \item{keep}{covariate to keep (see crprep for details). Optional here.}
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
}
%\details{}
\value{
\item{surv}{An estimate of the survival function}
\item{sigma}{An estimate of the variance covariance matrix}
\item{beta}{An estimate of the coeficients}
}
\references{
Gerds, TA and Scheike, T and Andersen, PK (2011) Absolute risk regression for competing risks: interpretation, link functions and prediction Research report 11/7. Department of Biostatistics, University of Copenhagen
}
\author{
Marc Delord
}
%\note{}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{FGR}}
\code{\link{riskRegression}}
\code{\link{cuminc}}
}
%\examples{}
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line