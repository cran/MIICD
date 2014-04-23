\name{get.est.cr.cox}
\alias{get.est.cr.cox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute estimates.}
\description{
This function uses the \emph{crprep} function from \code{mstate} package to compute coeficients, variance and Breslow baseline estimate of hazards using a cox proportionnal hazard model with counting process notation for competing risks.
For internal use.
}
\usage{get.est.cr.cox( x , data , t1 , status , trans , cens , keep , formula )}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{time to events.}
  \item{data}{original dataset.}
  \item{t1}{unique single times at which the survival is to be estimatied.}
  \item{status}{columns indicating the status of observations}
  \item{trans}{transission of interest}
  \item{cens}{censure code}
  \item{keep}{covariate to keep (see crprep for details)}
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
}
%\details{}
\value{
\item{surv}{An estimate of the survival function}
\item{sigma}{An estimate of the variance covariance matrix}
\item{beta}{An estimate of the coeficients}
}
\references{
GESKUS, Ronald B. Cause-Specific Cumulative Incidence Estimation and the Fine and Gray Model Under Both Left Truncation and Right Censoring. Biometrics, 2011, vol. 67, no 1, p. 39-49.
}
\author{
Marc Delord <mdelord@gmail.com>
}
%\note{}
\seealso{
\code{\link{crprep}}
\code{\link{mstate}}
}
%\examples{}