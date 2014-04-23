\name{get.set}
\alias{get.set}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute sets of data for PMDA. 
}
\description{
From a set of interval censored data, keep right censored data and drow from a uniform distribution within closed intervals.
}
\usage{
get.set( k , data )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
Number of set to be generated.
}
  \item{data}{Original set from which data sets are drown.}
}
\details{
data must be a data frame representing interval censored data. If one item is censored the column representing the upper bounds of the intervals MUST contain the value 'Inf'.
}
\value{
A data frame with the k so generated data sets. One per column.
}
\references{
PAN, Wei. A Multiple Imputation Approach to Cox Regression with Interval-Censored Data. Biometrics, 2000, vol. 56, no 1, p. 199-203.
}
\author{
Marc Delord
}
%\note{}
\seealso{
\code{\link{runif}}
}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{}
%\keyword{}% __ONLY ONE__ keyword per line
