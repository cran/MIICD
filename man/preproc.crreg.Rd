\name{preproc.crreg}
\alias{preproc.crreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{preproces the data for \code{MIICD.crreg}}

\description{For internal uses}

\usage{
preproc.crreg(data = data, k = k, trans = trans, status = status, lt = 3 * nrow(data))
}
\arguments{
  \item{data}{entry data}
  \item{k}{number of dataset to be imputeed at each iteration}
  \item{trans}{transition of interest}
  \item{status}{name of the column of status}
  \item{lt}{number of time points were the cumulative incidence function is to be estimated}
}
\details{data1 is a data frame were the same value is repeated k time for each row. This function has been made to reduce the size of realy interval sensored data what will be imputed iteratively. The idea is to treat right censored and exact time failure observation ones for all.
}

\value{
\item{data2}{subset of interval censored data}
\item{t1}{single times at wich the cumulative incidence function is estimated}
\item{data1}{data frame of imputated data for exact time observation or right censored data}
\item{or}{original order of the data}
\item{I}{indicator of interval censored observations}
}
%\references{}
\author{
Marc Delord <mdelord@gmail.com>
}
%\note{}
\seealso{
\code{\link{preproc.coxph}}
}
%\examples{}
