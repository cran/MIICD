\name{preproc.coxph}
\alias{preproc.coxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{preprocessed the data for \code{MIICD.coxph}}
\description{For internal uses}
\usage{preproc.coxph(data = data, k = k)}
\arguments{
  \item{data}{original dataset}
 \item{k}{number of dataset to be imputeed at each iteration}
}

\details{data1 is a data frame were the same value is repeated k time for each row. This function has been made to reduce the size of realy interval sensored data what will be imputed iteratively. The idea is to treat right censored and exact time failure observation ones for all.

}
\value{
\item{data2}{subset of interval censored data}
\item{t1}{single times at wich the survival function is estimated}
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
\code{\link{preproc.crreg}}
}
%\examples{}
