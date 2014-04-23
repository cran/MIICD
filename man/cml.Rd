\name{cml}
\alias{cml}
\docType{data}
\title{CML data}
\description{
Time to major molecular response
}
\usage{data(cml)}
\format{
  A data frame with 236 observations on the following 4 variables.
  \describe{
    \item{\code{left}}{a numeric vector, stands for the left hand side of intervals}
    \item{\code{right}}{a numeric vector, stands for the right hand side of intervals. If the value is \code{Inf}, the corresponding data is right censored}
    \item{\code{treatment}}{factor \code{400mg} \code{600mg} , indicates the daily dose of imatinib}
    \item{\code{status}}{a factor with levels \code{cens} \code{Echec} \code{MMR} \code{Tox}. \code{cens} is the censor indicator. \code{Echec} stands for a treatment failure (exact time failure), \code{Tox} stands for toxicity appearence (exact time failure). \code{MMR} stands for Major Molecular Response : the event of main interest which is interval censored.}
  \item{\code{SOKAL} }{factor : prognostic score for CML patients}
  \item{\code{Sexe} }{factor : sex of patients}
  }
}
%\details{}
\source{
Preudhomme, C., Guilhot, J., Nicolini, F. E., Guerci-Bresler, A., Rigal-Huguet, F., Maloisel, F., ... & Guilhot, F. (2010). Imatinib plus peginterferon alfa-2a in chronic myeloid leukemia. New England Journal of Medicine, 363(26), 2511-2521.
}
\references{DELORD, Marc, ROUSSELOT, Philippe, CAYUELA, Jean Michel, et al. High imatinib dose overcomes insufficient response associated with ABCG2 haplotype in chronic myelogenous leukemia patients. Oncotarget, 2013, vol. 4, no 10, p. 1582.}
\examples{
data(cml)
head(cml)}
\keyword{Interval Censored Data}
\keyword{datasets}