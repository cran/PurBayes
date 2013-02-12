\name{PB.plot}
\alias{PB.plot}
\title{PurBayes plotting subroutine}
\description{Provides plot functionality for the \code{plot} option in \code{PurBayes}}
\usage{PB.plot(N, Y, out.PB)}
\arguments{
  \item{N}{numeric vector of total reads for each somatic mutation from the tumor tissue NGS data}
  \item{Y}{numeric vector of mutant allele supporting read counts for each somatic mutation from the tumor tissue NGS data}
  \item{out.PB}{output object from \code{PurBayes}}
}
\value{Generates a scatterplot of \emph{Y} against \emph{N}, and superimposes posterior median values of \eqn{\lambda_j} and corresponding symmetric 95\% credible intervals as linear relations.}
\author{Nicholas B. Larson}

\seealso{\code{\link{PurBayes}}}

