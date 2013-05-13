\name{plot.PurBayes}
\alias{plot.PurBayes}
\title{Plotting Function for PurBayes}
\description{Generates plot of somatic mutation data and model fit}
\usage{
  \method{plot}{PurBayes}(x,\dots)
  }
\arguments{
  \item{x}{Output object from \code{PurBayes}}
  \item{\dots}{Other arguments passed to or from other methods}
}
\value{Produces a graphical representation of the somatic mutation data and estimated purity (and subclonal cellularities, if detected).  These are presented as linear relations against the read data, along with 95\% credible intervals.}
\author{Nicholas B. Larson}

\seealso{\code{\link{PurBayes}}}