\name{write.PB}
\alias{write.PB}
\title{PurBayes model file composer}
\description{Procedurally generates JAGS model files for \code{PurBayes}}
\usage{
write.PB(fn.jags, prior = NULL, het = FALSE)
}
\arguments{
  \item{fn.jags}{File name and location where JAGS model file is written}
  \item{prior}{Prior specification of \eqn{\lambda} under a homogeneous tumor model.  If NULL, defaults to \eqn{Uniform(0,1)}}
  \item{het}{Indicates whether a model file for a heterogeneous tumor model should be written}
}

\value{Writes JAGS model file to location \code{fn.jags}.}

\author{Nicholas B. Larson}
\seealso{\code{\link{PurBayes}}}