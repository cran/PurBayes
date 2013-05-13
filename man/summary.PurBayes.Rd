\name{summary.PurBayes}
\alias{summary.PurBayes}
\title{Analysis summary from PurBayes}
\description{Provides summary analysis of model fitting from \code{PurBayes}}
\usage{
  \method{summary}{PurBayes}(object,\dots)
  }
\arguments{
  \item{object}{Output from \code{PurBayes}}
  \item{\dots}{Other arguments passed to or from other methods}
  }
\value{Returns a list with object class \code{summary.PurBayes}:
  \item{purity}{Posterior distribution summary of tumor purity}
  \item{post.dist}{Posterior distribution summary of all model parameters}
  \item{n.pop}{Number of mutation populations detected by PurBayes}
  \item{dev.out}{Penalized expected deviance results from model selection procedure}
  }

\author{Nicholas B. Larson}
\seealso{\code{\link{PurBayes}},\code{\link{print.summary.PurBayes}}}