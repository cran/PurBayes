\name{print.summary.PurBayes}
\alias{print.summary.PurBayes}
\title{Print Summmary of PurBayes Analysis}
\description{S3 method for printing summary of PurBayes output}
\usage{
  \method{print}{summary.PurBayes}(x,\dots)
  }
\arguments{
  \item{x}{Object of class \code{summary.PurBayes}}
  \item{\dots}{Other arguments passed to or from other methods}
  }
  
\details{Prints out summary information from the PurBayes model fit}
\value{Prints the following information:
  \item{Purity Estimate}{Posterior distribution summary of tumor purity, given as the posteior median and 2.5\% and 97.5\% quantiles.}
  \item{Number of Populations}{Value of \eqn{J}, number of discrete populations of mutations detected by PurBayes}
  \item{Model Posterior Distributions}{Quantiles of all of the parameter posterior distributions in the final model}
  \item{Penalized Deviance Results}{Matrix of penalized deviance results used in determining the final model fit}
}

\author{Nicholas B. Larson}
\seealso{\code{\link{PurBayes}},\code{\link{summary.PurBayes}}}