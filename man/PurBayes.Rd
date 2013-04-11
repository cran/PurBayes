\name{PurBayes}
\alias{PurBayes}
\title{Bayesian Estimation of Tumor Purity and Clonality}
\description{PurBayes is an iterative Bayesian algorithm which simultaneously estimates tumor purity and clonality using finite mixture models, using the MCMC software JAGS to obtain posterior samples for inference.  Using a penalized deviance criterion, PurBayes iteratively fits models increasing in variant population count until an optimal fit is achieved.}
\usage{
PurBayes(N, Y, M = NULL, Z = NULL, pop.max = 5, prior = NULL, burn.in = 50000, n.post = 10000,
 fn.jags = "PB.jags", plot = FALSE)
}
\arguments{
  \item{N}{numeric vector of total reads for each somatic mutation from the tumor tissue NGS data}
  \item{Y}{numeric vector of mutant allele supporting read counts for each somatic mutation from the tumor tissue NGS data}
  \item{M}{optional numeric vector of total reads for germline heterogyous variants.  PurBayes uses these to estimate non-reference allele mapping rate to account for mapping bias}
  \item{Z}{optional numeric vector of alternate allele reads for germline heterozygous variants, corresponding to \code{M}}
  \item{pop.max}{Maximum number of variant populations allowed in the iterative modeling procedure.  Defaults to 5.}
  \item{prior}{Optional prior distribution for \eqn{\lambda_J} under the homogenenous tumor model.  If \code{NULL}, defaults to Uniform(0,1).  WARNING: This must be provided as a character string written within the JAGS modeling language.}
  \item{burn.in}{Number of MCMC draws that are excluded as a burn-in.  Defaults to 50000.}
  \item{n.post}{Number of MCMC draws that are sampled for posterior inference.  Defaults to 10000.}
  \item{fn.jags}{File location and name to which \code{write.PB} generates the appropriate JAGS model file.  Defaults to 'PB.jags' in the current working directory.}
  \item{plot}{If \code{plot=TRUE}, then \code{PB.plot} is called to generate a visual representation of the data along with the model fit by \code{PurBayes}.  Defaults to \code{FALSE}.}
}
\details{For a given tumor purity level \eqn{\lambda} PurBayes assumes a binomial-binomial mixture model for the tumor sequence reads which support the alternate allele, \eqn{Y_i^t \sim Bin(N_i,\lambda/2)}.  This model is fit to the data under the assumption of tumor homogeneity.  PurBayes also supports the possibility of intra-tumor heterogeneity, whereby the tumor tissue is comprised of additional subclonal variant populations, each with its own 'purity', \eqn{\lambda_j<\lambda}, for \eqn{j = 1,...,J-1} and \eqn{\lambda_J \equiv \lambda}.

The probability that a given variant corresponds to the \eqn{j^{th}} population is given by \eqn{\kappa_j}, and \eqn{\bm{\kappa}=(\kappa_1,\ldots,\kappa_J)} follows a dirichlet prior such that \eqn{\pi(\bm{\kappa})\sim Dirichlet(\alpha_1,\,\ldots,\alpha_J)} for a given variant population quantity \eqn{J}.  PurBayes applies a diffuse prior on \eqn{\bm{\kappa}}, such that \eqn{\alpha_1=\ldots=\alpha_J=1}.  While the user may specify a particular prior for \eqn{\lambda} under a homogeneous tumor, PurBayes defaults to \eqn{\pi(\lambda_j) \sim Uniform(0,1)} for all \emph{j}, and uses a sort function to avoid label switching.

The optimality criterion used for model selection with regard to size of \eqn{J} is based upon the penalized expected deviance (Plummer, 2008).
}
\value{
 \item{n.pop}{Numeric scalar corresponding to number of variant populations detected by PurBayes}
 \item{PB.post}{\code{mcmc.list} object corresponding to posterior samples of PurBayes model parameters.  This necessarily includes \emph{pur}, the tumor purity.  If \code{n.pop>1}, posterior samples of \eqn{\kappa_j} and \eqn{\lambda_j} for \eqn{j = 1,...,J} are also included.}
}
\references{
Plummer, M. (2008) Penalized loss functions for Bayesian model comparison. \emph{Biostatistics} doi:
10.1093/biostatistics/kxm049
}
\author{Nicholas B. Larson}
\examples{
#Homogeneous tumor example
N.var<-20
N<-round(runif(N.var,20,200))
lambda<-0.75
Y<-rbinom(N.var,N,lambda/2)
\dontrun{PB.hom<-PurBayes(N,Y)}

#Heterogeneous tumor example - 1 subclonal population
N.var<-20
N<-round(runif(N.var,20,200))
lambda.1<-0.75
lambda.2<-0.25
lambda<-c(rep(lambda.1,10),rep(lambda.2,10))
Y<-rbinom(N.var,N,lambda/2)
\dontrun{PB.het<-PurBayes(N,Y)}
}

