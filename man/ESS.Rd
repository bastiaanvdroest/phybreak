% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ESS.R
\name{ESS}
\alias{ESS}
\alias{ESS.phybreak}
\alias{ESS.factor}
\alias{ESS.numeric}
\title{Effective sample size}
\usage{
ESS(x)

\method{ESS}{phybreak}(x)

\method{ESS}{factor}(x)

\method{ESS}{numeric}(x)
}
\arguments{
\item{x}{Vector of class \code{numeric}, containing only integers. This will be treated as a categorical variable.}
}
\value{
Effective sample sizes.
}
\description{
Effective sample size
}
\details{
When applied to an object of class \code{phybreak}, \code{ESS} calculates for all parameters and continuous
  variables (infection times) the effective sample size (ESS) with the \code{\link[coda]{effectiveSize}} 
  function in \pkg{coda}. For the infectors, 
  a method is used that is similar to the method for the approximate ESS for phylogenetic trees, described in 
  Lanfaer et al (2016):
  \enumerate{
    \item Define as distance measure between two sampled infectors \eqn{D(i,j) = 0} if \eqn{i = j}, 
        and \eqn{D(i,j) = 1} if \eqn{i \neq= j}{i <> j}
    \item Calculate the mean squared distance f(k) between sampled infectors at intervals k = 1,2,... in the mcmc chain.
        The distance will increase with increasing interval k.
    \item Use the rate at which f(k) approaches the asymptote to calculate the ESS (see Lanfaer et al, 2016)
  }
  The latter method can also be directly called for single vectors of class \code{factor} or \code{integer}.
}
\section{Methods (by class)}{
\itemize{
\item \code{ESS(phybreak)}: Effective sample size of phybreak posterior.

\item \code{ESS(factor)}: Effective sample size of a categorical variable.

\item \code{ESS(numeric)}: Effective sample size of a categorical variable.

}}
\examples{
#First create a phybreak object
simulation <- sim_phybreak(obsize = 5)
MCMCstate <- phybreak(dataset = simulation)

MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
ESS(MCMCstate)
}
\references{
\href{http://dx.doi.org/10.1093/gbe/evw171}{Lanfaer et al. (2016)} Estimating 
  the effective sample size of tree topologies from Bayesian phylogenetic analyses. 
  \emph{Genome Biol Evol}, \strong{8}(8): 2319-2332.
}
\author{
Don Klinkenberg \email{don@xs4all.nl}
}
