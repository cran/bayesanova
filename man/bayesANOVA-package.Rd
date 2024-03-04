\name{bayesanova-package}
\alias{bayesanova-package}
\docType{package}
\title{
Bayesian ANOVA
}
\description{
Provides a Bayesian version of the analysis of variance (ANOVA) based on a three-component Gaussian mixture, for which a Gibbs sampler produces the posteriors of the means and standard deviation of each component. Also, model assumptions can be checked and results visualised.
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{bayesANOVA}
\packageIndices{bayesANOVA}
The core function is \code{bayes.anova} which provides the Bayesian version of the ANOVA. Also, assumptions can be checked via \code{assumption.check}, and \code{anovaplot} produces visualizations of the results.
}
\author{
Riko Kelter

Maintainer: Riko Kelter
}
\references{
For details, see: https://arxiv.org/abs/1906.07524v1
}
\keyword{ package }
\examples{
set.seed(42)
x1=rnorm(75,0,1)
x2=rnorm(75,1,1)
x3=rnorm(75,2,1)
assumption.check(x1,x2,x3,conf.level = 0.95)
result=bayes.anova(n=1000,first=x1,second=x2,third=x3)
anovaplot(result)
}
