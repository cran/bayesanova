\name{post.pred.check}
\alias{post.pred.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
post.pred.check
}
\description{
Provides a posterior predictive check for a fitted Bayesian ANOVA model.
}
\usage{
post.pred.check(anovafit, ngroups, out, reps = 50, eta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{anovafit}{
  A dataframe returned by \code{bayes.anova}
}
  \item{ngroups}{
  An integer which is the number of groups used in the ANOVA
}
  \item{out}{
  A numerical vector containing the originally observed data in all groups
}
  \item{reps}{
  An integer which is the number of posterior predictive distributions sampled from the ANOVA models posterior distribution. Defaults to 50 sampled parameters.
}
  \item{eta}{
  A numerical vector containing the weight values of the mixture.
}
}
\details{
  Provides a posterior predictive check for a fitted Bayesian ANOVA model.
}
\value{
  Produces a plot consisting of a density estimate of the original data and posterior predictive distributions sampled from the posterior of the Bayesian ANOVA model as density overlays.
}
\author{
  Riko Kelter
}
\examples{
set.seed(700)
x1=rnorm(1000,0,1)
x2=rnorm(1000,1,1)
x3=rnorm(1000,2,2)

result=bayes.anova(n=1000,first = x1, second=x2, third=x3)
post.pred.check(result, ngroups = 3, out = c(x1,x2,x3), reps = 25, eta = c(1/3,1/3,1/3))
}