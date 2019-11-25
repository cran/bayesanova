\name{post.fit}
\alias{post.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
post.fit
}
\description{
Provides a visualization of the posterior fit produces by a Bayesian ANOVA.
}
\usage{
post.fit(dataframe,x1,x2,x3,x4=NULL,x5=NULL,x6=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{
  A dataframe returned by \code{bayes.anova}
}
  \item{x1}{
  A numerical vector containing the values of the first group
}
  \item{x2}{
  A numerical vector containing the values of the second group
}
  \item{x3}{
  A numerical vector containing the values of the third group
}
  \item{x4}{
  A numerical vector containing the values of the fourth group. Default value is NULL.
}
  \item{x5}{
  A numerical vector containing the values of the fifth group. Default value is NULL.
}
  \item{x6}{
  A numerical vector containing the values of the sixth group. Default value is NULL.
}
}
\details{
  Provides a histogram and density overlays of each Gaussian mixture component obtained by the \code{bayes.anova} function. This way, the mixture density fit of the Gaussian mixture can be compared to the original data.
}
\value{
  Produces a plot consisting of a histogram of the original data and the Gaussian mixture components as density overlays.
}
\author{
  Riko Kelter
}
\examples{
set.seed(700)
x1=rnorm(1000,0,1)
x2=rnorm(1000,1,1)
x3=rnorm(1000,2,1)

result=bayes.anova(n=1000,first = x1, second=x2, third=x3)
post.fit(result,x1=x1,x2=x2,x3=x3)
}