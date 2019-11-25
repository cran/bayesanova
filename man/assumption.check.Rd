\name{assumption.check}
\alias{assumption.check}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
assumption.check
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function checks the assumption of normality for each of the groups \code{x1}, \code{x2}, \code{x3} (and optional \code{x4}, \code{x5} and \code{6}) used in the Bayesian ANOVA via Shapiro-Wilk tests with confidence level \code{conf.level}.
}
\usage{
assumption.check(x1,x2,x3,x4=NULL,x5=NULL,x6=NULL,conf.level=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
%%     ~~Describe \code{x1} here~~
  Numerical vector containing the values for the first group
}
  \item{x2}{
%%     ~~Describe \code{x2} here~~
  Numerical vector containing the values for the second group
}
  \item{x3}{
%%     ~~Describe \code{x3} here~~
  Numerical vector containing the values for the third group
}
\item{x4}{
%%     ~~Describe \code{x3} here~~
  Numerical vector containing the values for the fourth group. Defaults to NULL.
}
\item{x5}{
%%     ~~Describe \code{x3} here~~
  Numerical vector containing the values for the fifth group. Defaults to NULL.
}
\item{x6}{
%%     ~~Describe \code{x3} here~~
  Numerical vector containing the values for the sixth group. Defaults to NULL.
}
  \item{conf.level}{
%%     ~~Describe \code{conf.level} here~~
  Confidence level of the Shapiro-Wilk test used. Significance level equals \code{1-conf.level}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
  If a single Shapiro-Wilk test fails, the method returns a warning and recommends to use further diagnostics.
}
\value{
  Histograms and Quantile-Quantile plots for all groups are produced, either a warning or a confirmation of normality in all three groups is printed to the console.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
%%  ~~who you are~~
  Riko Kelter
}
\references{
For details, see: https://arxiv.org/abs/1906.07524v1
}
\examples{
set.seed(42)
x1=rnorm(75,0,1)
x2=rnorm(75,1,1)
x3=rnorm(75,2,1)

assumption.check(x1,x2,x3,conf.level = 0.95)
}
