\name{anovaplot}
\alias{anovaplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
anovaplot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plots the results of a Bayesian ANOVA
}
\usage{
anovaplot(dataframe, type="rope", sd="sd", ci=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataframe}{
%%     ~~Describe \code{dataframe} here~~
  A dataframe which is the result of a Bayesian ANOVA
}
  \item{type}{
%%     ~~Describe \code{type} here~~
  Selects the type of plot which should be produced. The default is \code{"rope"}, which produces a posterior region of practical equivalence (ROPE) analysis and posterior distributions for the three effect sizes of interest. Other options are to produce posteriors for each parameter via \code{type=="pars"}, or to produce posteriors for the difference of means and variances via \code{type=="diff"}.
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
  Selects if the results include posterior draws for the standard deviation (default) or the variance. \code{sd="sd"} is the default. \code{sd="var"} assumes that the dataframe includes posterior draws of the variances of each group.
}
  \item{ci}{
%%     ~~Describe \code{ci} here~~
The credible level used for producing credible intervals. The default is \code{ci=0.95}.
}
}
\value{
  Produces plots of the results depending on which type is selected.
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
result=bayes.anova(n=1000,first=x1,second=x2,third=x3)

anovaplot(result)
anovaplot(result, type="effect")

x4=rnorm(75,3,1)
result2=bayes.anova(n=1000,first=x1,second=x2,third=x3,fourth=x4)
anovaplot(result2)
}
