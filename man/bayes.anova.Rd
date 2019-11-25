\name{bayes.anova}
\alias{bayes.anova}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bayes.anova
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function runs a Bayesian analysis of variance (ANOVA) on the data. The Bayesian ANOVA model assumes normally distributed data in all three groups, and conducts inference based on a Gibbs sampler in a three-component Gaussian-mixture with unknown parameters.
}
\usage{
bayes.anova(n=10000,first,second,third,fourth=NULL,fifth=NULL,sixth=NULL,
hyperpars="custom",burnin=n/2,sd="sd",q=0.1,ci=0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
Number of posterior draws the Gibbs sampler produces. The default is \code{n=10000}.
}
  \item{first}{
%%     ~~Describe \code{first} here~~
  Numerical vector containing the values of the first group
}
  \item{second}{
%%     ~~Describe \code{second} here~~
  Numerical vector containing the values of the second group
}
  \item{third}{
%%     ~~Describe \code{third} here~~
  Numerical vector containing the values of the third group
}
  \item{fourth}{
%%     ~~Describe \code{third} here~~
  Numerical vector containing the values of the fourth group. Default value is NULL.
}
  \item{fifth}{
%%     ~~Describe \code{third} here~~
  Numerical vector containing the values of the fifth group. Default value is NULL.
}
  \item{sixth}{
%%     ~~Describe \code{third} here~~
  Numerical vector containing the values of the sixth group. Default value is NULL.
}
  \item{hyperpars}{
%%     ~~Describe \code{hyperpars} here~~
  Sets the hyperparameters on the prior distributions. Two options are provided. The default is \code{"custom"}, and the other is \code{"rafterys"}. For details, see the references.
}
  \item{burnin}{
%%     ~~Describe \code{burnin} here~~
  Burn-in samples for the Gibbs sampler
}
  \item{sd}{
%%     ~~Describe \code{sd} here~~
  Selects if posterior draws should be produced for the standard deviation (default) or the variance. The two options are \code{"sd"} and \code{"var"} respectively.
}
  \item{q}{
%%     ~~Describe \code{q} here~~
  Tuning parameter for the hyperparameters. The default is \code{q=0.1} and it is recommended not to change this.
}
  \item{ci}{
%%     ~~Describe \code{ci} here~~
  The credible level for the credible intervals produced. Default is \code{ci=0.95}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
  The Gibbs sampler is run with four Markov chains to run convergence diagnostics.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
  Returns a dataframe which includes four columns for each parameter of interest. Each column corresponds to the posterior draws of a single Markov chain obtained by the Gibbs sampling algorithm.
}
\references{
%% ~put references to the literature/web site here ~
For details, see: https://arxiv.org/abs/1906.07524v1
}
\author{
%%  ~~who you are~~
Riko Kelter
}
\examples{
set.seed(42)
x1=rnorm(75,0,1)
x2=rnorm(75,1,1)
x3=rnorm(75,2,1)
x4=rnorm(75,-1,1)
result=bayes.anova(first=x1,second=x2,third=x3)
result=bayes.anova(n=1000,first=x1,second=x2,third=x3,
hyperpars="custom",burnin=750,ci=0.99,sd="sd")
result2=bayes.anova(n=1000,first=x1,second=x2,third=x3,
fourth=x4)
}
