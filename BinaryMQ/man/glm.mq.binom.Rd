\name{glm.mq.binom}
\alias{glm.mq.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
glm.mq.binom
%%  ~~function to do ... ~~
}
\description{
Fitting the M-quantiles to Binary Data
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
glm.mq.binom(x, y, case.weights = rep(1, nrow(x)), maxit = 100, acc = 1e-04, var.weights = rep(1, nrow(x)), weights.x = FALSE, q = 0.5, k = 1.345)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix including the intercept.
%%     ~~Describe \code{x} here~~
}
  \item{y}{The response variable.
%%     ~~Describe \code{y} here~~
}
  \item{case.weights}{An optional vector of weights to be used in the fitting process.
%%     ~~Describe \code{case.weights} here~~
}
  \item{maxit}{An integer giving the number of maximum iterations (default maxit=100)
%%     ~~Describe \code{maxit} here~~
}
  \item{acc}{Positive convergence tolerance (default acc=1e-04).
%%     ~~Describe \code{acc} here~~
}
  \item{var.weights}{An optional vector of weights of the variance to be used in the fitting process.
%%     ~~Describe \code{var.weights} here~~
}
  \item{weights.x}{A character string (can be abbreviated) specifying how points (potential outliers) in x-space are downweighted. Weights on the design of the form sqrt(1-h_(ii)) are used, where h_(ii) are the diagonal elements of the hat matrix.
%%     ~~Describe \code{weights.x} here~~
}
  \item{q}{A vector of M-quantiles.
%%     ~~Describe \code{q} here~~
}
  \item{k}{Tuning constant k for Huber's psi-function (default is 1.345).
%%     ~~Describe \code{k} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{coefficients}{Estimated beta coefficients.}
\item{var.beta}{Estimated variance of the beta coefficients.}
\item{fitted.values}{A vector of fitted values for each M-quantile.}
\item{residuals}{A vector of residuals for each M-quantile.}
\item{q.values}{The grid of M-quantiles.}
\item{matQ, matM}{The matrices for computing the variance of the beta coefficients.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
