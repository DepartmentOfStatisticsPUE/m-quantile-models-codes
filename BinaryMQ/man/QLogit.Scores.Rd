\name{QLogit.Scores}
\alias{QLogit.Scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QLogit.Scores
%%  ~~function to do ... ~~
}
\description{Computing the M-quantile coefficient at unit level
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QLogit.Scores(x, y, case.weights = rep(1, nrow(x)), maxit = 50, acc = 1e-04, weights.x = FALSE, qgrid = c(0.1, 0.2, 0.25, 0.3, 0.4, seq(from = 0.45, to = 0.55, by = 0.005), 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9), plots = F, k.value = 1.345)
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
  \item{weights.x}{A character string (can be abbreviated) specifying how points (potential outliers) in x-space are downweighted. Weights on the design of the form sqrt(1-h_(ii)) are used, where h_(ii) are the diagonal elements of the hat matrix.
%%     ~~Describe \code{weights.x} here~~
}
  \item{qgrid}{A grid of M-quantile values (default grid=c(0.1,0.2,0.25,0.30,0.4,seq(from=0.45,to=0.55,by=0.005),0.60,0.65,0.70,0.75,0.8,0.85,0.9)
%%     ~~Describe \code{qgrid} here~~
}
  \item{plots}{Plot on QLogistic Coefficients (default is FALSE)
%%     ~~Describe \code{plots} here~~
}
  \item{k.value}{Tuning constant for Huber's psi-function.
%%     ~~Describe \code{k.value} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{coefficients}{Estimated beta coefficients.}
\item{fitted.values}{A vector of fitted values for each M-quantile.}
\item{residuals}{A vector of residuals for each M-quantile.}
\item{qscores}{M-quantile coefficients at unit level.}
\item{qscores.std}{Standardized M-quantile coefficients at unit level.}
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
