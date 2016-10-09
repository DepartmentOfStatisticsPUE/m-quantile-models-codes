\name{glm.rob.binom}
\alias{glm.rob.binom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
glm.rob.binom
%%  ~~function to do ... ~~
}
\description{
Fitting the robust glm for Binary Data
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
glm.rob.binom(X, y, weights.on.x = FALSE, chuber = 1.345)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Design matrix including the intercept.
%%     ~~Describe \code{X} here~~
}
  \item{y}{The response variable.
%%     ~~Describe \code{y} here~~
}
  \item{weights.on.x}{A character string (can be abbreviated) specifying how points (potential outliers) in x-space are downweighted. Weights on the design of the form sqrt(1-h_(ii)) are used, where h_(ii) are the diagonal elements of the hat matrix.
%%     ~~Describe \code{weights.on.x} here~~
}
  \item{chuber}{Tuning constant for Huber's psi-function.
%%     ~~Describe \code{chuber} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{coef}{Estimated beta coefficients.}
\item{fitted.values}{A vector of fitted values.}
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
