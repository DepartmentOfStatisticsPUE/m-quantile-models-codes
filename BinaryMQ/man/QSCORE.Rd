\name{QSCORE}
\alias{QSCORE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QSCORE
%%  ~~function to do ... ~~
}
\description{A function to compute the M-quantile coefficients
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
QSCORE(y, yhatq, qvals)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The response variable.
%%     ~~Describe \code{y} here~~
}
  \item{yhatq}{The fitted values by Binary MQ.
%%     ~~Describe \code{yhatq} here~~
}
  \item{qvals}{A grid of M-quantile values used in QLogit.Scores.
%%     ~~Describe \code{qvals} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{qscore}{M-quantile coefficients for each unit.}
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
