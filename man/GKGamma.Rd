\name{GKGamma}
\alias{GKGamma}
\title{Calculate Goodman-Kruskal's Gamma}
\description{
  Calculate Goodman-Kruskal's Gamma coefficient.  Code taken from \code{GoodmanKruskalGamma} in \pkg{DescTools} v0.99.22.
}
\usage{
  GKGamma(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric vector or a contingency table.  A matrix will be treated as a table.}

  \item{y}{NULL(default) or a vector with compatible dimensions to \code{x}.  If \code{y} is provided, \code{table(x,y)} is calculated.}

  \item{conf.level}{confidence level of the interval.  If set to \code{NA} (the default), no confidence intervals will be calculated.}

  \item{...}{Other arguments to be passed to \code{table}.}
}
\details{
See help file for \code{\link[DescTools]{GoodmanKruskalGamma}} for more details.

}
\value{A scalar value for the test statistic and optionally a confidence interval.}

\author{
  Dave Armstrong <dave@quantoid.net>
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\keyword{summaryStats}
