\name{plotCIgroup}
\alias{plotCIgroup}
\title{Plot Confidence Intervals by Group}
\description{
  Plot confidence intervals of a quantitative variable by a grouping variable
}
\usage{
  plotCIgroup(form, data, horiz=FALSE,
    arrowlen = 0, includeOverall=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{form}{A formula with quantitative variable on the left-hand side and grouping variable on the right-hand side}
  \item{data}{A data frame from which the variables in the formula can be obtained}
  \item{horiz}{Logical flag indicating whether the bars should be vertical (default) or horizontal, if \code{TRUE}}
  \item{arrowlen}{Lenth of the arrowheads at the end of the confidence intervals}
  \item{includeOverall}{Logical flag for whether or not to include the overall (un-grouped) confidence interval, \code{TRUE} by default}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{A plot}

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
