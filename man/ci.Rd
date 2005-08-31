% $Id: ci.Rd,v 1.9 2005/06/09 14:20:26 nj7w Exp $
%
\name{ci}
\alias{ci}
\alias{ci.default}
\alias{ci.lm}
\alias{ci.lme}
\title{Compute Confidence Intervals}
\description{ Compute and display confidence intervals for model
  estimates.  Methods for \code{lm} and \code{lme} objects are
  provided. }
\synopsis{
  ci(x, confidence = 0.95, alpha = 1 - confidence,...)
  \method{ci}{default}(x, confidence = 0.95, alpha = 1 - confidence, na.rm=FALSE)...)
  \method{ci}{lm}(x, confidence = 0.95, alpha = 1 - confidence,...)
  \method{ci}{lme}(x, confidence = 0.95, alpha = 1 - confidence,...) 
}
\arguments{
  \item{x}{ object from which to compute confidence intervals. }
  \item{confidence}{ confidence level. Defaults to 0.95. }
  \item{alpha}{type one error rate.  Defaults to 1.0-\code{confidence} }
  \item{na.rm}{boolean indicating whether missing values should be
    removed. Defaults to \code{FALSE}.}
  \item{\dots}{Arguments for methods}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  matrix with one row per model parameter and columns  \code{Estimate},
  \code{CI lower}, \code{CI upper}, \code{Std. Error}, \code{DF} (for
  lme objects only), and \code{p-value}.
}
\author{ Gregory R. Warnes \email{Gregory\_R\_Warnes\@groton.pfizer.com}
}
\seealso{ \code{\ref{MASS}\link{confint}}, \code{\link{lm}},
  \code{\link{summary.lm}} }

\examples{

data(state)
reg  <-  lm(Area ~ Population, data=as.data.frame(state.x77))
ci(reg) 

}
\keyword{ regression }
