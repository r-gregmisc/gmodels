% $Id: estimable.Rd,v 1.14 2005/06/09 14:20:26 nj7w Exp $
%
\name{estimable}
\alias{estimable}
\alias{.wald}
\title{Contrasts and estimable linear functions of model coefficients
  for lm, glm, lme, and geese objects}
\description{Compute and test contrasts and other estimable linear
  functions of model coefficients for for lm, glm, lme, and geese
  objects}
\usage{
estimable(obj, cm, beta0, conf.int=NULL, joint.test=FALSE, show.beta0)
.wald(obj, cm,beta0=rep(0, ifelse(is.null(nrow(cm)), 1, nrow(cm))))
}
\arguments{
   \item{obj}{Regression (lm,glm,lme) object. }
   \item{cm}{Matrix specifying estimable linear functions or contrasts (one per
     row). The number of columns must match the number of fitted
     coefficients in the model.}
   \item{beta0}{Vector of null hypothesis values}
   \item{conf.int}{Confidence level.  If provided, confidence intervals
     will be computed.}
   \item{joint.test}{Logical value. If TRUE a 'joint' Wald test for the
     hypothesis \eqn{L \beta=\beta_0}{L \%*\% beta=beta0} is performed.
     Otherwise 'row-wise' tests are performed, i.e. \eqn{(L
      \beta)_i=\beta_{0i}}{(L \%*\% beta)[i]=beta0[i]}
    }
   \item{show.beta0}{Logical value. If TRUE a column for beta0 will be
     included in the output table.  Defaults to TRUE when beta0 is
     specified, FALSE otherwise.}
}
\details{
  \code{estimable} computes an estimate, test statitic, significance
  test, and (optional) confidence interval for each linear functions of
  the model coefficients specified by the rows of \code{cm}.  The
  estimates and their variances are obtained by applying the matrix
  \code{cm} to the model estimates variance-covariance matrix.  Degrees
  of freedom are obtained from the appropriate model terms.
  
  The user is responsible for ensuring that the specified
  linear functions are meaningful.

  For computing contrasts among levels of a single factor,
  \code{\link{fit.contrast}} may be more convenient.  For computing
  contrasts between two specific combinations of model parameters, the
  \code{contrast} function in Frank Harrell's Design library may be more
  convenient.

  The \code{.wald} function is called internally by \code{estimable} and
  is not intended for direct use.
}
\note{
  The estimated fixed effect parameter of \code{lme} objects may have
  different degrees of freedom.  If a specified contrast includes
  nonzero coefficients for parameters with differing degrees of freedom,
  the smallest number of degrees of freedom is used and a warning
  message is issued.
  }
\value{
  Returns a matrix with one row per linear function.  Columns contain
  the beta0 value (optional, see \code{show.beta0} above), estimated
  coefficients, standard errors, t values, degrees of freedom, two-sided
  p-values, and the lower and upper endpoints of the
  1-alpha confidence intervals. 
}
\author{ BXC (Bendix Carstensen) \email{bxc\@novonordisk.com}, 
  Gregory R. Warnes \email{Gregory\_R\_Warnes\@groton.pfizer.com}, and
  Søren Højsgaard \email{sorenh@agrsci.dk}}
\seealso{
  \code{\link{fit.contrast}},
  \code{\link[stats]{lm}}, \code{\link[nlme]{lme}},
  \code{\link[stats]{contrasts}},
  \code{\link[Design]{contrast}},
  }

\examples{
# simple contrast and confidence interval
y <- rnorm(100)
x <-  cut(rnorm(100, mean=y, sd=0.25),c(-4,-1.5,0,1.5,4))
reg <- lm(y ~ x)
estimable(reg, c(    0,   1,    0,   -1) )

# Fit a spline with a single knot at 0.5 and plot the *pointwise*
# confidence intervals
library(gplots)
x2 <- rnorm(100,mean=y,sd=0.5)
reg2 <- lm(y ~ x + x2 + pmax(x2-0.5,0) )
range <- seq(-2,2,,50)
tmp <- estimable(reg2,cbind(1,0,0,1,range,pmax(range-0.5,0)), conf.int=0.95)
plotCI(x=range,y=tmp[,1],li=tmp[,6],ui=tmp[,7])

# Fit both linear and quasi-Poisson models to iris data, then compute
# conficence intervals on a contrast.
data(iris)
lm1  <- lm(Sepal.Length~Sepal.Width+Species+Sepal.Width:Species, data=iris)
glm1 <- glm(Sepal.Length~Sepal.Width+Species+Sepal.Width:Species, data=iris,
            family=quasipoisson("identity"))
cm <- rbind( lambda1 = c(1,0,1,0,0,0),
             lambda2 = c(1,0,0,1,0,0))
estimable(lm1,cm)
estimable(glm1,cm)
}
\keyword{ models }
\keyword{ regression }

