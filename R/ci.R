#' Compute Confidence Intervals
#' 
#' Compute and display confidence intervals for model estimates.  Methods are
#' provided for the mean of a numeric vector \code{ci.default}, the probability
#' of a binomial vector \code{ci.binom}, and for \code{lm}, \code{lme}, and
#' \code{mer} objects are provided.
#' 
#' 
#' @aliases ci ci.numeric ci.binom ci.lm ci.lme ci.estimable ci.fit_contrast
#' @param x object from which to compute confidence intervals.
#' @param confidence confidence level. Defaults to 0.95.
#' @param alpha type one error rate.  Defaults to 1.0-\code{confidence}
#' @param \dots Arguments for methods
#' @return vector or matrix with one row per model parameter and
#' elements/columns \code{Estimate}, \code{CI lower}, \code{CI upper},
#' \code{Std. Error}, \code{DF} (for lme objects only), and \code{p-value}.
#' @author Gregory R. Warnes \email{greg@@warnes.net}
#' @seealso \code{\link[stats]{confint}}, \code{\link[stats]{lm}},
#' \code{\link[stats]{summary.lm}}
#' @keywords regression
#' 
#' @examples
#' 
#' 
#' # mean and confidence interval
#' ci( rnorm(10) )
#' 
#' # binomial proportion and exact confidence interval
#' b <- rbinom( prob=0.75, size=1, n=20 )
#' ci.binom(b) # direct call
#' class(b) <- 'binom'
#' ci(b)       # indirect call
#' 
#' # confidence intervals for regression parameteres
#' data(state)
#' reg  <-  lm(Area ~ Population, data=as.data.frame(state.x77))
#' ci(reg)
#' 
#' @export
ci  <-  function(x, confidence=0.95,alpha=1-confidence,...)
  UseMethod("ci")

#' @rdname ci
#' @param na.rm `logical` indicating whether missing values should be removed.
#' @exportS3Method gmodels::ci
#' @importFrom stats qt
#' @importFrom stats sd
ci.numeric <- function(x, confidence=0.95,alpha=1-confidence,na.rm=FALSE,...)
{
  warning("No class or unkown class.  Using default calcuation.")
  est <- mean(x, na.rm=na.rm)
  stderr <-  sd(x, na.rm=na.rm)/sqrt(nobs(x));
  ci.low  <- est + qt(alpha/2, nobs(x)-1) * stderr
  ci.high <- est - qt(alpha/2, nobs(x)-1) * stderr
  retval  <- c(
    Estimate=est,
    "CI lower"=ci.low,
    "CI upper"=ci.high,
    "Std. Error"=stderr
  )
  
  retval
}

#' @export ci.binom
#' @exportS3Method gmodels::ci
#' @importFrom gdata nobs
#' @importFrom stats qbeta
#' @importFrom stats qt
ci.binom <- function(x, confidence=0.95,alpha=1-confidence,...)
{
  if( !(all(x %in% c(0,1))) ) stop("Binomial values must be either 0 or 1.")
  if( all(x==0) || all(x==1) )
    warning("All observed values are ", as.numeric(x[1]), ", so estimated Std. Error is 0.")
  
  est  <-  mean(x, na.rm=TRUE)
  n <- nobs(x)
  x <- sum(x)
  stderr <- sqrt(est*(1-est)/n)
  
  ci.low  <- qbeta(   alpha/2, x  , n + 1 - x)
  ci.high <- qbeta(1- alpha/2, x+1, n-x      )
  
  retval  <- cbind(Estimate=est,
                   "CI lower"=ci.low,
                   "CI upper"=ci.high,
                   "Std. Error"= stderr
  )
  retval
}

#' @exportS3Method gmodels::ci
#' @importFrom stats coef
ci.lm  <-  function(x,confidence=0.95,alpha=1-confidence,...)
{
  x  <-  summary(x)
  est  <-  coef(x)[,1] ;
  ci.low  <- est + qt(alpha/2, x$df[2]) * coef(x)[,2] ;
  ci.high <- est - qt(alpha/2, x$df[2]) * coef(x)[,2] ;
  retval  <- cbind(Estimate=est,
                   "CI lower"=ci.low,
                   "CI upper"=ci.high,
                   "Std. Error"= coef(x)[,2],
                   "p-value" = coef(x)[,4])
  
  retval
}

#' @exportS3Method gmodels::ci
#' @importFrom stats qt
ci.lme <- function(x,confidence=0.95,alpha=1-confidence,...)
{
  x  <-  summary(x)
  est  <-  x$tTable[,"Value"] ;
  ci.low  <- est + qt(alpha/2, x$tTable[,"DF"]) * x$tTable[,"Std.Error"] ;
  ci.high <- est - qt(alpha/2, x$tTable[,"DF"]) * x$tTable[,"Std.Error"] ;
  retval  <- cbind(Estimate=est,
                   "CI lower"=ci.low,
                   "CI upper"=ci.high,
                   "Std. Error"= x$tTable[,"Std.Error"],
                   "DF" = x$tTable[,"DF"],
                   "p-value" = x$tTable[,"p-value"])
  rownames(retval)  <-  rownames(x$tTable)
  retval
}

## ci.mer <- function (x,
##                     confidence = 0.95,
##                     alpha = 1 - confidence,
##                     n.sim = 1e4,
##                     ...)
## {
##     x.effects <- x@fixef
##     n <- length(x.effects)

##     retval <- gmodels::est.mer(obj = x,
##                                 cm = diag(n),
##                                 beta0 = rep(0, n),
##                                 conf.int = confidence,
##                                 show.beta0 = FALSE,
##                                 n.sim = n.sim)

##     retval <- retval[,
##                      c("Estimate", "Lower.CI", "Upper.CI", "Std. Error", "p value"),
##                      drop=FALSE
##                      ]
##     colnames(retval)[c(2:3, 5)] <- c("CI lower", "CI upper", "p-value")
##     rownames(retval) <- names(x.effects)

##     retval
## }


#' @exportS3Method gmodels::ci
#' @importFrom stats qt
ci.estimable  <-  function(x,confidence=0.95,alpha=1-confidence,...)
{
  ci.low  <- x$Estimate + qt(alpha/2, x$DF) * x$"Std. Error"
  ci.high <- x$Estimate - qt(alpha/2, x$DF) * x$"Std. Error"
  retval  <- cbind(Estimate=x$Estimate,
                   "CI lower"=ci.low,
                   "CI upper"=ci.high,
                   "Std. Error"= x$"Std. Error",
                   "p-value" = x$"Pr(>|t|)"
  )
  rownames(retval) <- rownames(x)
  
  retval
}


#' @exportS3Method gmodels::ci
ci.fit_contrast <- function (x, confidence = 0.95, alpha = 1 - confidence, ...)
{
  if( !all(c("lower CI", "upper CI") %in% colnames(x) ) )
    stop("object does not contain confidence interval information.")
  colnames(x) <- c("Estimate", "Std. Error", "Delete",
                   "p-value",
                   "CI lower", "CI upper")
  x[, c(1, 5:6, 2, 4), drop=FALSE]
}
