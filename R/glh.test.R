#' Test a General Linear Hypothesis for a Regression Model
#' 
#' Test, print, or summarize a general linear hypothesis for a regression model
#' 
#' Test the general linear hypothesis \eqn{C \hat{\beta} = d } for the regression model `reg`.
#' 
#' The test statistic is obtained from the formula: 
#' \deqn{
#'  f = \frac{(C \hat{\beta} - d)' ( C (X'X)^{-1} C' ) (C \hat{\beta} - d) / r }{ 
#'      SSE / (n-p) } 
#'  }
#' where
#' * `r` is the number of contrasts contained in `C`, and
#' * `n-p` is the model degrees of freedom.
#'  
#' Under the null hypothesis, `f` will follow a F-distribution with `r` and `n-p`
#' degrees of freedom
#' 
#' @aliases glh.test print.glh.test summary.glh.test
#' @param reg Regression model
#' @param cm contrast matrix `C` .  Each row specifies a linear combination of the
#' coefficients
#' @param d vector `d` specifying the null hypothesis values for each linear
#' combination
#' 
#' @return Object of class `c("glh.test","htest")` with elements:
#' \item{call }{Function call that created the object}
#' \item{statistic }{F statistic} 
#' \item{parameter}{vector containing the numerator (r) and 
#'  denominator (n-p) degrees of freedom} 
#' \item{p.value}{p-value}
#' \item{estimate}{computed estimate for each row of `cm`}
#' \item{null.value}{d} 
#' \item{method}{description of the method}
#' \item{data.name}{name of the model given for `reg`}
#' \item{matrix}{matrix specifying the general linear hypotheis (`cm`)}
#' 
#' @note When using treatment contrasts (the default) the first level of the
#' factors are subsumed into the intercept term.  The estimated model
#' coefficients are then contrasts versus the first level. This should be taken
#' into account when forming contrast matrixes, particularly when computing
#' contrasts that include 'baseline' level of factors.
#' 
#' See the comparison with `fit.contrast` in the examples below.
#' @author Gregory R. Warnes \email{greg@@warnes.net}
#' @seealso [fit.contrast()], [estimable()], [stats::contrasts()]
#' @references R.H. Myers, Classical and Modern Regression with Applications,
#' 2nd Ed, 1990, p. 105
#' @keywords models regression
#' @examples
#' 
#' 
#' # fit a simple model
#' y <- rnorm(100)
#' x <-  cut(rnorm(100, mean=y, sd=0.25),c(-4,-1.5,0,1.5,4))
#' reg <- lm(y ~ x)
#' summary(reg)
#' 
#' # test both group 1 = group 2  and group 3 = group 4
#' # *Note the 0 in the column for the intercept term*
#' 
#' C <- rbind( c(0,-1,0,0), c(0,0,-1,1) )
#' ret <- glh.test(reg, C)
#' ret  # same as 'print(ret) '
#' summary(ret)
#' 
#' # To compute a contrast between the first and second level of the factor
#' # 'x' using 'fit.contrast' gives:
#' 
#' fit.contrast( reg, x,c(1,-1,0,0) )
#' 	
#' # To test this same contrast using 'glh.test', use a contrast matrix
#' # with a zero coefficient for the intercept term.  See the Note section,
#' # above, for an explanation.
#' 
#' C <- rbind( c(0,-1,0,0) )
#' glh.test( reg, C )
#' 
#' @importFrom stats pf
#' @importFrom stats coef
#' @importFrom stats summary.lm
#' 
#' @export
glh.test <- function(reg, cm, d=rep(0, nrow(cm)) )
{

  if( !is.matrix(cm) && !is.data.frame(cm) )
    cm <- matrix(cm, nrow=1)

  if ( !( "lm" %in% class(reg) ) )
    stop("Only defined for lm,glm objects")

  bhat <- summary.lm(reg)$coefficients[,1,drop=FALSE]
  XpX <- summary.lm(reg)$cov.unscaled
  df <- reg$df.residual
  msr <- summary.lm(reg)$sigma  # == SSE / (n-p)
  r <- nrow(cm)


  if ( ncol(cm) != length(bhat) ) stop(
                   paste( "\n Dimension of ",
                         deparse( substitute( cm ) ), ": ",
                         paste( dim(cm), collapse="x" ),
                         ", not compatible with no of parameters in ",
                         deparse( substitute( reg ) ), ": ",
                         length(bhat), sep="" ) )


  #                        -1
  #     (Cb - d)' ( C (X'X)   C' ) (Cb - d) / r
  # F = ---------------------------------------
  #                 SSE / (n-p)
  #

  Fstat <- t(cm %*% bhat - d) %*% solve((cm %*% XpX %*% t(cm))) %*% (cm %*% bhat - d) / r / msr^2

  p <- 1-pf(Fstat,r,df)

  retval <- list()
  retval$call <- match.call()
  retval$statistic <- c(F=Fstat)
  retval$parameter <- c(df1=r,df2=df)
  retval$p.value <- p
  retval$conf.int <- NULL
  retval$estimate <- cm%*%bhat
  retval$null.value <- d
  retval$method <- "Test of General Linear Hypothesis"
  retval$data.name <- deparse(substitute(reg))
  retval$matrix <- cm
  colnames(retval$matrix) <- names(reg$coef)

  class(retval) <- c("glh.test","htest")

  retval
}

#' @exportS3Method base::print
#' @inheritParams base::print
print.glh.test <- function(x, digits = 4, ... )
{
    cat("\n")
    cat("\t",x$method, prefix = "\t")
    cat("\n")
    cat("Call:\n")
    print(x$call)

    if (!is.null(x$statistic))
        cat(names(x$statistic), " = ", format(round(x$statistic,
            4)), ", ", sep = "")
    if (!is.null(x$parameter))
        cat(paste(names(x$parameter), " = ", format(round(x$parameter,
            3)), ",", sep = ""), "")
    cat("p-value =",
        format.pval(x$p.value, digits = digits),
        "\n")
    cat("\n")
  }


#' @exportS3Method base::summary
#' @inheritParams base::summary
summary.glh.test <- function(object, digits = 4, ... )
{
    cat("\n")
    cat("\t",object$method, prefiobject = "\t")
    cat("\n")
    cat("Regression: ", object$data.name, "\n")
    cat("\n")
    cat("Null Hypothesis: C %*% Beta-hat = d \n")
    cat("\n")
    cat("C matrix: \n")
    print(object$matrix, digits=digits)
    cat("\n")
    cat("d vector: \n")
    print(object$null.value, digits=digits)
    cat("\n")
    cat("C %*% Beta-hat: \n")
    print(c(object$estimate))
    cat("\n")

    if (!is.null(object$statistic))
        cat(names(object$statistic), " = ", format(round(object$statistic,
            4)), ", ", sep = "")
    if (!is.null(object$parameter))
        cat(paste(names(object$parameter), " = ", format(round(object$parameter,
            3)), ",", sep = ""), "")
    cat("p-value =",
        format.pval(object$p.value, digits = digits),
        "\n")
    cat("\n")
  }



