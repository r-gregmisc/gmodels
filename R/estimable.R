#' Contrasts and estimable linear functions of model coefficients
#' 
#' Compute and test contrasts and other estimable linear functions of model
#' coefficients for for lm, glm, lme, mer, and geese objects
#' 
#' `estimable` computes an estimate, test statitic, significance test, and
#' (optional) confidence interval for each linear functions of the model
#' coefficients specified by `cm`.
#' 
#' The estimable function(s) may be specified via a vector, list, or matrix.
#' If `cm` is a vector, it should contained named elements each of which
#' gives the coefficient to be applied to the corresponding parameter. These
#' coefficients will be used to construct the contrast matrix, with unspecified
#' model parameters assigned zero coefficients. If `cm` is a list, it
#' should contain one or more coefficient vectors, which will be used to
#' construct rows of the contrast matrix.  If `cm` is a matrix, column
#' names must match (a subset of) the model parameters, and each row should
#' contain the corresponding coefficient to be applied.  Model parameters which
#' are not present in the set of column names of `cm` will be set to zero.
#' 
#' The estimates and their variances are obtained by applying the contrast
#' matrix (generated from) `cm` to the model estimates variance-covariance
#' matrix.  Degrees of freedom are obtained from the appropriate model terms.
#' 
#' The user is responsible for ensuring that the specified linear functions are
#' meaningful.
#' 
#' For computing contrasts among levels of a single factor, `fit.contrast`
#' may be more convenient.  For computing contrasts between two specific
#' combinations of model parameters, the `contrast` function in Frank
#' Harrell's 'rms' library (formerly 'Design') may be more convenient.
#' 
#' %The `.wald` function is called internally by `estimable` and %is
#' not intended for direct use.
#' 
#' @aliases estimable estimable.default estimable.mlm
#' @param obj Regression (lm, glm, lme, mer, mlm) object.
#' @param cm Vector, List, or Matrix specifying estimable linear functions or
#' contrasts.  See below for details.
#' @param beta0 Vector of null hypothesis values
#' @param conf.int Confidence level.  If provided, confidence intervals will be
#' computed.
#' @param show.beta0 Logical value. If TRUE a column for beta0 will be included
#' in the output table.  Defaults to TRUE when beta0 is specified, FALSE
#' otherwise.
#' @param ... ignored
#' @return Returns a matrix with one row per linear function.  Columns contain
#' the beta0 value (optional, see `show.beta0` above), estimated
#' coefficients, standard errors, t values, degrees of freedom, two-sided
#' p-values, and the lower and upper endpoints of the 1-alpha confidence
#' intervals.
#' @note The estimated fixed effect parameters of `lme` objects may have
#' different degrees of freedom.  If a specified contrast includes nonzero
#' coefficients for parameters with differing degrees of freedom, the smallest
#' number of degrees of freedom is used and a warning message is issued.
#' @author BXC (Bendix Carstensen) \email{b@@bxc.dk}, Gregory R. Warnes
#' \email{greg@@warnes.net}, Soren Hojsgaard \email{sorenh@@agrsci.dk}, and
#' Randall C Johnson \email{rjohnson@@ncifcrf.gov}
#' @seealso [fit.contrast()], [stats::lm()],
#' [nlme::lme()], [stats::contrasts()],
#' [rms::contrast()]
#' @keywords models regression
#' @examples
#' 
#' # setup example data
#' y <- rnorm(100)
#' x <-  cut(rnorm(100, mean=y, sd=0.25),c(-4,-1.5,0,1.5,4))
#' levels(x) <- c("A","B","C","D")
#' x2 <- rnorm(100, mean=y, sd=0.5)
#' 
#' # simple contrast and confidence interval
#' reg <- lm(y ~ x)
#' estimable(reg, c(    0,   1,    0,   -1) )  # full coefficient vector
#' estimable(reg, c("xB"=1,"xD"=-1) )          # just the nonzero terms
#' 
#' 
#' # Fit a spline with a single knot at 0.5 and plot the *pointwise*
#' # confidence intervals
#' library(gplots)
#' pm <- pmax(x2-0.5, 0) # knot at 0.5
#' reg2 <- lm(y ~ x + x2 + pm )
#' 
#' range <- seq(-2, 2, , 50)
#' tmp <- estimable(reg2,
#'                  cm=cbind(
#'                           '(Intercept)'=1,
#'                           'xC'=1,
#'                           'x2'=range,
#'                           'pm'=pmax(range-0.5, 0)
#'                            ),
#'                  conf.int=0.95)
#' plotCI(x=range, y=tmp[, 1], li=tmp[, 6], ui=tmp[, 7])
#' 
#' # Fit both linear and quasi-Poisson models to iris data, then compute
#' # joint confidence intervals on contrasts for the Species and
#' # Sepal.Width by Species interaction terms.
#' data(iris)
#' lm1  <- lm (Sepal.Length ~ Sepal.Width + Species + Sepal.Width:Species, data=iris)
#' glm1 <- glm(Sepal.Length ~ Sepal.Width + Species + Sepal.Width:Species, data=iris,
#'             family=quasipoisson("identity"))
#' 
#' cm <- rbind(
#'             'Setosa vs. Versicolor'   = c(0, 0, 1, 0, 1, 0),
#'             'Setosa vs. Virginica'    = c(0, 0, 0, 1, 0, 1),
#'             'Versicolor vs. Virginica'= c(0, 0, 1,-1, 1,-1)
#'             )
#' estimable(lm1, cm)
#' estimable(glm1, cm)
#' 
#' @export
estimable <- function (obj, cm, beta0, conf.int=NULL,  show.beta0, ...)
  {
    UseMethod("estimable")
  }

#' @rdname estimable
#' @param joint.test Logical value. If TRUE a 'joint' Wald test for the
#'   hypothesis \eqn{L \beta=\beta_0} is performed.  
#'   Otherwise 'row-wise' tests are performed, i.e. \eqn{(L \beta)[i]=\beta_0[i]}.
#' @exportS3Method gmodels::estimable
#' @importFrom stats coef
#' @importFrom stats family
#' @importFrom stats pchisq
#' @importFrom stats pt
#' @importFrom stats qt
#' @importFrom stats summary.lm
#' 
estimable.default <- function (obj, cm, beta0, conf.int=NULL,
                               show.beta0, joint.test=FALSE, ...)
{
  if (is.matrix(cm) || is.data.frame(cm))
    {
      cm <- t(apply(cm, 1, .to.est, obj=obj))
    }
  else if(is.list(cm))
    {
      cm <- matrix(.to.est(obj, cm), nrow=1)
    }
  else if(is.vector(cm))
    {
      cm <- matrix(.to.est(obj, cm), nrow=1)
    }
  else
    {
      stop("`cm' argument must be of type vector, list, or matrix.")
    }

  if(missing(show.beta0))
    {
      if(!missing(beta0))
        show.beta0=TRUE
      else
        show.beta0=FALSE
    }


  if (missing(beta0))
    {
      beta0 = rep(0, ifelse(is.null(nrow(cm)), 1, nrow(cm)))

    }


  if (joint.test==TRUE)
    {
      .wald(obj, cm, beta0)
    }
  else
    {
      if ("lme" %in% class(obj)) {
        stat.name <- "t.stat"
        cf <- summary(obj)$tTable
        rho <- summary(obj)$cor
        vcv <- rho * outer(cf[, 2], cf[, 2])
        tmp <- cm
        tmp[tmp==0] <- NA
        df.all <- t(abs(t(tmp) * obj$fixDF$X))
        df <- apply(df.all, 1, min, na.rm=TRUE)
        problem <- apply(df.all !=df, 1, any, na.rm=TRUE)
        if (any(problem))
          warning(paste("Degrees of freedom vary among parameters used to ",
                        "construct linear contrast(s): ",
                        paste((1:nrow(tmp))[problem], collapse=", "),
                        ". Using the smallest df among the set of parameters.",
                        sep=""))
      }
      else if ("lm" %in% class(obj))
        {
          stat.name <- "t.stat"
          cf <- summary.lm(obj)$coefficients
          vcv <- summary.lm(obj)$cov.unscaled * summary.lm(obj)$sigma^2
          df <- obj$df.residual
          if ("glm" %in% class(obj))
            {
              vcv <- summary(obj)$cov.scaled
              if(family(obj)[1] %in% c("poisson", "binomial"))
                {
                  stat.name <- "X2.stat"
                  df <- 1
                }
              else
                {
                  stat.name <- "t.stat"
                  df <- obj$df.residual
                }
            }
        }
      else if ("geese" %in% class(obj))
        {
          stat.name <- "X2.stat"
          cf <- summary(obj)$mean
          vcv <- obj$vbeta
          df <- 1
        }
      else if ("gee" %in% class(obj))
        {
          stat.name <- "X2.stat"
          cf <- summary(obj)$coef
          vcv <- obj$robust.variance
          df <- 1
        }
      else
        {
          stop("obj must be of class 'lm', 'glm', 'aov', 'lme', 'gee', 'geese' or 'nlme'")
        }
      if (is.null(cm))
        cm <- diag(dim(cf)[1])
      if (!dim(cm)[2]==dim(cf)[1])
        stop(paste("\n Dimension of ", deparse(substitute(cm)),
                   ": ", paste(dim(cm), collapse="x"),
                   ", not compatible with no of parameters in ",
                   deparse(substitute(obj)), ": ", dim(cf)[1], sep=""))
      ct <- cm %*% cf[, 1]
      ct.diff <- cm %*% cf[, 1] - beta0

      vc <- sqrt(diag(cm %*% vcv %*% t(cm)))
      if (is.null(rownames(cm)))
        rn <- paste("(", apply(cm, 1, paste, collapse=" "),
                    ")", sep="")
      else rn <- rownames(cm)
      switch(stat.name,
             t.stat={
               prob <- 2 * (1 - pt(abs(ct.diff/vc), df))
             },
             X2.stat={
               prob <- 1 - pchisq((ct.diff/vc)^2, df=1)
             })

      if (stat.name=="X2.stat")
        {
          retval <- cbind(hyp=beta0, est=ct, stderr=vc,
                          "X^2 value"=(ct.diff/vc)^2,
                          df=df, prob=1 - pchisq((ct.diff/vc)^2, df=1))
          dimnames(retval) <- list(rn, c("beta0", "Estimate", "Std. Error",
                                         "X^2 value", "DF", "Pr(>|X^2|)"))
        }
      else if (stat.name=="t.stat")
        {
          retval <- cbind(hyp=beta0, est=ct, stderr=vc, "t value"=ct.diff/vc,
                          df=df, prob=2 * (1 - pt(abs(ct.diff/vc), df)))
          dimnames(retval) <- list(rn, c("beta0", "Estimate", "Std. Error",
                                         "t value", "DF", "Pr(>|t|)"))
        }

      if (!is.null(conf.int))
        {
          if (conf.int <=0 || conf.int >=1)
            stop("conf.int should be between 0 and 1. Usual values are 0.95, 0.90")
          alpha <- 1 - conf.int
          switch(stat.name,
                 t.stat={
                   quant <- qt(1 - alpha/2, df)
                 },
                 X2.stat={
                   quant <- qt(1 - alpha/2, 100)
                 })
          nm <- c(colnames(retval), "Lower.CI", "Upper.CI")
          retval <- cbind(retval, lower=ct.diff - vc * quant, upper=ct.diff +
                          vc * quant)
          colnames(retval) <- nm
        }
      rownames(retval) <- make.unique(rownames(retval))
      retval <- as.data.frame(retval)
      if(!show.beta0) retval$beta0 <- NULL

      class(retval) <- c("estimable", class(retval))

      return(retval)
    }
}

#' @importFrom stats coef
.wald <- function (obj, cm,
                   beta0=rep(0, ifelse(is.null(nrow(cm)), 1, nrow(cm))))
{
  if (!is.matrix(cm) && !is.data.frame(cm))
    cm <- matrix(cm, nrow=1)
  df <- nrow(cm)
  if ("geese" %in% class(obj))
    {
      cf <- obj$beta
      vcv <- obj$vbeta
    }
  else if ("gee" %in% class(obj))
    {
      cf <- summary(obj)$coef
      vcv <- obj$robust.variance
    }
  else if ("lm" %in% class(obj))
    {
      cf <- summary.lm(obj)$coefficients[, 1]
      if ("glm" %in% class(obj))
        vcv <- summary(obj)$cov.scaled
      else
        vcv <- summary.lm(obj)$cov.unscaled * summary.lm(obj)$sigma^2
    }
  else if ("lme" %in% class(obj))
    {
      s.o <- summary(obj)
      cf <- s.o$tTable[,1]
      se <- s.o$tTable[, 2]
      rho <- s.o$cor
      vcv <- rho * outer(se, se)
    }
  else
    stop("obj must be of class 'lm', 'glm', 'aov', 'gee', 'geese', or 'lme'.")
  u <- (cm %*% cf)-beta0
  vcv.u <- cm %*% vcv %*% t(cm)
  W <- t(u) %*% solve(vcv.u) %*% u
  prob <- 1 - pchisq(W, df=df)
  retval <- as.data.frame(cbind(W, df, prob))
  names(retval) <- c("X2.stat", "DF", "Pr(>|X^2|)")
  print(as.data.frame(retval))
}

## estimable.mer <- function (obj, cm, beta0, conf.int=NULL, show.beta0,
##                            sim.mer=TRUE, n.sim=1000, ...)
## {
##   if (is.matrix(cm) || is.data.frame(cm))
##     {
##       cm <- t(apply(cm, 1, .to.est, obj=obj))
##     }
##   else if(is.list(cm))
##     {
##       cm <- matrix(.to.est(obj, cm), nrow=1)
##     }
##   else if(is.vector(cm))
##     {
##       cm <- matrix(.to.est(obj, cm), nrow=1)
##     }
##   else
##     {
##       stop("'cm' argument must be of type vector, list, or matrix.")
##     }

##   if(missing(show.beta0))
##     {
##       if(!missing(beta0))
##         show.beta0=TRUE
##       else
##         show.beta0=FALSE
##     }


##   if (missing(beta0))
##     {
##       beta0 = rep(0, ifelse(is.null(nrow(cm)), 1, nrow(cm)))

##     }

##   if ("mer" %in% class(obj)) {
##     if(sim.mer)
##       return(est.mer(obj=obj, cm=cm, beta0=beta0, conf.int=conf.int,
##                      show.beta0=show.beta0, n.sim=n.sim))

##     stat.name <- "mer"
##     cf <- as.matrix(fixef(obj))
##     vcv <- as.matrix(vcov(obj))
##     df <- NA
##   }
##   else {
##     stop("obj is not of class mer")
##   }

##   if (is.null(rownames(cm)))
##     rn <- paste("(", apply(cm, 1, paste, collapse=" "),
##                 ")", sep="")
##   else rn <- rownames(cm)

##   ct <- cm %*% cf[, 1]
##   ct.diff <- cm %*% cf[, 1] - beta0
##   vc <- sqrt(diag(cm %*% vcv %*% t(cm)))

##   retval <- cbind(hyp=beta0, est=ct, stderr=vc, "t value"=ct.diff/vc)
##   dimnames(retval) <- list(rn, c("beta0", "Estimate", "Std. Error",
##                                  "t value"))

##   rownames(retval) <- make.unique(rownames(retval))
##   retval <- as.data.frame(retval)
##   if(!show.beta0) retval$beta0 <- NULL

##   class(retval) <- c("estimable", class(retval))

##   return(retval)

## }
