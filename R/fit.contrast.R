#' Compute and test arbitrary contrasts for regression objects
#' 
#' Compute and test arbitrary contrasts for regression objects.
#' 
#' Computes the specified contrast(s) by re-fitting the model with the
#' appropriate arguments.  A contrast of the form \code{c(1,0,0,-1)} would
#' compare the mean of the first group with the mean of the fourth group.
#' 
#' @aliases fit.contrast fit.contrast.lm fit.contrast.lme
#' @param model regression (lm,glm,aov,lme) object for which the contrast(s)
#' will be computed.
#' @param varname variable name
#' @param coeff vector or matrix specifying contrasts (one per row).
#' @param showall return all regression coefficients. If \code{TRUE}, all model
#' cofficients will be returned.  If \code{FALSE} (the default), only the
#' coefficients corresponding to the specified contrast will be returned.
#' @param conf.int numeric value on (0,1) or NULL.  If a numeric value is
#' specified, confidence intervals with nominal coverage probability
#' \code{conf.int} will be computed.  If \code{NULL}, confidence intervals will
#' not be computed.
#' @param df boolean indicating whether to return a column containing the
#' degrees of freedom.
#' @param \dots optional arguments provided by methods.
#' @return Returns a matrix containing estimated coefficients, standard errors,
#' t values, two-sided p-values. If \code{df} is TRUE, an additional column
#' containing the degrees of freedom is included.  If \code{conf.int} is
#' specified lower and upper confidence limits are also returned.
#' @author Gregory R. Warnes \email{greg@@warnes.net}
#' @seealso \code{\link{lm}}, \code{\link{contrasts}},
#' \code{\link{contr.treatment}}, \code{\link{contr.poly}}, Computation and
#' testing of General Linear Hypothesis: \code{\link{glh.test}}, Computation
#' and testing of estimable functions of model coefficients:
#' \code{\link{estimable}}, \code{\link{make.contrasts}}
#' @references Venables & Ripley, Section 6.2
#' @keywords models regression
#' @examples
#' 
#' y <- rnorm(100)
#' x <-  cut(rnorm(100, mean=y, sd=0.25),c(-4,-1.5,0,1.5,4))
#' reg <- lm(y ~ x)
#' summary(reg)
#' 
#' # look at the group means
#' gm <- sapply(split(y,x),mean)
#' gm
#' 
#' 
#' # mean of 1st group vs mean of 4th group
#' fit.contrast(reg, x, c(    1,    0,    0,   -1) )
#' # estimate should be equal to:
#' gm[1] - gm[4]
#' 
#' # mean of 1st and 2nd groups vs mean of 3rd and 4th groups
#' fit.contrast(reg, x, c( -1/2, -1/2,  1/2,  1/2) )
#' # estimate should be equal to:
#' sum(-1/2*gm[1], -1/2*gm[2], 1/2*gm[3], 1/2*gm[4])
#' 
#' # mean of 1st group vs mean of 2nd, 3rd and 4th groups
#' fit.contrast(reg, x, c( -3/3,  1/3,  1/3,  1/3) )
#' # estimate should be equal to:
#' sum(-3/3*gm[1], 1/3*gm[2], 1/3*gm[3], 1/3*gm[4])
#' 
#' # all at once
#' cmat <- rbind( "1 vs 4"    =c(-1, 0, 0, 1),
#'                "1+2 vs 3+4"=c(-1/2,-1/2, 1/2, 1/2),
#'                "1 vs 2+3+4"=c(-3/3, 1/3, 1/3, 1/3))
#' fit.contrast(reg,x,cmat)
#' 
#' #
#' x2 <- rnorm(100,mean=y,sd=0.5)
#' reg2 <- lm(y ~ x + x2 )
#' fit.contrast(reg2,x,c(-1,0,0,1))
#' 
#' #
#' # Example for Analysis of Variance
#' #
#' 
#' set.seed(03215)
#' Genotype <- sample(c("WT","KO"), 1000, replace=TRUE)
#' Time <- factor(sample(1:3, 1000, replace=TRUE))
#' y <- rnorm(1000)
#' data <- data.frame(y, Genotype, Time)
#' 
#' 
#' # Compute Contrasts & obtain 95% confidence intervals
#' 
#' model <- aov( y ~ Genotype + Time + Genotype:Time, data=data )
#' 
#' fit.contrast( model, "Genotype", rbind("KO vs WT"=c(-1,1) ), conf=0.95 )
#' 
#' fit.contrast( model, "Time",
#'             rbind("1 vs 2"=c(-1,1,0),
#'                   "2 vs 3"=c(0,-1,1)
#'                   ),
#'             conf=0.95 )
#' 
#' 
#' cm.G <- rbind("KO vs WT"=c(-1,1) )
#' cm.T <- rbind("1 vs 2"=c(-1,1,0),
#'               "2 vs 3"=c(0,-1,1) )
#' 
#' # Compute contrasts and show SSQ decompositions
#' 
#' model <- aov( y ~ Genotype + Time + Genotype:Time, data=data,
#'               contrasts=list(Genotype=make.contrasts(cm.G),
#'                              Time=make.contrasts(cm.T) )
#'             )
#' 
#' summary(model, split=list( Genotype=list( "KO vs WT"=1 ),
#'                            Time = list( "1 vs 2" = 1,
#'                                         "2 vs 3" = 2 ) ) )
#' 
#' 
#' # example for lme
#' library(nlme)
#' data(Orthodont)
#' fm1 <- lme(distance ~ Sex, data = Orthodont,random=~1|Subject)
#' 
#' # Contrast for sex.  This example is equivalent to standard treatment
#' # contrast.
#' #
#' fit.contrast(fm1, "Sex", c(-1,1), conf.int=0.95 )
#' #
#' # and identical results can be obtained using lme built-in 'intervals'
#' #
#' intervals(fm1)
#' 
#' # Cut age into quantile groups & compute some contrasts
#' Orthodont$AgeGroup <- gtools::quantcut(Orthodont$age)
#' fm2 <- lme(distance ~ Sex + AgeGroup, data = Orthodont,random=~1|Subject)
#' #
#' fit.contrast(fm2, "AgeGroup", rbind("Linear"=c(-2,-1,1,2),
#'                                     "U-Shaped"=c(-1,1,1,-1),
#'                                     "Change-Point at 11"=c(-1,-1,1,1)),
#'                               conf.int=0.95)
#' 
#' 
#' @export
fit.contrast <- function(model, varname, coeff, ...)
  UseMethod("fit.contrast")

#' @exportS3Method stats::coef
coef.fit_contrast <- function(object, ...)
  object

#' @exportS3Method base::print
print.fit_contrast <- function(object, ...)
  print(unclass(object))

#' @exportS3Method gmodels::fit.contrast
#' @importFrom stats coef
#' @importFrom stats pt
#' @importFrom stats qt
#' @importFrom stats summary.glm
#' @importFrom stats summary.lm
fit.contrast.lm <- function(model, varname, coeff, showall=FALSE,
                            conf.int=NULL,
                            df=FALSE,
                            ...)
{
  # check class of model
  if( !(any(class(model) %in% c("lm", "aov", "lme") ) ))
    stop("contrast.lm can only be applied to objects inheriting from 'lm'",
         "and 'lme' (eg: lm,glm,aov,lme).")
  
  # make sure we have the NAME of the variable
  if(!is.character(varname))
    varname <- deparse(substitute(varname))
  
  # make coeff into a matrix
  if(!is.matrix(coeff))
  {
    coeff <- matrix(coeff, nrow=1)
  }
  
  # make sure columns are labeled
  if (is.null(rownames(coeff)))
  {
    rn <- vector(length=nrow(coeff))
    for(i in 1:nrow(coeff))
      rn[i] <- paste(" c=(",paste(coeff[i,],collapse=" "), ")")
    rownames(coeff) <- rn
  }
  
  # now convert into the proper form for the contrast matrix
  cmat <- make.contrasts(coeff, ncol(coeff) )
  cn <- paste(" C",1:ncol(cmat),sep="")
  cn[1:nrow(coeff)] <- rownames(coeff)
  colnames(cmat) <- cn
  
  # recall fitting method with the specified contrast
  m <- model$call
  
  if(is.null(m$contrasts))
    m$contrasts <- list()
  m$contrasts[[varname]] <- cmat
  
  if(is.R())
    r <- eval(m, parent.frame())
  else
    r <- eval(m)
  
  # now return the correct elements ....
  if( 'lme' %in% class(model) )
  {
    est <- r$coefficients$fixed
    se  <- sqrt(diag(r$varFix))
    tval <- est/se
    df.lme   <- r$fixDF$X
    retval <- cbind(
      "Estimate"= est,
      "Std. Error"= se,
      "t-value"= tval,
      "Pr(>|t|)"=  2 * (1 - pt(abs(tval), df.lme)),
      "DF"=df.lme
    )
    
  }
  else if('glm' %in% class(model))
  {
    smodel <- summary.glm(r)
    retval <- cbind(coef(smodel), "DF"=smodel$df[2])
  }
  else # lm, aov
  {
    smodel <- summary.lm(r)
    retval <- cbind(coef(smodel), "DF"=smodel$df[2])
  }
  
  if( !showall )
  {
    if( !is.R() && ncol(cmat)==1 )
    {
      retval <- retval[varname,,drop=FALSE]
      rownames(retval) <- rn
    }
    else
    {
      rn <- paste(varname,rownames(coeff),sep="")
      ind <- match(rn,rownames(retval))
      retval <- retval[ind,,drop=FALSE]
    }
    
  }
  
  if(!missing(conf.int) && !is.null(conf.int)) # add confidence intervals
  {
    alpha <- 1-conf.int
    retval <- cbind( retval,
                     "lower CI"=retval[,1] -
                       qt(1-alpha/2,retval[,5])*retval[,2],
                     "upper CI"=retval[,1] +
                       qt(1-alpha/2,retval[,5])*retval[,2] )
  }
  
  if(!df)
    retval <- retval[,-5,drop=FALSE]
  
  class(retval) <- "fit_contrast"
  
  retval
}

# fit.contrast.lme and fit.contrast.mer are necessary because
# 'lme' and 'mer' objects do not inherit from 'lm'.
#
# **Make sure that the argument list *exactly* matches the one
# for fit.contrast.lm() above.**
#
#' @exportS3Method gmodels::fit.contrast
fit.contrast.lme <- function(model, varname, coeff, showall=FALSE,
                             conf.int=NULL, df=FALSE, ...)
{
  fit.contrast.lm(model, varname, coeff, showall, conf.int, df)
}

## # I made rather dramatic changes here and do all calculations in fit.contrast.mer rather than
## # fit.contrast.lm because of the simulation extras ... added sim.mer and n.sim to the parameter list
## fit.contrast.mer <- function(model, varname, coeff, showall=FALSE,
##                             conf.int=NULL, sim.mer=TRUE, n.sim=1000, ...)
## {
##   require(lme4)

##   # make sure we have the NAME of the variable
##   if(!is.character(varname))
##      varname <- deparse(substitute(varname))

##   # make coeff into a matrix
##   if(!is.matrix(coeff))
##     {
##        coeff <- matrix(coeff, nrow=1)
##      }

##   # make sure columns are labeled
##   if (is.null(rownames(coeff)))
##      {
##        rn <- vector(length=nrow(coeff))
##        for(i in 1:nrow(coeff))
##           rn[i] <- paste(" c=(",paste(coeff[i,],collapse=" "), ")")
##        rownames(coeff) <- rn
##      }

##   # now convert into the proper form for the contrast matrix
##   cmat <- make.contrasts(coeff, ncol(coeff) )
##   cn <- paste(" C",1:ncol(cmat),sep="")
##   cn[1:nrow(coeff)] <- rownames(coeff)
##   colnames(cmat) <- cn

##   m <- model@call

##   if(is.null(m$contrasts))
##     m$contrasts <- list()
##   m$contrasts[[varname]] <- cmat

##   if(is.R())
##     r <- eval(m, parent.frame())
##   else
##     r <- eval(m)
##   # now return the correct elements ....
##   r.effects <- fixef(r)
##   n <- length(r.effects)

##   if(sim.mer)
##   {
##     retval <- est.mer(obj = r, cm = diag(n), beta0 = rep(0, n),
##                        conf.int = conf.int, show.beta0 = FALSE,
##                        n.sim=n.sim)
##     rownames(retval) <- names(r.effects)
##   }else{
##     if(!is.null(conf.int))
##       warning("Confidence interval calculation for mer objects requires simulation -- use sim.mer = TRUE")

##     est <- fixef(r)
##     se  <- sqrt(diag(as.matrix(vcov(r))))
##     tval <- est/se
##     retval <- cbind(
##                     "Estimate"= est,
##                     "Std. Error"= se,
##                     "t-value"= tval
##                     )
##   }

##   if( !showall )
##   {
##     if( !is.R() && ncol(cmat)==1 )
##     {
##       retval <- retval[varname,,drop=FALSE]
##       rownames(retval) <- rn
##     }else{
##       rn <- paste(varname,rownames(coeff),sep="")
##       ind <- match(rn,rownames(retval))
##       retval <- retval[ind,,drop=FALSE]
##     }
##   }

##   return(retval)
## }


