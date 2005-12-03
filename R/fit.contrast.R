# $Id: fit.contrast.R,v 1.6 2005/12/04 06:12:22 warnes Exp $

fit.contrast.lm <- function(model, varname, coeff, showall=FALSE,
                            conf.int=NULL, df=FALSE, ...)
{
  # check class of model
  if( !(any(class(model) %in% c("lm", "aov", "lme", "lmer") ) ))
    stop("contrast.lm can only be applied to objects inheriting from 'lm'",
         "and 'lme' (eg: lm,glm,aov,lme,lmer).")

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
  if(!("lmer" %in% class(model)))
     m <- model$call
  else
     m <- model@call # lmer is class 4

  if(is.null(m$contrasts))
    m$contrasts <- list()
  m$contrasts[[varname]] <- cmat

  if(is.R())
    r <- eval(m, parent.frame())
  else
    r <- eval(m)

  # now return the correct elements ....
  if( 'lmer' %in% class(model) )
    {
      est <- fixef(r) 
      se  <- sqrt(diag(vcov(r)))
      tval <- est/se
      df.lmer   <- getFixDF(r)
      retval <- cbind(
                      "Estimate"= est,
                      "Std. Error"= se,
                      "t-value"= tval,
                      "Pr(>|t|)"=  2 * (1 - pt(abs(tval), df.lmer)),
                      "DF"=df.lmer
                      )

    }
  else if( 'lme' %in% class(model) )
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
          if(!("lmer" %in% class(model))) # doesn't add this on in lmer
            rn <- paste(varname,rownames(coeff),sep="")
          else
            rn <- varname
          ind <- match(rn,rownames(retval))
          retval <- retval[ind,,drop=FALSE]
        }

    }

  if(!is.null(conf.int)) # add confidence intervals
    {
      alpha <- 1-conf.int
      retval <- cbind( retval,
                      "lower CI"=retval[,1] -
                      qt(1-alpha/2,retval[,5])*retval[,2],
                      "upper CI"=retval[,1] +
                      qt(1-alpha/2,retval[,5])*retval[,2] )
    }

  if(!df)
    return(retval[,-5,drop=FALSE])
  else
    return(retval)
}

# fit.contrast.lme and fit.contrast.lmer are necessary because
# 'lme' and 'lmer' objects do not inherit from 'lm'.
#
# **Make sure that the argument list *exactly* matches the one
# for fit.contrast.lm() above.**
#
fit.contrast.lme <- function(model, varname, coeff, showall=FALSE,
                            conf.int=NULL, df=FALSE, ...)
  {
    require(nlme)
    fit.contrast.lm(model, varname, coeff, showall, conf.int, df)
  }

fit.contrast.lmer <- function(model, varname, coeff, showall=FALSE,
                            conf.int=NULL, df=FALSE, ...)
  {
    require(lme4)
    fit.contrast.lm(model, varname, coeff, showall, conf.int, df)
  }


fit.contrast <- function(model, varname, coeff, ...)
  UseMethod("fit.contrast")

