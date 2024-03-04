#' Display estimate, confidence interval and p-value for one model term
#' 
#' @param model model object
#' @param term model term 
#' @param mult  scale (multiply) the parameter by this factor 
#' @param digits number of significant digits to display
#' @param ... optional arguments
#' 
#' @examples
#' 
#' set.seed(42)
#' 
#' # fit an example model with 3 groups
#' y <- rnorm(100)
#' x <-  cut(rnorm(100, mean=y, sd=0.25),c(-4,-1.5,0,1.5,4))
#' reg <- lm(y ~ x)
#' reg
#' 
#' # show model estimate, p-value, and confidence interval
#' # for the first group
#' est_p_ci(reg, 2)
#' 
#' # estimate some group contrasts
#' cmat <- rbind( "1 vs 4"    =c(-1, 0, 0, 1),
#'                "1+2 vs 3+4"=c(-1/2,-1/2, 1/2, 1/2),
#'                "1 vs 2+3+4"=c(-3/3, 1/3, 1/3, 1/3))
#' cont <- fit.contrast(reg, x, cmat, conf.int = 0.95)
#' cont
#' 
#' # show the contrast estimate, p-value, and confidence interval
#' # for the first contrast
#' est_p_ci(cont, 2:3)
#' 
#' @export
est_p_ci <- function(model, term, mult=1, digits=2, ...)
  UseMethod("est_p_ci")

#' @exportS3Method gmodels::est_p_ci
est_p_ci.lm <- function(model, term, mult=1, digits=2, ...)
{
  info <- ci(model)
  if(is.character(term) && !(term %in% rownames(info)))
    stop(term, " is not a coefficient in model.")
  info <- info[term,,drop=FALSE]
  info.ci <- trimws(format( round(mult * info[,1:3, drop=FALSE], digits=digits) ))
  if(mult < 0)
    colnames(info.ci)[2:3] <- colnames(info.ci)[3:2]
  paste("Est=", info.ci[,'Estimate'],
        " ",
        "p=",format.pval(info[,'p-value'], digits=digits),
        " ",
        "95% CI: ",
        info.ci[,'CI lower'],
        " to ",
        info.ci[,'CI upper'],
        sep=""
  )
}

#' @exportS3Method gmodels::est_p_ci
est_p_ci.fit_contrast <- function(model, term, mult=1, digits=2, ...)
{
  if( !all(c("lower CI", "upper CI") %in% colnames(model) ) )
    stop("object does not contain confidence interval information.")

  if(is.character(term) && !(term %in% rownames(model)))
    stop(term, " is not a coefficient in model.")

  info.ci <- trimws(format( round(mult * model[term, c("Estimate", "lower CI", "upper CI")],
                           digits=digits) )  )
  if(mult < 0)
    colnames(info.ci)[2:3] <- colnames(info.ci)[3:2]
  
  paste("Est=", info.ci[,"Estimate"],
        " ",
        "p=",format.pval(model[term, 'Pr(>|t|)'], digits=digits),
        " ",
        "95% CI: ",
        info.ci[,"lower CI"],
        " to ",
        info.ci[,"upper CI"],
        sep=""
  )
}
