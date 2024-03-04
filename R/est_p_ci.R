#' Display Confidence interval and p-value for one model term
#' 
#' @param model model object
#' @param term model term 
#' @param mult  scale (multiply) the parameter by this factor 
#' @param digits number of significant digits to display
#' @param ... optional arguments
#' @export
est_p_ci <- function(model, term, mult=1, digits=2, ...)
  UseMethod("est_p_ci")

#' @exportS3Method gmodels::est_p_ci
est_p_ci.lm <- function(model, term, mult=1, digits=2, ...)
{
  info <- ci(model)
  if(is.character(term) && !(term %in% rownames(info)))
    stop(term, " is not a coefficient in model.")
  info <- info[term,]
  info.ci <- trimws(format( round(mult * info[1:3], digits=digits) ))
  if(mult < 0)
    names(info.ci) <- rev(names(info.ci))
  paste("Est=", info.ci[1],
        " ",
        "p=",format.pval(info['p-value'], digits=digits),
        " ",
        "95% CI: ",
        info.ci['CI lower'],
        " to ",
        info.ci['CI upper'],
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

  info.ci <- trimws(format( round(mult * model[term, c("lower CI", "upper CI")],
                           digits=digits) )  )
  if(mult < 0)
    names(info.ci) <- rev(names(info.ci))
  paste("Est=", info.ci[1],
        " ",
        "p=",format.pval(model[term, 'Pr(>|t|)'], digits=digits),
        " ",
        "95% CI: ",
        info.ci["lower CI"],
        " to ",
        info.ci["upper CI"],
        sep=""
  )
}
