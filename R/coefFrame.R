#' Return model parameters in a data frame
#' 
#' Fits a model to each subgroup defined by \code{by}, then returns a data
#' frame with one row for each fit and one column for each parameter.
#' 
#' 
#' @param mod a model formula, to be passed to by \code{fitfun}.
#' @param data a data frame, row subsets of which will be used as the
#' \code{data} argument to \code{fitfun}.
#' @param by names of columns in \code{x} that will be used to define the
#' subgroups.
#' @param fit.on a logical vector indicating which rows of \code{x} are to be
#' used to fit the model (like the \code{subset} argument in a lot of other
#' functions). Can be given in terms of variables in \code{x}
#' @param fitfun a model fitting function (e.g. lm, nls). More specifically, a
#' function that expects at least a formula object (as the first argument) and
#' a data.frame object (passed as an argument named \code{data}) and returns a
#' model object for which a \code{coef} method has been defined (e.g. coef.lm,
#' coef.nls) to extract fit values of model parameters.
#' @param keep.unused.levels Include rows in output for all unique values of
#' \code{by}, even those which were excluded by \code{fit.on}. The default
#' value \code{TRUE} should be left alone if you are going to go on to pass the
#' result to \code{backFit}.
#' @param byvar.sep passed to frameApply, used to form the subsets of the data.
#' @param ... other arguments to pass to \code{fitfun}.
#' @return a data frame with a row for each unique row of \code{x[by]}, and
#' column for each model paramter, as well as columns specified in \code{by}.
#' @author Jim Rogers \email{james.a.rogers@@pfizer.com}
#' @keywords models
#' @examples
#' 
#' # load example data
#' library(gtools)
#' data(ELISA)
#' 
#' # Coefficients for four parameter logistic fits:
#' coefFrame(log(Signal) ~ SSfpl(log(Concentration), A, B, xmid, scal),
#'            data = ELISA, fitfun = nls,
#'            by = c("PlateDay", "Read"),
#'            fit.on = Description == "Standard" & Concentration != 0)
#' 
#' # Coefficients for linear fits:
#' coefFrame(log(Signal) ~ log(Concentration), 
#'            data = ELISA, fitfun = lm, 
#'            by = c("PlateDay", "Read"),
#'            fit.on = Description == "Standard" & Concentration != 0 )
#' 
#' # Example passing arguments to fitfun, and example of
#' # error handling during model fitting:
#' ELISA$Signal[1] <- NA
#' coefFrame(log(Signal) ~ log(Concentration), 
#'            data = ELISA, fitfun = lm, na.action = na.fail,
#'            by = c("PlateDay", "Read"),
#'            fit.on = Description == "Standard" & Concentration != 0 )
#' 
#' 
#' 
#' 
#' @importFrom gdata frameApply
#' @importFrom stats coef
#'
#' @export
coefFrame <-  function (
    mod, 
    data, 
    by = NULL, 
    fit.on = TRUE, 
    fitfun, 
    keep.unused.levels = TRUE, 
    byvar.sep = "\001" , 
    ...
    ) 
{
    fit.on <- eval(substitute(fit.on), data, parent.frame())
    out <- frameApply(data, on = intersect(all.vars(mod), names(data)), 
        by = by, subset = fit.on, byvar.sep = byvar.sep, fun = function(sub.dat, 
            ...) {
            fit <- try(fitfun(mod, data = sub.dat, ...), silent = TRUE)
            if (inherits(fit, "try-error")) 
                return(fit)
            outi <- coef(fit)
            outi
        }, ...)
    if (keep.unused.levels) 
        out <- unique(merge(data[by], out, all.x = TRUE))
    out
}


