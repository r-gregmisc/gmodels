`estimable.mlm` <-
  function (obj, cm, beta0, conf.int=NULL,  show.beta0, ...)
{
  coef <- coef(object)
  ny <- ncol(coef)
  effects <- object$effects
  resid <- object$residuals
  fitted <- object$fitted
  ynames <- colnames(coef)
  if (is.null(ynames)) {
    lhs <- object$terms[[2]]
    if (mode(lhs) == "call" && lhs[[1]] == "cbind") 
      ynames <- as.character(lhs)[-1]
    else ynames <- paste("Y", seq(ny), sep = "")
  }
  value <- vector("list", ny)
  names(value) <- paste("Response", ynames)
  cl <- oldClass(object)
  class(object) <- cl[match("mlm", cl):length(cl)][-1]
  for (i in seq(ny)) {
    object$coefficients <- coef[, i]
    object$residuals <- resid[, i]
    object$fitted.values <- fitted[, i]
    object$effects <- effects[, i]
    object$call$formula[[2]] <- object$terms[[2]] <- as.name(ynames[i])
    value[[i]] <- estimable(obj, cm, beta0, conf.int=NULL,  show.beta0, ...)
  }
  class(value) <- "listof"
  value
}
