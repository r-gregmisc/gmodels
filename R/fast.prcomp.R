#' Efficient computation of principal components and singular value
#' decompositions.
#' 
#' The standard \code{\link[stats]{prcomp}} and \code{\link{svd}} function are
#' very inefficient for wide matrixes. \code{fast.prcomp} and \code{fast.svd}
#' are modified versions which are efficient even for matrixes that are very
#' wide.
#' 
#' The current implementation of the function \code{\link{svd}} in S-Plus and R
#' is much slower when operating on a matrix with a large number of columns
#' than on the transpose of this matrix, which has a large number of rows. As a
#' consequence, \code{\link[stats]{prcomp}}, which uses \code{\link{svd}}, is
#' also very slow when applied to matrixes with a large number of rows.
#' 
#' The simple solution is to use \code{\link{La.svd}} instead of
#' \code{\link{svd}}.  A suitable patch to \code{\link[stats]{prcomp}} has been
#' submitted.  In the mean time, the function \code{fast.prcomp} has been
#' provided as a short-term work-around.
#'  
#' \describe{ 
#'  \item{list("fast.prcomp")}{is a modified versiom of
#'    \code{\link[stats]{prcomp}} that calls \code{\link{La.svd}} instead of
#'      \code{\link{svd}} } 
#'    \item{list("fast.svd")}{is simply a wrapper around
#'      \code{\link{La.svd}}.  } 
#' }
#' 
#' @aliases fast.prcomp fast.svd
#' 
#' @param x data matrix
#' @inheritParams stats::prcomp
#' @inheritParams base::svd
#' @return See the documetation for \code{\link[stats]{prcomp}} or
#' \code{\link{svd}} .
#' @author Modifications by Gregory R. Warnes \email{greg@@warnes.net}
#' @seealso \code{\link[stats]{prcomp}}, \code{\link{svd}},
#' \code{\link{La.svd}}
#' @keywords multivariate algebra array
#' @examples
#' 
#' 
#'   # create test matrix
#'   set.seed(4943546)
#'   nr <- 50
#'   nc <- 2000
#'   x  <- matrix( rnorm( nr*nc), nrow=nr, ncol=nc )
#'   tx <- t(x)
#' 
#'   # SVD directly on matrix is SLOW:
#'   system.time( val.x <- svd(x)$u )
#' 
#'   # SVD on t(matrix) is FAST:
#'   system.time( val.tx <- svd(tx)$v )
#' 
#'   # and the results are equivalent:
#'   max( abs(val.x) - abs(val.tx) )
#' 
#'   # Time gap dissapears using fast.svd:
#'   system.time( val.x <- fast.svd(x)$u )
#'   system.time( val.tx <- fast.svd(tx)$v )
#'   max( abs(val.x) - abs(val.tx) )
#' 
#' 
#'   library(stats)
#' 
#'   # prcomp directly on matrix is SLOW:
#'   system.time( pr.x <- prcomp(x) )
#' 
#'   # prcomp.fast is much faster
#'   system.time( fast.pr.x <- fast.prcomp(x) )
#' 
#'   # and the results are equivalent
#'   max( pr.x$sdev - fast.pr.x$sdev )
#'   max( abs(pr.x$rotation[,1:49]) - abs(fast.pr.x$rotation[,1:49]) )
#'   max( abs(pr.x$x) - abs(fast.pr.x$x)  )
#' 
#'   # (except for the last and least significant component):
#'   max( abs(pr.x$rotation[,50]) - abs(fast.pr.x$rotation[,50]) )
#' 
#' @export
fast.prcomp <- function (x, retx = TRUE, center = TRUE, scale. = FALSE,
                         tol = NULL)
{
  x <- as.matrix(x)
  x <- scale(x, center = center, scale = scale.)
  s <- La.svd(x, nu = 0)
  if (!is.null(tol)) {
    rank <- sum(s$d > (s$d[1] * tol))
    if (rank < ncol(x))
      s$vt <- s$vt[, 1:rank, drop = FALSE]
  }
  s$d <- s$d/sqrt(max(1, nrow(x) - 1))
  
  dimnames(s$vt) <- list(paste("PC", seq(len = nrow(s$vt)), sep = ""),
                         colnames(x) )
  r <- list(sdev = s$d, rotation = t(s$vt) )
  if (retx)
    r$x <- x %*% t(s$vt)
  class(r) <- "prcomp"
  r
}

#' @export
fast.svd <- function( x, nu = min(n, p), nv = min(n, p), ...)
{
  x <- as.matrix(x)
  dx <- dim(x)
  n <- dx[1]
  p <- dx[2]
  
  retval <- La.svd(x, nu=nu, nv=nv,  ... )
  retval$v <- t(retval$vt)
  retval$vt <- NULL
  retval
}

