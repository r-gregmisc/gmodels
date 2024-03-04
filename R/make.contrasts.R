#' Construct a User-Specified Contrast Matrix
#' 
#' This function converts human-readable contrasts into the form that R
#' requires for computation.
#' 
#' Specifying a contrast row of the form `c(1,0,0,-1)` creates a contrast
#' that will compare the mean of the first group with the mean of the fourth
#' group.
#' 
#' @param contr vector or matrix specifying contrasts (one per row).
#' @param how.many dimensions of the desired contrast matrix.  This must equal
#' the number of levels of the target factor variable.
#' @return `make.contrasts` returns a matrix with dimensions
#' (`how.many`, `how.many`) containing the specified contrasts
#' augmented (if necessary) with orthogonal "filler" contrasts.
#' 
#' This matrix can then be used as the argument to [contrasts()] or
#' to the `contrasts` argument of model functions (eg, [lm()]).
#' @author Gregory R. Warnes \email{greg@@warnes.net}
#' @seealso 
#'  * [stats::lm()], [stats::contrasts()], [stats::contr.treatment()], 
#'    [stats::contr.poly()], 
#'  * Computation and testing of General Linear Hypothesis: [glh.test()], 
#'  * Computation and testing of estimable functions of model coefficients:
#'    [estimable()], 
#'  * Estimate and Test Contrasts for a previously fit linear model: 
#'    [fit.contrast()]
#'    
#' @keywords models regression
#' 
#' @examples
#' 
#' set.seed(4684)
#' y <- rnorm(100)
#' x.true <- rnorm(100, mean=y, sd=0.25)
#' x <-  factor(cut(x.true,c(-4,-1.5,0,1.5,4)))
#' reg <- lm(y ~ x)
#' summary(reg)
#' 
#' # Mirror default treatment contrasts
#' test <- make.contrasts(rbind( c(-1,1,0,0), c(-1,0,1,0), c(-1,0,0,1) ))
#' lm( y ~ x, contrasts=list(x = test ))
#' 
#' # Specify some more complicated contrasts
#' #   - mean of 1st group vs mean of 4th group
#' #   - mean of 1st and 2nd groups vs mean of 3rd and 4th groups
#' #   - mean of 1st group vs mean of 2nd, 3rd and 4th groups
#' cmat <- rbind( "1 vs 4"    =c(-1, 0, 0, 1),
#'                "1+2 vs 3+4"=c(-1/2,-1/2, 1/2, 1/2),
#'                "1 vs 2+3+4"=c(-3/3, 1/3, 1/3, 1/3))
#' 
#' summary(lm( y ~ x, contrasts=list(x=make.contrasts(cmat) )))
#' # or
#' contrasts(x) <- make.contrasts(cmat)
#' summary(lm( y ~ x ) )
#' 
#' # or use contrasts.lm
#' reg <- lm(y ~ x)
#' fit.contrast( reg, "x", cmat )
#' 
#' # compare with values computed directly using group means
#' gm <- sapply(split(y,x),mean)
#' gm %*% t(cmat)
#' 
#' 
#' #
#' # Example for Analysis of Variance
#' #
#' 
#' set.seed(03215)
#' Genotype <- sample(c("WT","KO"), 1000, replace=TRUE)
#' Time <- factor(sample(1:3, 1000, replace=TRUE))
#' data <- data.frame(y, Genotype, Time)
#' y <- rnorm(1000)
#' 
#' data <- data.frame(y, Genotype, as.factor(Time))
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
#' model <- model <- aov( y ~ Genotype + Time + Genotype:Time, data=data,
#'                       contrasts=list(Genotype=make.contrasts(cm.G),
#'                                      Time=make.contrasts(cm.T) )
#'                       )
#' 
#' summary(model, split=list( Genotype=list( "KO vs WT"=1 ),
#'                            Time = list( "1 vs 2" = 1,
#'                                         "2 vs 3" = 2 ) ) )
#' 
#' @importFrom MASS ginv
#' 
#' @export
make.contrasts <-  function (
    contr, 
    how.many=ncol(contr)
    ) 
{
  if(!is.matrix(contr))
    contr <- matrix(contr,ncol=length(contr))

  if(nrow(contr)+1 > how.many)
    stop("Too many contrasts specified. Must be less than the number of factor levels (columns).")
  
  value <- as.matrix(ginv(contr))  # requires library(MASS)
  if (nrow(value) != how.many) 
    stop("wrong number of contrast matrix rows")
  n1 <- if (missing(how.many)) 
    how.many - 1
  else how.many
  nc <- ncol(value)
  if (nc < n1) {
    cm <- qr(cbind(1, value))
    if (cm$rank != nc + 1) 
      stop("singular contrast matrix")
    cm <- qr.qy(cm, diag(how.many))[, 2:how.many, drop=FALSE]
    cm[, 1:nc] <- value
  }
  else cm <- value[, 1:n1, drop = FALSE]

  colnames(cm) <- paste( "C", 1:ncol(cm), sep="")
  rownames(cm) <- paste( "V", 1:nrow(cm), sep="")
  
  if(!is.null(rownames(contr)))
    {
      namelist <- rownames(contr)
      colnames(cm)[1:length(namelist)] <- namelist
    }

  if(!is.null(colnames(contr)))
    rownames(cm) <- colnames(contr)

  cm
}
