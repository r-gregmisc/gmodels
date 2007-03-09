# est.lmer.R
# generate estimable output for lmer objects using mcmcsamp()
# Randall Johnson
# Laboratory of Genomic Diversity at NCI Frederick
# SAIC Frederick, Inc
# Created April 25, 2006

est.lmer <- function(obj, cm, beta0, conf.int, show.beta0, n.sim)
{
  if(!require(coda, quietly=TRUE))
    stop("coda package required when sim.lmer == TRUE")

  samp <- mcmcsamp(obj, n.sim)
  samp.summ <- summary(samp)

  if(is.null(dim(cm)))
    n <- length(cm)
  else
    n <- dim(cm)[2]
                               # drop extra information on end
  samp.cm <- as.matrix(samp)[, 1:n] %*% t(cm)

  # calculate requested statistics
  est <- drop(cm %*% samp.summ$statistics[1:n,1])
  stderr <- sd(samp.cm)
  
  pval <- sapply(1:length(beta0),
                 function(i){percentile(beta0[i], samp.cm[,i])})
  pval <- ifelse(pval <= .5, 2*pval, 2*(1-pval))

  if(is.null(conf.int))
  {
    lower.ci <- NULL
    upper.ci <- NULL
  }else{
    alpha <- 1-conf.int
    samp.ci <- sapply(1:length(beta0),
                      function(i)
                        {
                          quantile(samp.cm[,i], probs=c(alpha/2, 1-alpha/2))
                        }
                      )

    lower.ci <- samp.ci[1,]
    upper.ci <- samp.ci[2,]
  }

  # return results
  if(!show.beta0)
    beta0 <- NULL
  
  samp.stats <- cbind('beta0' = beta0,
                      'Estimate' = est,
                      'Std. Error' = stderr,
                      'p value' = pval,
                      'Lower.CI' = lower.ci,
                      'Upper.CI' = upper.ci)

  row.names(samp.stats) <- paste('(', apply(cm, 1, paste, collapse=" "),
                                 ')', sep='')
    
  return(samp.stats)
}

percentile <- function(x, distn)
{
  n <- length(distn)
    
  return(findInterval(x, distn[order(distn)]) / n)
}
