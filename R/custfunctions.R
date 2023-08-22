#' ETC: Function to run equivalence testing on single correlation = 1
#'
#' @param dat Matrix object representing the data, subjects as rows, variables as columns.
#' @param varA Column number corresponding to one variable of the pair being tested for redundancy.
#' @param varB Column number corresponding to the other variable of the pair being tested for redundancy.
#' @param ei Single value indicating the width of the equivalence interval, where 1-\code{"ei"} represents the lower bound of this interval and 1 represents the upper bound.
#' @param alpha Single value representing the desired alpha level.
#' @details For a single pair of variables, this function runs equivalence testing on a single zero-order correlation. Specifically, the test evaluates whether the variables' correlation is close enough to 1 as to declare the variables effectively redundant. The inspiration for this function came from work by Goertzen & Cribbie (2010).
#'
#' @references
#'
#' Goertzen, J. R., & Cribbie, R. A. (2010). Detecting a lack of association: An equivalence testing approach. British Journal of Mathematical and Statistical Psychology, 63(3), 527-537.
#'
#' Some older refs too
#'
#' @return A list that includes a z-statistic, p-value, and redundancy decision, in addition to select input parameters for convenience.
#'
#' @export
#' @importFrom stats cor cov pnorm qnorm
#' @examples
#' \dontrun{
#' #### 1. Example with bfi dataset from psych package
#'
#'    ## Clean dataset to get complete data only (n = 2436, and 25 items)
#'    # Needs psych package loaded in
#'    library("psych")
#'    data(bfi)
#'    bfisub = bfi[,1:25]
#'    # Remove any rows with missing data
#'    bfisub <- na.omit(bfisub)
#'
#'    ## Pick alpha level and equivalence interval
#'    alpha = .05
#'    ei = .1
#'
#'    ## Run ETC test on single pair of items (N1 and N2)
#'    ETC <- ETC(dat=bfisub, varA=16, varB=17, alpha=alpha, ei=ei)
#' }
#'
ETC <- function(dat, varA, varB, ei, alpha = 0.05) {
  x1 <- dat[, varA]
  x2 <- dat[, varB]
  corx1x2 <- cor(x1, x2)
  n <- nrow(dat)

  # Run a t-test procedure with Fisher's z transformation
  zei <- log((1 + (1 - ei))/(1 - (1 - ei)))/2
  zcorx1x2 <- log((1 + corx1x2)/(1 - corx1x2))/2
  equivt_fz <- (zei - zcorx1x2)/(1/sqrt(n - 3)) # Designed for testing null that rho < rhotilde
  pvalue_fz <- pnorm(equivt_fz)
  ifelse(pvalue_fz <= alpha,
         decis_fz <- "Redundant",
         decis_fz <- "Not Redundant")

  # Summarize
  stats_fz <- matrix(c(round(corx1x2,3), round(1-ei,2), ei, round(equivt_fz,2), round(pvalue_fz,3)),1,5)
  colnames(stats_fz) <- c("r", "rhotilde", "EI", "z", "p-value")
  output <- list(stats_fz, ei, decis_fz, corx1x2, "ETC")
  names(output) <- c("Results",
                     "EI",
                     "Decision",
                     paste0("r",varA,varB),
                     "Method")
  return(output)
}

#' ETT: Function to run equivalence testing on difference between two dependent correlations
#'
#' @param dat Matrix object representing the data, subjects as rows, variables as columns.
#' @param varA Column number corresponding to one variable of the pair being tested for redundancy.
#' @param varB Column number corresponding to the other variable of the pair being tested for redundancy.
#' @param ei Single value indicating half the width of the equivalence interval, where -\code{"ei"} represents the lower bound of this interval and \code{"ei"} represents the upper bound (i.e., 0 is the midpoint).
#' @param alpha Single value representing the desired alpha level.
#' @param threshold Single value between 0 and 1, representing the minimum proportion of rejections required for the two variables to be declared redundant.
#' @param corMin Single value between 0 and 1, representing the minimum zero-order correlation required (between \code{"varA"} and \code{"varB"}) for the two variables to be declared redundant.
#' @details For a single pair of variables, this function runs equivalence testing on the topology of the variables - in other words, the difference between the correlations of \code{"varA"} and \code{"varB"} with each of the other variables in the dataset. Specifically, the test evaluates whether this difference is within an user-defined interval around zero such that the variables can be declared redundant. Note that multiple tests are conducted for a given variable pair when the number of variables in the dataset exceeds three. The inspiration for this function came from work by Counsell & Cribbie (2015).
#'
#' @references
#'
#' Counsell, A., & Cribbie, R. A. (2015). Equivalence tests for comparing correlation and regression coefficients. British Journal of Mathematical and Statistical Psychology, 68(2), 292-309.
#'
#' Some older refs too
#'
#' @return A list that includes results for all tests including p-values, as well as the proportion of tests resulting in rejection of the null hypothesis, and the redundancy decision, in addition to select input parameters for convenience.
#'
#' @export
#' @importFrom stats cor cov pnorm qnorm
#' @examples
#' \dontrun{
#' #### 1. Example with bfi dataset from psych package
#'
#'    ## Clean dataset to get complete data only (n = 2436, and 25 items)
#'    # Needs psych package loaded in
#'    library("psych")
#'    data(bfi)
#'    bfisub = bfi[,1:25]
#'    # Remove any rows with missing data
#'    bfisub <- na.omit(bfisub)
#'
#'    ## Pick alpha level, equivalence interval, proportion rejected threshold (ETT),
#'    ## and minimum zero-order correlation (ETT)
#'    alpha = .05
#'    ei = .1
#'    threshold = .75
#'    corMin = .5
#'
#'    ## Run ETC test on single pair of items (N1 and N2)
#'    ETT <- ETT(dat=bfisub, varA=16, varB=17, alpha=alpha, ei=ei,
#'      threshold=threshold, corMin=corMin)
#' }
#'
ETT <- function(dat, varA, varB, ei, alpha = 0.05, threshold = 0.75, corMin = 0.50) {
  # Define vector containing all variable indices not in the pair of interest
  nvar <- dim(cor(dat))[1]
  vthird <- 1:nvar
  pairdrop <- c(varA, varB)
  vthird <- vthird[-pairdrop]
  # Define initial results matrix
  results.mat <- matrix(NA,0,4)
  # Loop the function over all "third" variables
  for(varC in vthird){
    x1 <- dat[, varA]
    x2 <- dat[, varB]
    x3 <- dat[, varC]
    r13 <- cor(x1, x3) # First correlation
    r23 <- cor(x2, x3) # Second correlation
    r12 <- cor(x1, x2) # Pair-of-interest correlation
    n <- nrow(dat)

    # Compute the determinant of the 3x3 correlation matrix for varA, varB, varC
    detR <- (1 - r13^2 - r23^2 - r12^2) + (2 * r13 * r23 * r12)
    # Compute std. error term
    serr <- 1/sqrt(((n - 1) * (1 + r12))/((2 * ((n - 1)/(n - 3)) * detR) +
                                            (((r13 + r23)^2)/4) * ((1 - r12)^3)))
    # Compute p-value
    p1 <- pnorm((abs(r13 - r23) - ei) * (1/serr))
    p2 <- pnorm((-abs(r13 - r23) - ei) * (1/serr))
    pval <- p1 - p2
    # Compute CIs
    upper <- (r13 - r23) + qnorm(alpha) * serr
    lower <- (r13 - r23) - qnorm(alpha) * serr
    if (lower < upper) {
      lower2 <- lower
      upper2 <- upper
    }
    if (lower > upper) {
      lower2 <- upper
      upper2 <- lower
    }
    CI <- c(lower2, upper2)

    results.single <- matrix(c(varC, round(pval,3), round(CI[1],3), round(CI[2],3)),
                             1, 4)
    results.mat <- rbind(results.mat,results.single)

  }

  # Name results matrix columns
  colnames(results.mat) <- c("VarC", "p-value", "CI(L)","CI(U)")

  # Get proportion of comparisons rejected
  prop.rej <- sum(results.mat[,2] <= alpha)/(nvar-2)
  # Get redundancy decision
  redundancy.decision <- if ((prop.rej >= threshold)&(r12>corMin)) {
    "Redundant"
  } else {
    "Not Redundant"
  }

  # Generate output
  output <- list(results.mat, corMin, ei, prop.rej, r12, redundancy.decision, "ETT")
  names(output) <- c("Results",
                     "corMin",
                     "EI",
                     "PropRej",
                     paste0("r",varA,varB),
                     "Decision",
                     "Method")
  return(output)

}









