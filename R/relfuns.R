# Model-based functions to take into account error (reliability)
# 2023-07-11 Carl F. Falk
# 2023-07-17 CFF - adding documentation
# 2023-08-22 CFF - change some defaults for model fitting, bugfix

#' Create (b)lavaan syntax to fit two alternative models
#'
#' One model is saturated; other assumes two items measure same construct
#' @param v1 name of variable 1
#' @param v2 name of variable 2
#' @param dat dataset; should be a data frame with column labels. Should also
#' only contain the variables of interest.
#' @param method which type of model should we estimate?
#'  lv = latent variable model to test topology. pconst = implement
#'    proportionality constraints on topology.
#'  sat = just does saturated model (e.g., in case one wants only to check for
#'    some constraint on the correlation b/w variables)
#' @param warmstart use sample means and covariances to help w/ starting values?
#' @param corconst if there should be any constraint on the model-implied
#'   correlation between the two variables, put the value here. Defaults to no
#'   such constraint (i.e., test only topology)
#'
#' @export
# TODO: is it possible to parameterize corconst in terms of standardized loadings?
# TODO: input checking so that variable names don't conflict with parameter labels
syn_redmod <- function(dat, v1, v2, method=c("lv","pconst","sat"),
                       warmstart=FALSE, corconst=NULL){

  # input checking
  method <- match.arg(method)

  stopifnot(is.data.frame(dat))
  stopifnot(!is.null(colnames(dat)))
  stopifnot(v1 %in% colnames(dat))
  stopifnot(v2 %in% colnames(dat))
  if(!is.null(corconst)) stopifnot(corconst > 0 & corconst < 1)

  if(warmstart){
    S <- round(cov(dat),3) # sample covariance matrix; can use for starting values
    m <- round(colMeans(dat),3) # sample means; can use for starting values
  }

  vn <- colnames(dat) # all variable names
  vo <- vn[!vn %in% c(v1,v2)] # other variables not in v1, v2

  ## Set up stuff that all models will have
  # blavaan will expect intercepts, so include them
  mod <- "# Intercepts"
  if(warmstart){
    mod <- paste(mod, paste0(vn, " ~ start(", m, ")*1", collapse="\n"), sep="\n")
  } else {
    mod <- paste(mod, paste0(vn, " ~ 1", collapse="\n"), sep="\n")
  }

  # (error) variances
  mod <- paste(mod, "# (error) variances", sep="\n")
  if(warmstart){
    mod <- paste(mod, paste0(vn, " ~~ start(", diag(S), ")*", vn, " + e", vn, "*", vn,  collapse="\n"), sep="\n")
  } else {
    mod <- paste(mod, paste0(vn, " ~~ e", vn, "*", vn, collapse="\n"), sep="\n")
  }

  # covariance b/w other variables
  mod <- paste(mod, "# Other covariances", sep="\n")
  for(j in 1:(length(vo)-1)){
    for(jj in ((j+1):length(vo))){
      if(warmstart){
        mod <- paste(mod, paste0(vo[j], " ~~ start(",S[vo[j],vo[jj]], ")*", vo[jj]), sep="\n")
      } else {
        mod <- paste(mod, paste0(vo[j], " ~~ ", vo[jj]), sep="\n")
      }
    }
  }

  ## create syntax for latent variable approach
  if(method=="lv"){

    modsyn <- paste(mod, "# Latent variable", sep="\n")
    modsyn <- paste(modsyn, paste0("F1 =~ l1*", v1, " + l2*", v2), sep="\n")
    modsyn <- paste(modsyn, "F1 ~~ vF1*F1", sep="\n")

    # covariance b/w other variables and F1
    modsyn <- paste(modsyn, "# Cov b/w F1 and other variables", sep="\n")
    for(j in 1:length(vo)){
      modsyn <- paste(modsyn, paste0("F1 ~~ ", vo[j]), sep="\n")
    }

  }

  ## Saturated model
  # Between v1 and v2
  mod <- paste(mod, "# Cov b/w variables of interest", sep="\n")
  if(warmstart){
    mod <- paste(mod, paste0(v1, " ~~ start(", S[v1,v2] , ")*", v2, "+ cov*", v2), sep="\n")
  } else {
    mod <- paste(mod, paste0(v1, " ~~ ", v2, "+ cov*", v2), sep="\n")
  }

  # Between v1, v2, and other variables
  mod <- paste(mod, "# Cov b/w variables of interest and other variables", sep="\n")
  for(j in 1:length(vo)){
    # includes parameter labels for constraints
    if(warmstart){
      mod <- paste(mod, paste0(v1, " ~~ c1",j, "*", vo[j], " + start(", S[v1,vo[j]] ,")*", vo[j]), sep="\n")
      mod <- paste(mod, paste0(v2, " ~~ c2",j, "*", vo[j], " + start(", S[v2,vo[j]] ,")*", vo[j]), sep="\n")
    } else {
      mod <- paste(mod, paste0(v1, " ~~ c1",j, "*", vo[j]), sep="\n")
      mod <- paste(mod, paste0(v2, " ~~ c2",j, "*", vo[j]), sep="\n")
    }
  }

  ## create syntax for proportionality constraint approach
  if (method=="pconst"){

    # constraints
    modsyn <- paste(mod, "# Proportionality constraints", sep="\n")
    for(j in 2:length(vo)){
      # more than one way to do this; not all options yet provided
      modsyn <- paste(modsyn, paste0("c11/c21 == c1", j, "/c2", j), sep="\n")
    }
  }

  if(method=="sat"){
    modsyn <- mod
  }

  if(!is.null(corconst)){

    modsyn <- paste(modsyn, "# Correlation constraints", sep="\n")

    # correlation b/w variables of interest
    if(method=="lv"){

      # define correlation
      #modsyn <- paste(modsyn, paste0("r12 := l1*l2*vF1 / sqrt((l1^2*vF1 + e",v1,")*(l2^2*vF1 + e",v2,"))"), sep="\n")
      # lavaan does not like vF1 in definition if it is fixed
      modsyn <- paste(modsyn, paste0("r12 := l1*l2 / sqrt((l1^2 + e",v1,")*(l2^2 + e",v2,"))"), sep="\n")

      #r12 := c12/sqrt(v1*v2)
    } else if (method=="pconst"){

      # define correlation
      modsyn <- paste(modsyn, paste0("r12 := cov/sqrt(e", v1, "*e", v2, ")"), sep="\n")

    }

    # apply inequality constraint
    #modsyn <- paste(modsyn, paste0(corconst, "< abs(r12)"), sep="\n")
    modsyn <- paste(modsyn, paste0("abs(r12) > ", corconst), sep="\n")

  }

  out <- list(mod = modsyn,
              modsat = mod,
              v1 = v1,
              v2 = v2,
              method = method,
              warmstart = warmstart,
              corconst = corconst)

  return(out)
}

#' Fits two alternative (b)lavaan models for purposes of comparing item overlap
#'
#' @param syn syntax of two models, as created by \code{syn_redmod}
#' @param dat dataset; should be a data frame with column labels. Should also
#'  only contain the variables of interest.
#' @param software which software package to use? Note that use of
#'  \code{blavaan} is still experimental.
#' @param likelihood for \code{lavaan}, which likelihood to use? (Used to
#'  override default behavior)
#' @param std.lv for \code{lavaan} and \code{blavaan}, whether to standardize
#'  latent variables (Used to override default behavior)
#' @param auto.fix.first for \code{lavaan} and \code{blabvaan}, whether to fix
#'  loading for first item of each latent variable to one (Used to override
#'  default behavior)
#' @param h1 See documentation for \code{lavOptions}.
#' @param baseline See documentation for \code{lavOptions}.
#' @param ... Additional arguments passed to \code{lavaan} or \code{blavaan}
#' @import lavaan
#' @import blavaan
#' @export
# TODO: convergence of blavaan models
fit_redmod <- function(syn, dat, software = c("lavaan","blavaan"),
                       likelihood="wishart", std.lv=TRUE, auto.fix.first=FALSE,
                       h1=FALSE, baseline=FALSE,
                       ...){

  software <- match.arg(software)

  if(software =="lavaan"){
    fitsat <- lavaan(syn$modsat, dat, likelihood=likelihood, h1=h1,
                     baseline=baseline, ...)
    fit <- lavaan(syn$mod, dat, likelihood=likelihood, std.lv=std.lv,
                  auto.fix.first=auto.fix.first, h1=h1, baseline=baseline, ...)
    out <- list(
      Fit = fit,
      Fitsat = fitsat,
      software = software,
      conv = lavInspect(fit, "converged"),
      convsat = lavInspect(fitsat, "converged")
    )
  } else if (software =="blavaan"){
    fitsat <- blavaan(syn$modsat, dat, ...)
    fit <- blavaan(syn$mod, dat, std.lv=std.lv, auto.fix.first=auto.fix.first,
                   h1=h1, baseline=baseline, ...)
    out <- list(
      Fit = fit,
      Fitsat = fitsat,
      software = software
    )
  }

  return(out)

}

#' Does model selection
#'
#' @param fits two fitted models, as created by \code{syn_redmod} then
#'  \code{fit_redmod}
#' @param index which fit index to examine? Options depend on whether
#'  \code{lavaan} or \code{blavaan} was used to estimate the models.
#' @param thresh threshold for determining whether the saturated model wins.
#'  Positive values result in a stronger preference for parsimony (i.e., items
#'  are redundant).
#' @importFrom stats AIC BIC
#' @importFrom blavaan blavCompare
#' @importFrom semTools moreFitIndices
#' @export
sel_redmod <- function(fits, index=c("aic","bic","sic","waic","loo"), thresh=0){
  index <- match.arg(index)

  if(fits$software=="lavaan"){
    stopifnot(index %in% c("aic","bic","sic"))
    if(index == "aic"){
      idx <- AIC(fits$Fit)
      idxsat <- AIC(fits$Fitsat)
    } else if (index == "bic"){
      idx <- BIC(fits$Fit)
      idxsat <- BIC(fits$Fitsat)
    } else if (index == "sic"){
      idx <- moreFitIndices(fits$Fit, "sic")
      idxsat <- moreFitIndices(fits$Fitsat, "sic")
    }

  } else if (fits$software=="blavaan"){
    stopifnot(index %in% c("waic","loo"))
    idxs <- blavCompare(fits$Fit, fits$Fitsat)
    if(index=="loo"){
      idx <- idxs$loo[[1]]["looic","Estimate"]
      idxsat <- idxs$loo[[2]]["looic","Estimate"]
    } else if (index=="waic"){
      idx <- idxs$waic[[1]]["waic","Estimate"]
      idxsat <- idxs$waic[[2]]["waic","Estimate"]
    }
  }

  # does the saturated model win?
  satpreferred <- ifelse((idxsat - idx) < thresh, TRUE, FALSE)
  
  # if lavaan used and either model not converged, declare NA (i.e., inconclusive)
  if(fits$software=="lavaan" & !fits$conv | !fits$convsat){satpreferred <- NA}

  out <- list(satpreferred = satpreferred,
              idx = idx,
              idxsat = idxsat,
              thresh = thresh)
  return(out)
}


#' Wrapper to compare two variables at a time using model-based approach
#'
#' @param dat dataset; should be a data frame with column labels. Should also
#'  only contain the variables of interest.
#' @param v1 name of variable 1
#' @param v2 name of variable 2
#' @param method which type of model should we estimate?
#'  lv = latent variable model to test topology.
#'  pconst = implement proportionality constraints on topology.
#'  sat = just does saturated model (e.g., in case one wants only to check for
#'  some constraint on the correlation b/w variables).
#' @param software which software package to use? Note that use of
#'  \code{blavaan} is experimental.
#' @param index which fit index to examine? Options depend on whether
#'  \code{lavaan} or \code{blavaan} was used to estimate the models.
#' @param thresh threshold for determining whether the saturated model wins.
#'  Positive values result in a stronger preference for parsimony (i.e., items
#'  are redundant).
#' @param warmstart use sample means and covariances to help w/ starting values?
#' @param corconst if there should be any constraint on the model-implied
#'  correlation between the two variables, put the value here. Defaults to no
#'  such constraint (i.e., test only topology)
#' @param ... Additional arguments passed to \code{lavaan} or \code{blavaan}
#'
#' @details
#' This function attempts to automate the procedure used by Starr and Falk
#' (see references) for detecting redundant items using a model-based approach.
#' A saturated model is compared to one in which the two items in question
#' are represented as either a single latent variable, or in which the items'
#' covariances are subject to proportionality constraints. If the latter model
#' fits better than the saturated model (according to information criteria, 
#' and optionally some threshold on how much improved model fit must be), the
#' items may be flagged as possibly redundant.
#' 
#' Note that although the function has been evaluated in simulations, it is very
#' new and any issues in model fitting should be reported to the authors. Note
#' that use of blavaan has not yet been evaluated in simulations.
#' @references Starr, J., \& Falk, C.F. On the testing of equivalent items:
#' Perfect correlations and correlational topology. Preprint:
#' https://doi.org/10.31234/osf.io/vhgfk
#' @export
#' @examples
#' \dontrun{
#' # Example with bfi dataset from psych package
#' library(psych)
#' bfisub <- bfi[,1:25]
#'
#' # remove any rows with missing data
#' bfisub <- na.omit(bfisub)
#'
#' # Check for redundancy for two variables
#' # Use latent variable representation, BIC for model selection,
#' # and saturated model needs to have BIC 3 lower than alternative model
#' # to be selected.
#' bic.res1 <- redmod_pair(bfisub,"A2","A3",method="lv",index="bic", thresh=3,
#'                         optim.method="BFGS")
#' 
#' # Additionally include threshold on correlation between two variables.
#' # Also, slightly experimental as alternative to saturated model may
#' # be difficult to estimate.
#' bic.res2 <- redmod_pair(bfisub,"A2","A3",method="lv",index="bic", thresh=3,
#'                         corconst = .5, optim.method="BFGS")
#'                         
#' # Decision in either text or boolean format
#' bic.res2$Decision
#' bic.res2$Redundant      
#' 
#' # BIC for redundant and saturated (not redundant) models
#' bic.res2$Selection$idx
#' bic.res2$Selection$idxsat
#' 
#' # Difference in BIC (here, redundant model wins by over 4)
#' bic.res2$Selection$idxsat - bic.res2$Selection$idx 
#' 
#' # possible error codes in fitting and model selection
#' # bic.res2$fiterr
#' # bic.res2$selerr
#' 
#' # Where syntax is for fitting both models
#' # bic.res2$Syntax
#' 
#' # Where fitted models are
#' # bic.res2$Fit
#' 
#' # With different variables
#' bic.res3 <- redmod_pair(bfisub,"A1","A2",method="lv",index="bic", thresh=3,
#'                         optim.method="BFGS")
#' bic.res4 <- redmod_pair(bfisub,"A2","A1",method="lv",index="bic", thresh=3,
#'                         corconst = .5, optim.method="BFGS")
#' 
#' # Use model-based approach based on proportionality of covariances,
#' # sample covariances for starting values. Anecdotally, appears to have
#' # more estimation difficulty than latent variable representation.
#' #bic.res5 <- redmod_pair(bfisub,"A2","A3",method="pconst",index="bic", thresh=3,
#' #                        warmstart=TRUE, optim.method="BFGS")
#'
#' 
#' # A Bayesian example (very experimental)
#' library(blavaan)
#' # future::plan("multisession")
#' #loo.res1 <- redmod_pair(bfisub,"A1","A2",software="blavaan",
#' #               method="lv",index="loo",thresh=3)
#'
#'
#' }
#TODO: warning handling
redmod_pair <- function(dat, v1, v2, method=c("lv","pconst","sat"), software=c("lavaan","blavaan"),
                        index=c("aic","bic","sic","waic","loo"), thresh=0,
                        warmstart=FALSE, corconst=NULL, ...){

  # output container and error and warning flags
  out <- list(fiterr = FALSE,
              selerr = FALSE)

  # generate syntax
  out$Syntax <- syn_redmod(dat, v1, v2, method=method, warmstart=warmstart, corconst=corconst)

  # fit models
  out$Fit <- try(fit_redmod(out$Syntax, dat, software=software, ...))

  # do model selection w/ some error handling
  if(!inherits(out$Fit, "try-error")){
    out$Selection <- try(sel_redmod(out$Fit, index=index, thresh=thresh))

    if(!inherits(out$Selection, "try-error")){
      if(is.na(out$Selection$satpreferred)){
        out$Fiterr <- TRUE
        out$Decision <- "Inconclusive"
        out$Redundant <- NA       
      } else {
        out$Decision <- ifelse(out$Selection$satpreferred, "Not Redundant", "Redundant")
        out$Redundant <- ifelse(out$Selection$satpreferred, FALSE, TRUE)
      }
    } else {
      out$Selectionerr <- TRUE
      out$Decision <- "Inconclusive"
      out$Redundant <- NA
    }

  } else {
    out$Fiterr <- TRUE
    out$Selection <- NA
    out$Decision <- "Inconclusive"
    out$Redundant <- NA
  }

  return(out)

}
