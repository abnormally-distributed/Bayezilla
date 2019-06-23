#' Get WAIC and/or LOO-IC
#'
#' @description A common method of model comparison in Bayesian methods is to use the ratio
#' of the marginal likelihoods, the Bayes Factor. 
#' likelihood function itself, p(Y | theta ; M), is conditioned on the parameters, telling you the likelihood
#' of the present sample based on the present parameter estimates. The marginal likelihood, however, is
#' dependent on the structure of the model itself and not the parameter estimates. Because the 
#' parameters are integrated out of the model, it is an
#' ideal way in theory to perform model comparison without the risk of overfitting.
#'
#' Unfortunately, practical Bayesian analyses often without a way to get the marginal likelihoods.
#' Marginal likelihoods are difficult
#' to calculate, and sometimes analytically intractable.
#' An alternative to marginal likelihood based model comparison are information criteria,
#' which minimize the Kullback-Leibler divergence between the estimated distribution of the observed data and
#' the unknown population distribution. Classical information criteria such as the AIC and SIC (aka the BIC or SBIC)
#' penalize the deviance of the model evaluated at the maximum likelihood estimate, i.e.
#' -2*log( p(Y | theta_MLE ; M) ) + penalty.
#' \cr
#' Full Bayesian analysis with MCMC affords an opportunity to use extensions of the AIC that are not reliant on
#' maximum liklihood point estimates. These utilize the log pointwise predictive density, or lppd. This is given by taking the
#' product of log-likelihood function at each data point, p(Y_i | theta ; M)
#' and the log-posterior function of the parameters p(theta | Y ; M) and integrating this product
#' over the parameters. Finally, this is summed across all N_i data points.
#' \cr
#' \if{html}{\figure{lppd.png}{the log pointwise predictive density.}}
#' \if{latex}{\figure{lppd.png}{the log pointwise predictive density.}}
#' \cr
#' The Watanabe-Akaike Information Criterion substitutes the lppd for the deviance, and use a penalty 
#' based on the number of effective parameters. The number of effective
#' parameters is given by calculating the variance of each data point's log-likelihood function
#' across the MCMC samples, and then summing these variances. This is very useful since Bayesian models
#' often utilize a number of hyperparameters intended to penalize and shrink the estimates, so using the
#' number of effective parameters avoids the issue of overpenalzing such models.
#' \cr
#' \if{html}{\figure{pWAIC.png}{the WAIC penalty term.}}
#' \if{latex}{\figure{pWAIC.png}{the WAIC penalty term.}}
#' \cr
#'
#' Finally, the pWAIC term is subtracted from the lppd to obtain the WAIC, which is an estimator
#' of the expected log predictive density on the deviance scale (- 2  * elpd).
#'
#' The other alternative, which is often preferable, uses leave-one-out (LOO) cross validation to
#' obtain the expected log predictive density function, which serves as an estimator of the the log density
#' function of the population distribution from which the observed data came. In other words, this estimates the
#' predictive performance of the model on unseen data that might be obtained by sampling from the same population
#' from which the current sample came. Note that asymptotically the WAIC is equivalent to the LOOIC (Watanabe, 2010).
#'
#' The LOO-IC is found by calculating loo weights and re-weighing the lppd (which is calculated the same way as in WAIC)
#' and multiplying the reweighted lppd by -2 to put it on the deviance scale. Fortunately this can be done in
#' an extremely efficient manner without having to actually refit the model N times.
#' This reweighted lppd is the LOO-IC, an estimator of the expected log predictive density.
#' Note that the penalty term is not calculated directly as in the case of pWAIC. However, the
#' pWAIC can be extracted by taking the difference of the lppd and LOO-IC.
#' The steps of calculating the LOO-IC are shown below.
#' \cr
#' \if{html}{\figure{looic.png}{the steps in obtaining the LOO-IC.}}
#' \if{latex}{\figure{looic.png}{the steps in obtaining the LOO-IC.}}
#' \cr
#' The difference between the marginal likelihood and the expected log predictive density is that
#' the marginal likelihood integrates over the prior distribution of the parameters in the model, p(theta, M),
#' while the expected log predictive density integrates over the posterior distribution of the parameters given
#' the data and model, p(theta | Y, M). This reflects different inferential goals. The goal of the lppd
#' and elpd quantities is assessing predictive accuracy. Gelman, Wang, and Vehtari (2016) state that, "..the prior is
#' relevant in estimating the parameters but not in assessing a model’s accuracy. We are not saying that the prior
#' cannot be used in assessing a model’s fit to data; rather we say that the prior density is
#' not relevant in computing predictive accuracy." Hence, marginal likelihood based model comparison can be best
#' understood as being most appropriate when the true model, or a very good approximation to it, is among the candidate
#' set of models (M-closed and M-complete scenarios) while expected log predictive density based quantities are suitable
#' for M-open scenarios where the best one can hope for is good predictive accuracy. Regardless of inferential preference,
#' a pragmatic reason to use the WAIC or LOO-IC is that it can be calculated for any Bayesian model fit with ease.
#'
#' Below is a figure displaying the log-marginal likelihood at the top, along with the WAIC and LOO relevant quantities
#' re-expressed with full conditional probability notation for comparison.
#' \cr
#' \if{html}{\figure{likelihoods.png}{A lineup of the different likelihood based quantities.}}
#' \if{latex}{\figure{likelihoods.png}{A lineup of the different likelihood based quantities.}}
#' \cr
#'
#' @references 
#' Gelman, A., Hwang, J., and Vehtari, A. (2013). Understanding predictive information criteria for Bayesian models. Statistics and Computing. Volume 24, Issue 6, pp 997–1016 \cr
#' \cr
#' Vehtari, A., Gelman, A. (2014). WAIC and cross-validation in Stan. Online manuscript. http://www.stat.columbia.edu/~gelman/research/unpublished/waic_stan \cr
#' \cr
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable
#' information criterion in singular learning theory. Journal of Machine Learning Research 11,
#' 3571–3594. \cr
#' \cr
#' 
#'
#' @param out the stanfit or runjags objects. Must contain a log likelihood for each observation, not the
#' total log likelihood of the model.
#' @param loo Should the LOO also be returned? Defaults to TRUE.
#' @param summarize Should only the expected values be returned? Defaults to TRUE. If set to FALSE, a data frame
#' containing the pointwise WAIC and LOO-IC will be returned.
#'
#' @return
#' either a point summary or data frame
#'
#' @export
#'
#' @examples
#' IC()
#' 
IC <- function (out, loo = TRUE, summarize = TRUE) 
  {
    stan <- inherits(out, "stanfit")
    if (stan == TRUE) {
      LL <- as.matrix(out)
      LL <- LL[, -which(colnames(LL) == "lp__")]
      wch = which(regexpr("log_lik", colnames(LL)) == 1)
      if (length(wch) != 0){
        LL <- LL[,wch]
      } else if (length(wch) == 0){
        stop("Variables matching log_lik not found. Please enable log_lik = TRUE in the model function.")
      }
    }
    else if (class(out) == "runjags") {
      LL <- combine.mcmc(out, collapse.chains = TRUE)
      wch = which(regexpr("log_lik", colnames(LL)) == 1)
      if (length(wch) != 0){
        LL <- LL[,wch]
      } else if (length(wch) == 0){
        stop("Variables matching log_lik not found. Please enable log_lik = TRUE in the model function.")
      }
    }
    
    
    S <- nrow(LL)
    n <- ncol(LL)
    
    # A trick I found online to prevent numerical overflow.
    offset <- LL[cbind(max.col(abs(t(LL))), 1:n)] 
    LPPD  <- log(1./S)+log(colSums(exp(sweep(LL, 2, offset))))+offset
    
    pWAIC <- apply(LL, 2, var)
    WAIC <- -2 * (LPPD - pWAIC)
    IC <- data.frame(WAIC = WAIC, lppd = LPPD, pWAIC = pWAIC)
    
    if (loo == TRUE) {
      loo_weights_raw <- 1/exp(LL - max(LL))
      loo_weights_normalized <- loo_weights_raw/matrix(colMeans(loo_weights_raw), 
                                                       nrow = S, ncol = n, byrow = TRUE)
      loo_weights_regularized <- pmin(loo_weights_normalized, 
                                      sqrt(S))
      elpd_loo <- log(colMeans(exp(LL) * loo_weights_regularized)/colMeans(loo_weights_regularized))
      LOOIC <- -2 * elpd_loo
      p_loo <- LPPD - elpd_loo
      IC <- cbind.data.frame(WAIC = WAIC, `LOO-IC` = LOOIC, 
                             lppd = LPPD, pWAIC = pWAIC, pLOO = p_loo)
    }
    if (summarize == TRUE) {
      return(colSums(IC))
    }
    else {
      return(IC)
    }
  }
