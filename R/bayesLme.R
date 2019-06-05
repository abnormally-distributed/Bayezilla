#' Bayesian Linear Mixed Effects Models
#'
#'
#' @description This is an extremely heavily modified variant of the template.jags function in the runjags package. I have added response families,
#' enabled greater control over the priors, and also set the default priors to more sensible options. I have also added the ability for the model
#' to execute automatically, much like what one would experience with the rstanarm package. I have also added automatic generation of the full
#' log-likelihood function for use with this package's IC() function (or the loo package from the stan team) as well as automatic generation of
#' the full posterior predictive distribution. Furthermore, the model Deviance is automatically calculated without need for the dic module. This is 
#' a very useful addition because the dic module will not function when using parallel cores or clusters for sampling in JAGS 4.3.0.
#' 
#' PLEASE NOTE: This function only supports random intercept models. This does not support random slopes. To use a random slopes model
#' consider using rstanarm, brms, or the blme package. Alternatively, you can generate the random intercepts model and set autorun to FALSE,
#' which will save a text file with the model. You can customize the JAGS code in that file to include random slopes and then source 
#' the file using run.jags("yourRandomSlopesModel.txt")
#' 
#' While I have made a strong effort to choose sensible weakly informative default priors, do make sure that they are sensible for your data. 
#' These defaults generally work best if your numeric covariates are centered and scaled, which also generally improves MCMC sampling 
#' performance. This package includes a scale() function that is superior to the base R scale() function in that it allows factor 
#' columns or character columns to be in the data frame and simply does not try to scale non-numeric columns. 
#' 
#' @param formula a formula representation of the desired model, using lme4 style syntax.  Two-way interactions for all variables are permitted, as are random intercepts but NOT random slopes.
#' @param data a data frame containing the variables specified in formula.  This must be specified.
#' @param scale should the data be scaled for you? Defaults to FALSE.
#' @param file the filename of the model to output. Defaults to "Bayezilla_model.txt". This will be over-written if it already exists. 
#' @param family the response family - one of "gaussian", "student_t", "laplace", "binomial", "poisson", "negative binomial", "ZIB", "ZIP", "ZINB" (the latter denote zero-inflated distributions).
#' @param write.data option to write the data to file with the model.  If the data is very large it may be better not to write this to file, but the same data frame must be given to the subsequent run.jags call that runs the model.
#' @param write.inits option to write the initial values to file with the model.
#' @param intercept.prior the prior on the model intercept. Defaults to "dnorm(0,1)".
#' @param coef.prior the prior distribution to be used for linear and fixed effect terms, as well as interactions. Defaults to "dt(0, 1, 1)".
#' @param ranef.prior the prior distribution on random effects variance components. Defaults to "dscaled.gamma(1, 1)".
#' @param kappa.prior the prior for the negative binomial dispersion parameter. Defaults to "dmouch(1)"
#' @param tau.prior the prior distribution to be used for precision parameters. Defaults to "dgamma(.001, .001)".
#' @param nu.prior the prior distribution on the normality parameter for the student_t family. Defaults to "dgamma(2, 0.01)".
#' @param n.chains the number of chains to use.
#' @param tau.inits a numeric vector of initial values from which the precision parameters in the model will be randomly chosen.  It is recommended to make these over-dispersed, but if the values are too extreme the model may not compile.
#' @param coef.inits  a numeric vector of initial values from which the effect parameters in the model will be randomly chosen.  It is recommended to make these over-dispersed, but if the values are too extreme the model may not compile.
#' @param inits an optional list of named lists to specify initial values for one or more parameters in each chain.  The number of named lists must match n.chains.
#' @param autorun should the model automatically execute (defaults to TRUE) or would you like to first save the model for custom modification (select FALSE for this)..
#' @param iter How many post-warmup samples? Defaults to 10000. (for autorun)
#' @param warmup How many warmup samples? Defaults to 1000. (for autorun)
#' @param adapt How many adaptation steps? Defaults to 2000. (for autorun)
#' @param chains How many chains? Defaults to 4. Max allowed is 4. (for autorun)
#' @param thin Thinning interval. Defaults to 3. (for autorun)
#' @param method Defaults to "parallel". For an alternative parallel option, choose "rjparallel". Otherwise, "rjags" (single core run). (for autorun)
#' @param cl Use parallel::makeCluster(# clusters) to specify clusters for the parallel methods. Defaults to two cores. (for autorun)
#' 
#' 
#' @return
#' a runjags object
#' @export
#'
#' @examples(
#' bayesLme()
bayesLme <- function (formula, 
                      data, 
                      scale = FALSE,
                      family = "gaussian", 
                      file = "Bayezilla_model.txt", 
                      intercept.prior = "dnorm(0, 1)", 
                      tau.prior = "dgamma(.001, .001)",
                      ranef.prior = "dscaled.gamma(1, 8)", 
                      nu.prior = "dgamma(2, 0.01)",
                      kappa.prior = "dmouch(1)", 
                      coef.prior = "dt(0, 1, 1)",
                      write.data = TRUE, 
                      write.inits = TRUE,
                      tau.inits = c(0.5, 1), 
                      coef.inits = c(runif(10, -1, 1)), 
                      inits = NULL, 
                      autorun = TRUE, 
                      iter = 15000, 
                      warmup = 5000, 
                      adapt = 10000, 
                      chains = 4, 
                      thin = 10,
                      method = "rjparallel", 
                      cl = makeCluster(2)) 
{
  
  n.chains <- chains
  formula <- as.formula(formula)
  if (scale == TRUE){
    data <- scale(data)
  }
  data <- as.data.frame(data)
  if (length(as.character(formula)) != 3) 
    stop("Unsupported formula expression")
  if (length(tau.inits) < 2 || length(coef.inits) < 
      2) 
    stop("At least two precision and effect initial values must be provided")
  if (n.chains < 2) 
    stop("Two or more chains are required for standard methods of assessing convergence")
  if (is.null(inits)) 
    inits <- lapply(1:n.chains, function(x) return(list()))
  if (is.list(inits) && !any(sapply(inits, is.list))) 
    inits <- lapply(1:n.chains, function(x) return(inits))
  if (!is.list(inits) || length(inits) != n.chains || !all(sapply(inits, 
                                                                  is.list))) 
    stop("If initial values are provided to the bayesLme function this must be as a named list of length equal to the number of chains", 
         call. = FALSE)
  passedinits <- inits
  convertstrings <- FALSE
  if (class(family) != "character") 
    stop("Invalid family specification - a character string must be supplied")
  possibles <- c("gaussian", "laplace", "binomial","student_t", "poisson", "nb", "negative binomial", 
                 "zib", "zip", "zinb")
  family <- possibles[pmatch(tolower(family), possibles)]
  if (length(family) != 1 || is.na(family)) 
    stop("Invalid family specification - consult the help file for the possibilites.  Additional families (and alternative link functions) can be used by manually editing a template file created with a supported family.")
  if (family == "negative binomial") 
    family <- "nb"
  zifamily <- FALSE
  if (family == "zib") {
    zifamily <- TRUE
    family <- "binomial"
  }
  if (family == "zip") {
    zifamily <- TRUE
    family <- "poisson"
  }
  if (family == "zinb") {
    zifamily <- TRUE
    family <- "nb"
  }
  terms <- gsub("[[:space:]]", "", attr(terms(formula), "term.labels"))
  if (any(grepl("+", terms, fixed = TRUE) | grepl("*", terms, 
                                                  fixed = TRUE))) 
    stop("There was a problem parsing the formula - did you forget the parentheses around a random effects term e.g. (1 | random) ?")
  Intercept <- attr(terms(formula), "intercept")
  response <- as.character(formula)[2]
  if (grepl("+", response, fixed = TRUE) || grepl("(", response, 
                                                  fixed = TRUE) || grepl("[", response, fixed = TRUE)) 
    stop("Unsupported response expression")
  if (any(nchar(gsub("[^:]", "", terms)) > 1)) 
    stop("Unsupported 3+ way interaction term")
  if (any(grepl("I(", terms, fixed = TRUE))) 
    stop("The I( ) construct is not supported - provide variables in a pre-calculated form, or manually edit the template file")
  if (!any(response == names(data))) 
    stop("The response variable was not found in the data")
  offsets <- gsub("[[:space:]]", "", strsplit(as.character(formula)[3], 
                                              "+", fixed = TRUE)[[1]])
  offsets <- gsub("offset(", "", offsets[grepl("offset(", offsets, 
                                               fixed = TRUE)], fixed = TRUE)
  if (any(grepl("I(", offsets, fixed = TRUE))) 
    stop("The I( ) construct is not supported - provide variables in a pre-calculated form, or manually edit the template file")
  if (any(grepl("(", offsets, fixed = TRUE))) 
    stop("Functions of offset terms (and parentheses) are not supported - provide variables in a pre-calculated form, or manually edit the template file")
  offsets <- gsub(")", "", offsets, fixed = TRUE)
  notfound <- !offsets %in% names(data)
  if (any(notfound)) 
    stop(paste("The following offset term(s) was/were not found in the data: ", 
               paste(offsets[notfound], collapse = ", ")))
  missingwarn <- FALSE
  randoms <- terms[grepl("|", terms, fixed = TRUE)]
  if (any(sapply(strsplit(randoms, "|", fixed = TRUE), length) != 
          2)) 
    stop("Random slope terms are not supported by this function - but you can add these by manually editing the template file created with only random Intercept terms")
  if (any(sapply(strsplit(randoms, "|", fixed = TRUE), function(x) return(gsub("[[:space:]]", 
                                                                               "", x[1]))) != "1")) 
    stop("Random slope terms are not supported by this function - but you can add these by manually editing the template file created with only random Intercept terms")
  madefactors <- character(0)
  if (length(randoms) > 0) {
    if (grepl("+", randoms, fixed = TRUE) || grepl("/", randoms, 
                                                   fixed = TRUE) || grepl(":", randoms, fixed = TRUE) || 
        grepl("*", randoms, fixed = TRUE) || grepl("-", randoms, 
                                                   fixed = TRUE)) 
      stop("Unsupported random Intercept expression - random effects must take the form (1 | group) - consult the help file for more information")
    randoms <- sapply(strsplit(randoms, "|", fixed = TRUE), 
                      function(x) return(gsub("[[:space:]]", "", x[2])))
    notfound <- !randoms %in% names(data)
    if (any(notfound)) 
      stop(paste("The following random effects term(s) was/were not found in the data: ", 
                 paste(randoms[notfound], collapse = ", ")))
    for (i in 1:length(randoms)) {
      if (class(data[[randoms[i]]]) != "factor") {
        data[[randoms[i]]] <- factor(data[[randoms[i]]])
        madefactors <- c(madefactors, randoms[i])
      }
      if (any(is.na(data[[randoms[i]]]))) 
        missingwarn <- TRUE
    }
  }
  if (length(madefactors) > 0) 
    warning("One or more random effects terms were coerced into factors")
  termvars <- lapply(strsplit(terms[!grepl("|", terms, fixed = TRUE)], 
                              ":", fixed = TRUE), sort)
  allvars <- unique(unlist(termvars))
  if (any(grepl("I(", allvars, fixed = TRUE))) 
    stop("The I( ) construct is not supported - provide variables in a pre-calculated form")
  if (any(grepl("(", allvars, fixed = TRUE))) 
    stop("Functions of terms (and parentheses) are not supported - provide variables in a pre-calculated form")
  notfound <- !allvars %in% names(data)
  if (any(notfound)) 
    stop(paste("The following term(s) is/are not in the data: ", 
               paste(allvars[notfound], collapse = ", ")))
  missingwarn <- missingwarn || any(sapply(data[allvars], function(x) return(any(is.na(x)))))
  classes <- sapply(data[allvars], class)
  classes[classes == "integer"] <- "numeric"
  for (i in which(classes == "character")) {
    data[allvars[i]] <- factor(data[allvars[i]])
    madefactors <- c(madefactors, allvars[i])
    classes[i] <- "factor"
    stop("Some of the linear and/or fixed effects in the data frame are character variables - please convert these to factors manually, choosing the most appropriate category as the reference level")
  }
  supported <- classes %in% c("numeric", "factor")
  if (any(!supported)) 
    stop(paste("The following unsupported term classes were found in the data: ", 
               paste(classes[!supported], collapse = ", ")))
  names(classes) <- allvars
  variances <- numeric(0)
  centerwarn <- FALSE
  for (i in which(classes == "numeric")) {
    mean <- mean(data[[allvars[i]]], na.rm = TRUE)
    variance <- var(data[[allvars[i]]], na.rm = TRUE)
    variances <- c(variances, variance)
    if (any(is.na(data[[allvars[i]]]))) 
      missingwarn <- TRUE
  }
   if (length(variances) > 1 && any(variances/max(variances, 
                                                 na.rm = TRUE) < 0.01, na.rm = TRUE)) 
    warning("There is a marked discrepancy in the variance of the numeric predictor variables - it may help convergence to re-scale predictor variables")
  if (missingwarn) 
    warning("One or more of the predictor variables contains missing values - a prior distribution must be placed on these within the model before it will run. Alternatively, you can re-run the bayesLme call after running na.omit on your data to remove the missing values.")
  termtypes <- sapply(termvars, function(x) {
    stopifnot(length(x) %in% 1:2)
    if (length(x) == 1) 
      return(classes[x])
    if (length(x) == 2 && all(c("numeric", "factor") %in% 
                              classes[x])) 
      return("b_int")
    if (length(x) == 2 && all(classes[x] == "numeric")) 
      return("n_int")
    if (length(x) == 2 && all(classes[x] == "factor")) 
      return("f_int")
    return(NA)
  })
  for (i in which(termtypes == "f_int")) {
    labels <- termvars[[i]]
    matches <- termtypes == "f_int" & sapply(termvars, function(x) return(any(grepl(labels[1], 
                                                                                    x)) || any(grepl(labels[2], x))))
    matches[i] <- FALSE
    if (!any(matches)) {
      termtypes[i] <- "f_matrix"
      termtypes[sapply(termvars, length) == 1 & sapply(termvars, 
                                                       function(x) return(x[1])) == labels[1]] <- "dropped_f"
      termtypes[sapply(termvars, length) == 1 & sapply(termvars, 
                                                       function(x) return(x[1])) == labels[2]] <- "dropped_f"
    }
  }
  for (i in which(termtypes == "b_int")) {
    labels <- termvars[[i]]
    linvar <- labels[which(classes[labels] == "numeric")]
    matches <- termtypes == "b_int" & sapply(termvars, function(x) return(any(grepl(linvar, 
                                                                                    x))))
    matches[i] <- FALSE
    if (!any(matches)) {
      termtypes[i] <- "n_matrix"
      termtypes[sapply(termvars, length) == 1 & sapply(termvars, 
                                                       function(x) return(x[1])) == linvar] <- "dropped_n"
    }
  }
  varvalues <- vector("list", length = n.chains)
  varnames <- character(0)
  if (zifamily) {
    extraline <- "\n\tlinpred[i] <- regression_positive[i] * non_zero_group[i]"
    respline <- "regression_positive[i]"
  }
  else {
    extraline <- ""
    respline <- "linpred[i]"
  }
  extradata <- list(N = nrow(data))
  
  if (family == "gaussian") {
    if (!class(data[[response]]) %in% c("numeric", "integer")) 
      stop("The response variable class must be either numeric or integer for the Gaussian family")
    respline <- paste("ySim[i] ~ dnorm(linpred[i], tau)\n", 
                      "log_lik[i] <- dnorm(", response, "[i], linpred[i], tau)\n", 
                      "\t", response, "[i] ~ dnorm(linpred[i], tau)\n\tregression_residual[i] <- ", response,
                      "[i] - linpred[i]\n\tlinpred[i] <- ", 
                      sep = "")
    priorline <- paste("tau ~ ", tau.prior, 
                       "\n", sep = "")
    varnames <- c(varnames, "tau")
    signs <- sample(tau.inits, n.chains, replace = TRUE)
    for (c in 1:n.chains) varvalues[[c]] <- c(varvalues[[c]], 
                                              list(signs[c]))
  }
  if (family == "laplace") {
    if (!class(data[[response]]) %in% c("numeric", "integer")) 
      stop("The response variable class must be either numeric or integer for the Laplace family")
    respline <- paste("ySim[i] ~ ddexp(linpred[i], tau)\n", 
                      "log_lik[i] <- ddexp(", response, "[i], linpred[i], tau)\n", 
                      "\t", response, "[i] ~ ddexp(linpred[i], tau)\n\tregression_residual[i] <- ", response,
                      "[i] - linpred[i]\n\tlinpred[i] <- ", 
                      sep = "")
    priorline <- paste("tau ~ ", tau.prior, 
                       "\n", sep = "")
    varnames <- c(varnames, "tau")
    signs <- sample(tau.inits, n.chains, replace = TRUE)
    for (c in 1:n.chains) varvalues[[c]] <- c(varvalues[[c]], 
                                              list(signs[c]))
  }
  if (family == "student_t") {
    if (!class(data[[response]]) %in% c("numeric", "integer")) 
      stop("The response variable class must be either numeric or integer for the Student-t family")
    respline <- paste("ySim[i] ~ dt(linpred[i], tau, nu)\n", 
                      "log_lik[i] <- dt(", response, "[i], linpred[i], tau, nu)\n", "\t", response ,
                      "[i] ~ dt(linpred[i], tau, nu)\n\tregression_residual[i] <- ", 
                      response, "[i] - linpred[i]\n\tlinpred[i] <- "
                      , 
                      sep = "") 
    
     priorline <- paste("\ntau ~ ", tau.prior, 
                       "\n", "\n \nnu ~ ", nu.prior, "\n", sep = "")
    
    varnames <- c(varnames, c("tau", "nu"))
    varvalues <- lapply(varvalues, function(x) c(x, 3))
    signs <- c(sample(tau.inits, n.chains, replace = TRUE))
    for (c in 1:n.chains) varvalues[[c]] <- c(varvalues[[c]], 
                                              list(signs[c]))
  }
  if (family == "binomial") {
    if (length(offsets) > 0) 
      warning("Using an offset() with a logistic regression model is not recommended - did you mean to specify the response variable as a matrix with columns for successes and failures?")
    ok <- FALSE
    if (class(data[[response]]) == "matrix") {
      if (ncol(data[[response]]) != 2) 
        stop("If the response is a matrix, it must have exactly 2 columns")
      data$Binomial_Total <- data[[response]][, 1] + data[[response]][, 
                                                                      2]
      if (any(is.na(data$Binomial_Total))) 
        stop("Missing values are not allowed in the total number of trials")
      if (!all(is.numeric(data$Binomial_Total)) || any(abs(as.integer(data$Binomial_Total) - 
                                                           data$Binomial_Total) > 0.001)) 
        stop("Unexpected non integer value in the numeric response variable (or Binomial_Total variable)")
      data$Binomial_Total <- as.integer(data$Binomial_Total)
      data[[response]] <- data[[response]][, 1]
      if (!write.data) 
        stop("The data must be written to file when supplying a Binomial response as a matrix")
    }
    else {
      if (is.null(data$Binomial_Total)) 
        data$Binomial_Total <- 1
      if (any(is.na(data$Binomial_Total))) 
        stop("Missing values are not allowed in the total number of trials")
      if (!all(is.numeric(data$Binomial_Total)) || any(abs(as.integer(data$Binomial_Total) - 
                                                           data$Binomial_Total) > 0.001)) 
        stop("Unexpected non integer value in the supplied Binomial_Total variable")
      data$Binomial_Total <- as.integer(data$Binomial_Total)
    }
    if (all(data$Binomial_Total == 1) && zifamily) 
      stop("The ZIB model is only available for data with multiple trials - try specifying the data as a matrix, or using Binomial_Total to denote the total number of trials")
    if (!write.data && class(data[[response]]) %in% c("factor", 
                                                      "logical")) 
      stop("The data must be written to file when supplying a Binomial response as a factor or logical")
    if (class(data[[response]]) == "factor") {
      data[[response]] <- as.numeric(data[[response]]) - 
        1
      if (any(data[[response]] > 1, na.rm = TRUE)) 
        warning("Grouping factor levels 2 and above in the response variable")
      ok <- TRUE
    }
    if (class(data[[response]]) %in% c("logical", "numeric", 
                                       "integer")) {
      if (any(abs(as.integer(data[[response]]) - data[[response]]) > 
              0.001, na.rm = TRUE)) 
        stop("Warning: Unexpected non integer value in the numeric response variable (or Binomial_Total variable)")
      data[[response]] <- as.integer(data[[response]])
      if (all(data$Binomial_Total == 1) && !all(data[[response]] %in% 
                                                c(0, 1), na.rm = TRUE)) 
        stop("Warning: non 0-1 value in the numeric response variable")
      if (any(data$Binomial_Total < data[[response]], na.rm = TRUE)) 
        stop("Warning: non 0-1 value in the numeric response variable")
      ok <- TRUE
    }
    priorline <- ""
    if (!ok) 
      stop("Unrecognised response variable format - possibilities are a matrix with columns for successes and failures, or a factor or numeric variable assuming one trial per observation")
    respline <- paste("\t", response, "[i] ~ dbin(regression_prob[i], Binomial_Total[i])\n\tregression_residual[i] <- ", 
                      response, "[i] - linpred[i]\n\tlinpred[i] <- regression_prob[i] * Binomial_Total[i]\n", 
                      if (zifamily) 
                        "\tregression_prob[i] <- non_zero_regression[i] * non_zero_group[i]\n\tlogit(non_zero_regression[i]) <- "
                      else "\tlogit(regression_prob[i]) <- ", sep = "")
    extradata <- c(extradata, list(Binomial_Total = data$Binomial_Total))
  }
  if (family == "poisson") {
    if (!class(data[[response]]) %in% c("numeric", "integer")) 
      stop("The response variable class must be either numeric or integer for the poisson family")
    if (any(data[[response]] < 0, na.rm = TRUE) || any(as.integer(data[[response]]) != 
                                                       data[[response]], na.rm = TRUE)) 
      stop("Only positive integers are allowed in the response variable for the poisson family")
    respline <- paste("\t", response, "[i] ~ dpois(linpred[i])\n\tregression_residual[i] <- ", 
                      response, "[i] - linpred[i]\n\t", if (zifamily) 
                        "linpred[i] <- non_zero_regression[i] * non_zero_group[i]\n\tlog(non_zero_regression[i]) <- "
                      else "log(linpred[i]) <- ", sep = "")
    priorline <- ""
  }
  if (family == "nb") {
    if (!class(data[[response]]) %in% c("numeric", "integer")) 
      stop("The response variable class must be either numeric or integer for the Negative Binomial family")
    if (any(data[[response]] < 0, na.rm = TRUE) || any(as.integer(data[[response]]) != 
                                                       data[[response]], na.rm = TRUE)) 
      stop("Only positive integers are allowed in the response variable for the poisson family")
    respline <- paste("\t", response, "[i] ~ dpois(linpred[i])\n\tregression_residual[i] <- ", 
                      response, "[i] - linpred[i]\n\tdispersion[i] ~ dgamma(k, k)\n\tlinpred[i] <- regression_mean[i] * dispersion[i]", 
                      if (zifamily) 
                        " * non_zero_group[i]", "\n\t# Note: this formulation of a gamma-poisson is exactly equivalent to a Negative Binomial\n\tlog(regression_mean[i]) <- ", 
                      sep = "")
    priorline <- paste("kappa ~", kappa.prior, sep = "")
    varnames <- c(varnames, "kappa")
    signs <- sample(tau.inits, n.chains, replace = TRUE)
    for (c in 1:n.chains) varvalues[[c]] <- c(varvalues[[c]], 
                                              list(signs[c]))
  }
  if (Intercept != 0) {
    respline <- paste(respline, "Intercept + ", sep = "")
    priorline <- paste(priorline, "Intercept ~ ", intercept.prior, 
                       "\n", sep = "")
    varnames <- c(varnames, "Intercept")
    signs <- sample(coef.inits, n.chains, replace = TRUE)
    for (c in 1:n.chains) varvalues[[c]] <- c(varvalues[[c]], 
                                              list(signs[c]))
  }
  for (i in which(termtypes == "numeric")) {
    
    respline <- paste(respline, termvars[[i]][1], '_coefficient * ', termvars[[i]][1], '[i] + ', sep='')
    priorline <- paste(priorline, termvars[[i]][1], '_coefficient ~ ', coef.prior, '\n', sep='')
    
    varnames <- c(varnames, paste(termvars[[i]][1], '_coefficient', sep=''))
    signs <- sample(coef.inits, n.chains, replace=TRUE)
  
      for(c in 1:n.chains)
      varvalues[[c]] <- c(varvalues[[c]], list(signs[c]))		
  
    for (c in 1:n.chains) varvalues[[c]] <- c(varvalues[[c]], 
                                              list(signs[c]))
  }
  for (i in which(termtypes == "n_int")) {
    respline <- paste(respline, termvars[[i]][1], "_", termvars[[i]][2], 
                      "_interaction * ", termvars[[i]][1], "[i] * ", termvars[[i]][2], 
                      "[i] + ", sep = "")
    priorline <- paste(priorline, termvars[[i]][1], "_", 
                       termvars[[i]][2], "_interaction ~ ", coef.prior, 
                       "\n", sep = "")
    varnames <- c(varnames, paste(termvars[[i]][1], "_", 
                                  termvars[[i]][2], "_interaction", sep = ""))
    signs <- sample(coef.inits, n.chains, replace = TRUE)
    for (c in 1:n.chains) varvalues[[c]] <- c(varvalues[[c]], 
                                              list(signs[c]))
  }
  for (i in which(termtypes == "n_matrix")) {
    linvar <- termvars[[i]][which(classes[termvars[[i]]] == 
                                    "numeric")]
    fixvar <- termvars[[i]][which(classes[termvars[[i]]] == 
                                    "factor")]
    respline <- paste(respline, linvar, "_coefficient_", 
                      fixvar, "_level[", fixvar, "[i]] * ", linvar, "[i] + ", 
                      sep = "")
    factpriors <- rep(paste(" ~ ", coef.prior, sep = ""), 
                      length(levels(data[[fixvar]])))
    factpriors <- paste(linvar, "_coefficient_", fixvar, 
                        "_level[", 1:length(factpriors), "]", factpriors, 
                        "    # Factor level \"", levels(data[[fixvar]]), 
                        "\"", sep = "")
    priorline <- paste(priorline, paste(factpriors, collapse = "\n"), 
                       "\n", sep = "")
    varnames <- c(varnames, paste(linvar, "_coefficient_", 
                                  fixvar, "_level", sep = ""))
    for (c in 1:n.chains) {
      signs <- sample(coef.inits, length(levels(data[[fixvar]])), 
                      replace = TRUE)
      varvalues[[c]] <- c(varvalues[[c]], list(signs))
    }
  }
  for (i in which(termtypes == "factor")) {
    respline <- paste(respline, termvars[[i]][1], "_effect[", 
                      termvars[[i]][1], "[i]] + ", sep = "")
    factpriors <- rep(paste(" ~ ", coef.prior, sep = ""), 
                      length(levels(data[[termvars[[i]][1]]])))
    #factpriors[1] <- " <- 0"
    factpriors <- paste(termvars[[i]][1], "_effect[", 1:length(factpriors), 
                        "]", factpriors, "    # Factor level \"", levels(data[[termvars[[i]][1]]]), 
                        "\"", sep = "")
    priorline <- paste(priorline, paste(factpriors, collapse = "\n"), 
                       "\n", sep = "")
    varnames <- c(varnames, paste(termvars[[i]][1], "_effect", 
                                  sep = ""))
    for (c in 1:n.chains) {
      signs <- sample(coef.inits, length(levels(data[[termvars[[i]][1]]])), 
                      replace = TRUE)
      signs[1] <- NA
      varvalues[[c]] <- c(varvalues[[c]], list(signs))
    }
  }
  for (i in which(termtypes == "f_matrix")) {
    respline <- paste(respline, termvars[[i]][1], "_", termvars[[i]][2], 
                      "_effect[", termvars[[i]][1], "[i],", termvars[[i]][2], 
                      "[i]] + ", sep = "")
    factindices <- expand.grid(1:length(levels(data[[termvars[[i]][1]]])), 
                               1:length(levels(data[[termvars[[i]][2]]])))
    factpriors <- rep(paste(" ~ ", coef.prior, sep = ""), 
                      nrow(factindices))
    factpriors[1] <- " <- 0"
    factpriors <- paste(termvars[[i]][1], "_", termvars[[i]][2], 
                        "_effect[", factindices[, 1], ",", factindices[, 
                                                                       2], "] ", factpriors, "    # Factor level \"", 
                        levels(data[[termvars[[i]][1]]])[factindices[, 1]], 
                        "\", \"", levels(data[[termvars[[i]][2]]])[factindices[, 
                                                                               2]], "\"", sep = "")
    priorline <- paste(priorline, paste(factpriors, collapse = "\n"), 
                       "\n", sep = "")
    varnames <- c(varnames, paste(termvars[[i]][1], "_", 
                                  termvars[[i]][2], "_effect", sep = ""))
    for (c in 1:n.chains) {
      signs <- matrix(sample(coef.inits, nrow(factindices), 
                             replace = TRUE), nrow = length(levels(data[[termvars[[i]][1]]])))
      signs[1, ] <- NA
      signs[, 1] <- NA
      varvalues[[c]] <- c(varvalues[[c]], list(signs))
    }
  }
  for (i in which(termtypes == "f_int")) {
    respline <- paste(respline, termvars[[i]][1], "_", termvars[[i]][2], 
                      "_interaction[", termvars[[i]][1], "[i],", termvars[[i]][2], 
                      "[i]] + ", sep = "")
    factindices <- expand.grid(1:length(levels(data[[termvars[[i]][1]]])), 
                               1:length(levels(data[[termvars[[i]][2]]])))
    factpriors <- rep(paste(" ~ ", coef.prior, sep = ""), 
                      nrow(factindices))
    factpriors[factindices[, 1] == 1] <- " <- 0"
    factpriors[factindices[, 2] == 1] <- " <- 0"
    factpriors <- paste(termvars[[i]][1], "_", termvars[[i]][2], 
                        "_interaction[", factindices[, 1], ",", factindices[, 
                                                                            2], "] ", factpriors, "    # Factor level \"", 
                        levels(data[[termvars[[i]][1]]])[factindices[, 1]], 
                        "\", \"", levels(data[[termvars[[i]][2]]])[factindices[, 
                                                                               2]], "\"", sep = "")
    priorline <- paste(priorline, paste(factpriors, collapse = "\n"), 
                       "\n", sep = "")
    varnames <- c(varnames, paste(termvars[[i]][1], "_", 
                                  termvars[[i]][2], "_effect", sep = ""))
    signs <- sample(tau.inits, n.chains, replace = TRUE)
    for (c in 1:n.chains) {
      signs <- matrix(sample(coef.inits, nrow(factindices), 
                             replace = TRUE), nrow = length(levels(data[[termvars[[i]][1]]])))
      signs[1, ] <- NA
      signs[, 1] <- NA
      varvalues[[c]] <- c(varvalues[[c]], list(signs))
    }
  }
  
  
  # Then random effects:
  for(r in randoms){
    respline <- paste(respline, r, '_randomeffect[', r, '[i]] + ', sep='')
    priorline <- paste(priorline, 'for(', r, '_iterator in 1:', length(levels(data[[r]])), '){\n\t', r,
                       '_randomeffect[', r, '_iterator] ~ dnorm(0, ', r, '_precision)\n}\n', r, 
                       '_precision ~ ', ranef.prior, '\n', sep='') 
    
    varnames <- c(varnames, paste(r, '_precision', sep=''))
    signs <- sample(tau.inits, n.chains, replace=TRUE)
    for(c in 1:n.chains)
      varvalues[[c]] <- c(varvalues[[c]], list(signs[c]))		
  }
  
  # And offsets:
  for(o in offsets){
    respline <- paste(respline, o, '[i] + ', sep='')
  }
  
  # Horrible hack to get rid of the trailing +:
  respline <- paste(respline, '_ + _ +\n', sep='')
  respline <- gsub('+ _ + _ +', '', respline, fixed=TRUE)
  
  if(zifamily){
    respline <- paste(respline, '\tnon_zero_group[i] ~ dbern(non_zero_prob[i])\n\tlogit(non_zero_prob[i]) <- -(zero_inflation_intercept)\n\t\t# Note: this line (inside the parentheses) could specify a separate linear regression\n\t\t# To make this the probability of zero-inflation, the - symbol is required!\n', sep='')
    priorline <- paste(priorline, 'zero_inflation_intercept ~ ', coef.prior, '\nnon_zero_propotion <- ilogit(-zero_inflation_intercept)\n', sep='')
    signs <- sample(coef.inits, n.chains, replace=TRUE)
    zistarts <- rep(1, nrow(data))
    varnames <- c(varnames, 'zero_inflation_intercept', 'non_zero_group')
    for(c in 1:n.chains)
      varvalues[[c]] <- c(varvalues[[c]], list(signs[c], zistarts))		
  }
  
  for(c in 1:n.chains)
    names(varvalues[[c]]) <- varnames
  
  
  if (write.data) {
    magicline <- ""
    data <- dump.format(c(data[unique(c(response, allvars, 
                                        offsets, randoms))], extradata))
  }
  else {
    magicline <- paste("#data# ", paste(unique(c(response, 
                                                 allvars, offsets, randoms)), collapse = ", "), "\n", 
                       sep = "")
    data <- dump.format(extradata)
  }
  if (zifamily) {
    modules <- c("glm")
  }
  else {
    modules <- c("glm")
  }
  factories <- ""
  monitor <- c(varnames[!varnames %in% c("non_zero_group", 
                                         "zero_inflation_Intercept")], if (zifamily) c("non_zero_propotion", "Deviance", "ySim", "log_lik") else c("sigma", "Deviance", "ySim", "log_lik"))
  
  monitor = monitor[-which(monitor == "tau")]
  if (family == "student_t"){
    monitor = monitor[-which(monitor == "nu")]
    monitor = c(monitor, "nu")
  }
  
  
  for (c in 1:n.chains) {
    for (n in names(passedinits[[c]])) {
      varvalues[[c]][[n]] <- passedinits[[c]][[n]]
    }
  }
  if (write.inits) {
    end.state <- sapply(varvalues, dump.format)
  }
  else {
    end.state <- ""
  }
  model <- paste("\n\nmodel{\n\n\nfor(i in 1:N){\n\n", 
                 respline, "}\n\n# These lines give the prior distributions for the parameters to be estimated:\n", 
                 priorline, "\nsigma <- sqrt(1/tau)\n", "\nDeviance <- -2 * sum(log_lik[1:N])\n",
                 magicline, "\n}\n\n# These lines are hooks to be read by runjags (they are ignored by JAGS):", 
                 sep = "")
  rjo <- list(model = model, data = data, end.state = end.state, 
              monitor = monitor, modules = modules, factories = factories, 
              response = response, fitted = "linpred", residual='regression_residual')
  class(rjo) <- "runjags"
  write.jagsfile(runjags.object = rjo, file = file, remove.tags = FALSE)
  
  if (autorun == FALSE){
    cat("Your model template was created at", file)
    if (write.data) 
      cat("You can customize it further and then run the model using run.jags(", file, ")\n", sep = "")
    else cat("You can then run the model using run.jags(\"",  
             file, "\", data=data) - where \"data\" is the same data frame specified to the bayesLme function\n", 
             sep = "")
    invisible(file)
  }
  else if (autorun == TRUE){
    run.jags(file , burnin = warmup, sample = iter, modules = "glm", n.chains = n.chains, thin = thin, adapt = adapt, method = method, cl = cl, summarise = FALSE)
  }
  
}