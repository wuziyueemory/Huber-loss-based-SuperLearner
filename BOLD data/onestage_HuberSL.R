# Load necessary packages
library(SuperLearner)
library(CVXR)
library(e1071)

## Nested-CV vs. partial-CV

# one stage super-learner with Square Loss and Huber Loss 
############################################### Square Loss #########################################################
# function for generating weights (coefficients) under Square Loss: scaled quadratic programming
method.CC_LS.scale <- function() {
  computeCoef = function(Z, Y, libraryNames, verbose,
                         obsWeights=rep(1, length(Y)),
                         errorsInLibrary = NULL, ...) {
    # compute cvRisk
    cvRisk <- apply(Z, 2, function(x) mean(obsWeights*(x-Y)^2))
    names(cvRisk) <- libraryNames
    # compute coef
    compute <- function(x, y, wt=rep(1, length(y))) {
      wX <- sqrt(wt) * x
      wY <- sqrt(wt) * y
      D <- crossprod(wX)
      d <- crossprod(wX, wY)
      A <- cbind(rep(1, ncol(wX)), diag(ncol(wX)))
      bvec <- c(1, rep(0, ncol(wX)))
      sc <- norm(D,"2")
      # scale D matrix & d vector to aviod inconsistent constraints
      fit <- quadprog::solve.QP(Dmat=D/sc, dvec=d/sc, Amat=A, bvec=bvec, meq=1,
                                factorized = F)
      invisible(fit)
    }
    modZ <- Z
    # check for columns of all zeros. assume these correspond
    # to errors that SuperLearner sets equal to 0. not a robust
    # solution, since in theory an algorithm could predict 0 for
    # all observations (e.g., SL.mean when all Y in training = 0)
    naCols <- which(apply(Z, 2, function(z){ all(z == 0 ) }))
    anyNACols <- length(naCols) > 0
    if(anyNACols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[naCols],collapse = ", "), " have NAs.",
                     "Removing from super learner."))
    }
    # check for duplicated columns
    # set a tolerance level to avoid numerical instability
    tol <- 8
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    if(anyDupCols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[dupCols],collapse = ", "),
                     " are duplicates of previous learners.",
                     " Removing from super learner."))
    }
    # remove from Z if present
    if(anyDupCols | anyNACols){
      rmCols <- unique(c(naCols,dupCols))
      modZ <- Z[,-rmCols]
    }
    # compute coefficients on remaining columns
    fit <- compute(x = modZ, y = Y, wt = obsWeights)
    coef <- fit$solution
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] = 0
    }
    # add in coefficients with 0 weights for algorithms with NAs
    if(anyDupCols | anyNACols){
      ind <- c(seq_along(coef), rmCols - 0.5)
      coef <- c(coef, rep(0, length(rmCols)))
      coef <- coef[order(ind)]
    }
    # Set very small coefficients to 0 and renormalize.
    coef[coef < 1.0e-4] <- 0
    coef <- coef / sum(coef)
    if(!sum(coef) > 0) warning("All algorithms have zero weight", call. = FALSE)
    list(cvRisk = cvRisk, coef = coef, optimizer = fit)
  }
  
  computePred = function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out <- list(require = "quadprog",
              computeCoef = computeCoef,
              computePred = computePred)
  invisible(out)
}


############################################### Huber Loss #########################################################
# function for generating weights (coefficients) under Huber Loss: scaled quadratic programming
method.CC_HUBER <- function() {
  computeCoef = function(Z, Y, libraryNames, verbose,lambda,
                         obsWeights=rep(1, length(Y)),
                         errorsInLibrary = NULL, ...) {
    # compute cvRisk
    cvRisk <- apply(Z, 2, function(x) 
      mean(ifelse((abs(Y-x) > 10000*lambda),
                  10000*lambda*(obsWeights*abs(Y-x) - 0.5*10000*lambda),
                  0.5*(obsWeights*(x-Y)^2))))
    names(cvRisk) <- libraryNames
    
    modZ <- Z
    # check for columns of all zeros. assume these correspond
    # to errors that SuperLearner sets equal to 0. not a robust
    # solution, since in theory an algorithm could predict 0 for
    # all observations (e.g., SL.mean when all Y in training = 0)
    naCols <- which(apply(Z, 2, function(z){ all(z == 0 ) }))
    anyNACols <- length(naCols) > 0
    if(anyNACols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[naCols],collapse = ", "), " have NAs.",
                     "Removing from super learner."))
    }
    # check for duplicated columns
    # set a tolerance level to avoid numerical instability
    tol <- 8
    dupCols <- which(duplicated(round(Z, tol), MARGIN = 2))
    anyDupCols <- length(dupCols) > 0
    if(anyDupCols){
      # if present, throw warning identifying learners
      warning(paste0(paste0(libraryNames[dupCols],collapse = ", "),
                     " are duplicates of previous learners.",
                     " Removing from super learner."))
    }
    # remove from Z if present
    if(anyDupCols | anyNACols){
      rmCols <- unique(c(naCols,dupCols))
      modZ <- Z[,-rmCols]
    }
    
    # compute coefficients on remaining columns -- Use CVXR 
    # Variables minimized over
    beta <- Variable(ncol(modZ))
    # Problem definition
    objective <- Minimize(sum(huber((Y - modZ%*%beta)/10000, lambda)))
    constraints <- list(beta >= 0, sum(beta)==1)
    prob <- Problem(objective,constraints)
    # Problem solution
    result <- solve(prob,solver="ECOS", MAXIT=as.integer(2000))
    coef <- as.vector(result$getValue(beta))
    
    if (anyNA(coef)) {
      warning("Some algorithms have weights of NA, setting to 0.")
      coef[is.na(coef)] = 0
    }
    # add in coefficients with 0 weights for algorithms with NAs
    if(anyDupCols | anyNACols){
      ind <- c(seq_along(coef), rmCols - 0.5)
      coef <- c(coef, rep(0, length(rmCols)))
      coef <- coef[order(ind)]
    }
    # Set very small coefficients to 0 and renormalize.
    coef[coef < 1.0e-4] <- 0
    coef <- coef / sum(coef)
    if(!sum(coef) > 0) warning("All algorithms have zero weight", call. = FALSE)
    list(cvRisk = cvRisk, coef = coef, optimizer = result)
  }
  
  computePred = function(predY, coef, ...) {
    predY %*% matrix(coef)
  }
  out <- list(require = "CVXR",
              computeCoef = computeCoef,
              computePred = computePred)
  invisible(out)
}



############################## Change the Print function in SuperLearner ####################################
print.SuperLearner <- function(x, ...) {
  cat("\nCall: ", deparse(x$call, width.cutoff = .9*getOption("width")), "\n\n", fill = getOption("width"))
  print(list('Standard (MSE): one-stage super learner' = 
               cbind(Risk = x$cvRisk$`one-stage super learner`$`Standard (MSE)`,
                     Coef = x$coef$`one-stage`$`Standard (MSE)`),
             'Huber-partial CV: one-stage super learner' = 
               cbind(Risk = x$cvRisk$`one-stage super learner`$`Huber-partial CV`,
                     Coef = x$coef$`one-stage`$`Huber-partial CV`),
             'Huber-nested CV: one-stage super learner' = 
               cbind(Risk = x$cvRisk$`one-stage super learner`$`Huber-nested CV`,
                     Coef = x$coef$`one-stage`$`Huber-nested CV`))
  )
}







## One-Stage Huber loss-based Super Learner
##############################################################################################################

onestage.HuberSL <- function(Y, X, newX = NULL, SL.library, lambda,
                             family = gaussian(),
                             id=NULL, verbose=FALSE, control = list(), 
                             cvControl = list(), obsWeights = NULL, env = parent.frame()){
  
  # Begin timing how long two-stage SuperLearner takes to execute
  time_start = proc.time()
  
  # Get details of estimation algorithm for the algorithm weights (coefficients)
  method <- list('Square Loss' = method.CC_LS.scale(),
                 'Huber Loss' = method.CC_HUBER())
  
  # get defaults for controls and make sure in correct format
  control <- do.call('SuperLearner.control', control)
  cvControl <- do.call('SuperLearner.CV.control', cvControl)
  V <- cvControl$V
    
  # put together the library
  # should this be in a new environment?
  library <- SuperLearner:::.createLibrary(SL.library)
  SuperLearner:::.check.SL.library(library = c(unique(library$library$predAlgorithm), library$screenAlgorithm))
  
  call <- match.call(expand.dots = TRUE)
  # should we be checking X and newX for data.frame?
  # data.frame not required, but most of the built-in wrappers assume a data.frame
  if(!inherits(X, 'data.frame')) message('X is not a data frame. Check the algorithms in SL.library to make sure they are compatible with non data.frame inputs')
  varNames <- colnames(X)
  N <- dim(X)[1L]
  p <- dim(X)[2L]
  k <- nrow(library$library)
  kScreen <- length(library$screenAlgorithm)
  n.lambda <- length(lambda)
  Z <- matrix(NA, N, k)
  libraryNames <- paste(library$library$predAlgorithm, library$screenAlgorithm[library$library$rowScreen], sep="_")
  
  if(p < 2 & !identical(library$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }
  
  # put fitLibrary in it's own environment to locate later
  fitLibEnv <- new.env()
  assign('fitLibrary', vector('list', length = k), envir = fitLibEnv)
  assign('libraryNames', libraryNames, envir = fitLibEnv)
  evalq(names(fitLibrary) <- libraryNames, envir = fitLibEnv)
  
  # errors* records if an algorithm stops either in the CV step and/or in full data
  errorsInCVLibrary <- rep(0, k)
  errorsInLibrary <- rep(0, k)
  
  # if newX is missing, use X
  if(is.null(newX)) {
    newX <- X
  }
  # Are these checks still required?
  if(!identical(colnames(X), colnames(newX))) {
    stop("The variable names and order in newX must be identical to the variable names and order in X")
  }
  if (sum(is.na(X)) > 0 | sum(is.na(newX)) > 0 | sum(is.na(Y)) > 0) {
    stop("missing data is currently not supported. Check Y, X, and newX for missing values")
  }
  if (!is.numeric(Y)) {
    stop("the outcome Y must be a numeric vector")
  }
  # family can be either character or function, so these lines put everything together (code from glm())
  if(is.character(family))
    family <- get(family, mode="function", envir=parent.frame())
  if(is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if (family$family != "binomial" & isTRUE("cvAUC" %in% method$require)){
    stop("'method.AUC' is designed for the 'binomial' family only")
  }
  
  # create CV folds
  validRows <- CVFolds(N = N, id = id, Y = Y, cvControl = cvControl)
  
  # test id
  if(is.null(id)) {
    id <- seq(N)
  }
  if(!identical(length(id), N)) {
    stop("id vector must have the same dimension as Y")
  }
  # test observation weights
  if(is.null(obsWeights)) {
    obsWeights <- rep(1, N)
  }
  if(!identical(length(obsWeights), N)) {
    stop("obsWeights vector must have the same dimension as Y")
  }
  
  ########################################################################################################
  # step 1: Model fit via cross-validation
  # create function for the cross-validation step:
  .crossValFUN <- function(valid, Y, dataX, id, obsWeights, library, 
                           kScreen, k, p, libraryNames, saveCVFitLibrary) {
    tempLearn <- dataX[-valid, , drop = FALSE]
    tempOutcome <- Y[-valid]
    tempValid <- dataX[valid, , drop = FALSE]
    tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
    tempId <- id[-valid]
    tempObsWeights <- obsWeights[-valid]
    
    # should this be converted to a lapply also?
    for(s in seq(kScreen)) {
      screen_fn = get(library$screenAlgorithm[s], envir = env)
      testScreen <- try(do.call(screen_fn, list(Y = tempOutcome, X = tempLearn, family = family, id = tempId, obsWeights = tempObsWeights)))
      if(inherits(testScreen, "try-error")) {
        warning(paste("replacing failed screening algorithm,", library$screenAlgorithm[s], ", with All()", "\n "))
        tempWhichScreen[s, ] <- TRUE
      } else {
        tempWhichScreen[s, ] <- testScreen
      }
      if(verbose) {
        message(paste("Number of covariates in ", library$screenAlgorithm[s], " is: ", sum(tempWhichScreen[s, ]), sep = ""))
      }
    } #end screen
    
    # should this be converted to a lapply also?
    out <- matrix(NA, nrow = nrow(tempValid), ncol = k)
    if(saveCVFitLibrary){
      model_out <- vector(mode = "list", length = k)
    }else{
      model_out <- NULL
    }
    
    for(s in seq(k)) {
      pred_fn = get(library$library$predAlgorithm[s], envir = env)
      testAlg <- try(do.call(pred_fn, list(Y = tempOutcome, X = subset(tempLearn, select = tempWhichScreen[library$library$rowScreen[s], ], drop=FALSE), newX = subset(tempValid, select = tempWhichScreen[library$library$rowScreen[s], ], drop=FALSE), family = family, id = tempId, obsWeights = tempObsWeights)))
      if(inherits(testAlg, "try-error")) {
        warning(paste("Error in algorithm", library$library$predAlgorithm[s], "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n" ))
        # errorsInCVLibrary[s] <<- 1
      } else {
        out[, s] <- testAlg$pred
        if(saveCVFitLibrary){
          model_out[[s]] <- testAlg$fit
        }
      }
      if (verbose) message(paste("CV", libraryNames[s]))
    } #end library
    if(saveCVFitLibrary){
      names(model_out) <- libraryNames
    }
    invisible(list(out = out, model_out = model_out))
  }
  # the lapply performs the cross-validation steps to create Z
  # additional steps to put things in the correct order
  # rbind unlists the output from lapply
  # need to unlist folds to put the rows back in the correct order
  time_train_start = proc.time()
  
  crossValFUN_out <- lapply(validRows, FUN = .crossValFUN, 
                            Y = Y, dataX = X, id = id, 
                            obsWeights = obsWeights, 
                            library = library, kScreen = kScreen, 
                            k = k, p = p, libraryNames = libraryNames,
                            saveCVFitLibrary = control$saveCVFitLibrary)
  Z[unlist(validRows, use.names = FALSE), ] <- do.call('rbind', lapply(crossValFUN_out, "[[", "out"))
  if(control$saveCVFitLibrary){
    cvFitLibrary <- lapply(crossValFUN_out, "[[", "model_out")
  }else{
    cvFitLibrary <- NULL
  }
  # Check for errors. If any algorithms had errors, replace entire column with
  # 0 even if error is only in one fold.
  errorsInCVLibrary <- apply(Z, 2, function(x) anyNA(x))
  if (sum(errorsInCVLibrary) > 0) {
    Z[, as.logical(errorsInCVLibrary)] <- 0
  }
  if (all(Z == 0)) {
    stop("All algorithms dropped from library")
  }
  
  ########################################################################################################
  # step 2: calculate super learner weights for square loss (MSE)
  print("Calculate super learner weights - MSE")
  
  # using loss functions: Squared Loss
  # get optimum weights for each algorithm
  
  # One-stage Super Learner
  getCoef.LS.one <- method$`Square Loss`$computeCoef(Z = Z, Y = Y, libraryNames = libraryNames,
                                obsWeights = obsWeights, control = control,
                                verbose = verbose,
                                errorsInLibrary = errorsInCVLibrary)
  
  coef.LS.one <- getCoef.LS.one$coef
  
  # One-stage Discrete Learner
  discrete.LS.one <- libraryNames[which.min(onestage.fit$cvRisk)]
  coef.LS.discrete.one <- as.numeric(libraryNames == discrete.LS.one)
  names(coef.LS.discrete.one) <- libraryNames
  
  # Set a default in case the method does not return the optimizer result.
  if (!("optimizer" %in% names(getCoef.LS.one))) {
    getCoef.LS.one["optimizer"] <- NA
  }
  
  ########################################################################################################
  # step 3: calculate super learner weights for Huber loss
  print("Calculate super learner weights for Huber loss")
  
  # Calculate super learner weights for all possible lambda
  
  # Huber One-stage - Different Lambda
  discrete.HUB.one <- NULL
  lamcand.one <- NULL
  coef.HUB.one <- NULL
  coef.HUB.one.discrete <- NULL
  cvRisk.HUBER.one <- NULL
  
  for (i in 1:n.lambda){
    lam <- lambda[i]
    # Super Learner
    getCoef.HUBER <- method$`Huber Loss`$computeCoef(Z=Z, Y=Y, libraryNames=libraryNames, lambda = lam,
                                                     verbose=verbose)
    # Set a default in case the method does not return the optimizer result.
    if (!("optimizer" %in% names(getCoef.HUBER))) {
      getCoef.HUBER["optimizer"] <- NA
    }
    coef.HUBtemp <- getCoef.HUBER$coef
    cvRisk.HUBtemp <- getCoef.HUBER$cvRisk
    # Discrete Learner
    discrete.HUBtemp <- libraryNames[which.min(getCoef.HUBER$cvRisk)]
    coef.discrete.HUBtemp <- as.numeric(libraryNames == discrete.HUBtemp)
    discrete.HUB.one <- c(discrete.HUB.one,discrete.HUBtemp)
    coef.HUB.one <- cbind(coef.HUB.one,coef.HUBtemp)
    coef.HUB.one.discrete <- cbind(coef.HUB.one.discrete,coef.discrete.HUBtemp)
    lamcand.one <- c(lamcand.one,lam)
    cvRisk.HUBER.one <- cbind(cvRisk.HUBER.one,cvRisk.HUBtemp)
  }
  lamname <- paste('lambda =',10000*lambda)                           
  colnames(coef.HUB.one) <- lamname
  rownames(coef.HUB.one) <- libraryNames
  colnames(coef.HUB.one.discrete) <- lamname
  rownames(coef.HUB.one.discrete) <- libraryNames
  discrete.HUBER.one <- data.frame('Lambda'=lamname,'Optimal Learner'=discrete.HUB.one)
  colnames(cvRisk.HUBER.one) <- lamname
  
  
  ########################################################################################################
  ## Huber Loss: Partial-CV 
  ########################################################################################################
  print("Huber partial-CV")
  
  ## Split cv prediction matrix Z into V folds (same split) -- use same v-fold cv in selecting optimal lambda
  ## generate empty matrix to store hold-out predictions
  z.lambda.one.discrete.partial <- matrix(NA, nrow = N, ncol = n.lambda)
  z.lambda.one.partial <- matrix(NA, nrow = N, ncol = n.lambda)
  
  # Only one CV
  for(i in 1:V){
    tempLearn <- Z[-folds[[i]], , drop = FALSE] # training set
    tempOutcome <- Y[-folds[[i]]] # training outcome Y
    tempValid <- Z[folds[[i]], , drop = FALSE] # individial predictions on validation set
    
    # Huber Loss - Different Lambda (just fit on the training set)
    coef.HUB.one <- NULL
    coef.HUB.one.discrete <- NULL
    # different lambda
    for (j in 1:n.lambda){
      lam.1 <- lambda[j]
      # one-stage
      getCoef.HUBER <- method$`Huber Loss`$computeCoef(Z=tempLearn, Y=tempOutcome, libraryNames=libraryNames,
                                                         lambda = lam.1, verbose=verbose)
      # Set a default in case the method does not return the optimizer result.
      if (!("optimizer" %in% names(getCoef.HUBER))) {
        getCoef.HUBER["optimizer"] <- NA
      }
      
      coef.HUBtemp <- getCoef.HUBER$coef
      
      # Discrete Learner -- Huber Loss
      # one-stage
      discrete.HUBtemp <- libraryNames[which.min(getCoef.HUBER$cvRisk)]
      coef.discrete.HUBtemp <- as.numeric(libraryNames == discrete.HUBtemp)
      
      # Super Learner -- Huber Loss
      # one-stage
      coef.HUB.one <- cbind(coef.HUB.one,coef.HUBtemp)
      coef.HUB.one.discrete <- cbind(coef.HUB.one.discrete, coef.discrete.HUBtemp)
    }
    
    # one-stage
    colnames(coef.HUB.one) <- lamname
    rownames(coef.HUB.one) <- libraryNames
    colnames(coef.HUB.one.discrete) <- lamname
    rownames(coef.HUB.one.discrete) <- libraryNames
    
    # Get predicitons on the validation set (from the model just fit on the training set)
    # one-stage discrete learner
    z.lambda.one.discrete.partial[folds[[i]],] = tempValid %*% coef.HUB.one.discrete
    # one-stage super learner
    z.lambda.one.partial[folds[[i]],] = tempValid %*% coef.HUB.one
  }
  
  # one-stage discrete learner
  cvRisk.lambda.one.discrete.partial <- apply(z.lambda.one.discrete.partial, 2, function(x) mean(obsWeights*(x-Y)^2))
  huber.optimal.one.discrete.partial <- lamname[which.min(cvRisk.lambda.one.discrete.partial)]
  coef.huber.optimal.one.discrete.partial <- coef.HUB.one.discrete[,which.min(cvRisk.lambda.one.discrete.partial)]
  cvRisk.huber.optimal.one.discrete.partial <- cvRisk.HUBER.one[,which.min(cvRisk.lambda.one.discrete.partial)]
  
  # one-stage super learner
  cvRisk.lambda.one.partial <- apply(z.lambda.one.partial, 2, function(x) mean(obsWeights*(x-Y)^2))
  huber.optimal.one.partial <- lamname[which.min(cvRisk.lambda.one.partial)]
  coef.huber.optimal.one.partial <- coef.HUB.one[,which.min(cvRisk.lambda.one.partial)]
  cvRisk.huber.optimal.one.partial <- cvRisk.HUBER.one[,which.min(cvRisk.lambda.one.partial)]
  
  
  ########################################################################################################
  ## Huber Loss: Nested-CV 
  ########################################################################################################
  print("Huber nested-CV")
  
  ## nested CV (two cv)
  ## Inner loop CV: calculate optimal weights alpha
  ## Outer loop CV: calculate optimal threshold lambda
  
  ## generate empty matrix to store hold-out predictions
  z.lambda.one.discrete.nested <- matrix(NA, nrow = N, ncol = n.lambda)
  z.lambda.one.nested <- matrix(NA, nrow = N, ncol = n.lambda)
  
  ## Outer CV (fold = V): use same v-fold cv as square loss (MSE)
  for(i in 1:V){
    nest.X <- X[-folds[[i]], , drop = FALSE] # training set X
    nest.Outcome <- Y[-folds[[i]]] # training outcome Y
    nest.Valid <- Z[folds[[i]], , drop = FALSE] # individual predictions on validation set
    
    ## Inner CV (fold = D): create a new cross-validation
    # D-fold cross-validation on training set
    D <- 10 # number of folds
    tempN <- length(nest.Outcome)
    tempid <- seq(1,tempN,1)
    temp.Weights <- rep(1,tempN)
    cvControl.2 <- do.call('SuperLearner.CV.control', list(V=D))
    cvControl.2$validRows = CVFolds(N = tempN, id = tempid, Y = nest.Outcome, cvControl = cvControl.2)
    
    # fit the whole model using one stage option (rather than two stages)
    onestage.fit.temp <- SuperLearner(Y=nest.Outcome, X=nest.X, family=family.single,
                                      SL.library=SL.library, verbose=verbose,
                                      method=method.CC_LS.scale,
                                      control=list(saveCVFitLibrary=F),
                                      cvControl=cvControl.2)
    
    # get the cross-validated predicted values for each algorithm in SL.library
    z.temp <- onestage.fit.temp$Z
    
    # Huber Loss - Different Lambda (just fit on the training set)
    coef.HUB.nest.one <- NULL
    coef.discrete.HUB.nest.one <- NULL
    
    # different lambda
    for (j in 1:n.lambda){
      lam.1 <- lambda[j]
      
      # one-stage
      getCoef.HUBER.nest.one <- method$`Huber Loss`$computeCoef(Z=z.temp, Y=tempOutcome, libraryNames=libraryNames,
                                                                lambda = lam.1, verbose=verbose)
      
      # Set a default in case the method does not return the optimizer result.
      if (!("optimizer" %in% names(getCoef.HUBER.nest.one))) {
        getCoef.HUBER.nest.one["optimizer"] <- NA
      }
      
      # coefficients
      coef.HUBtemp.nest.one <- getCoef.HUBER.nest.one$coef
      
      # Discrete Learner
      # one stage
      discrete.HUBtemp.nest.one <- libraryNames[which.min(getCoef.HUBER.nest.one$cvRisk)]
      coef.discrete.HUBtemp.nest.one <- as.numeric(libraryNames == discrete.HUBtemp.nest.one)
      
      # Super Learner
      coef.HUB.nest.one <- cbind(coef.HUB.nest.one, coef.HUBtemp.nest.one)
      coef.discrete.HUB.nest.one <- cbind(coef.discrete.HUB.nest.one, coef.discrete.HUBtemp.nest.one)
    }
    
    # Get predicitons for the validation set (from the model just fit on the training set)
    # one-stage discrete learner
    z.lambda.one.discrete.nested[folds[[i]],] = nest.Valid %*% coef.discrete.HUB.nest.one
    # one-stage super learner
    z.lambda.one.nested[folds[[i]],] = nest.Valid %*% coef.HUB.nest.one
  }
  
  # one-stage discrete learner
  cvRisk.lambda.one.discrete.nested <- apply(z.lambda.one.discrete.nested, 2, function(x) mean(obsWeights*(x-Y)^2))
  huber.optimal.one.discrete.nested <- lamname[which.min(cvRisk.lambda.one.discrete.nested)]
  coef.huber.optimal.one.discrete.nested <- coef.HUB.one.discrete[,which.min(cvRisk.lambda.one.discrete.nested)]
  cvRisk.huber.optimal.one.discrete.nested <- cvRisk.HUBER.one[,which.min(cvRisk.lambda.one.discrete.nested)]
  
  # one-stage super learner
  cvRisk.lambda.one.nested <- apply(z.lambda.one.nested, 2, function(x) mean(obsWeights*(x-Y)^2))
  huber.optimal.one.nested <- lamname[which.min(cvRisk.lambda.one.nested)]
  coef.huber.optimal.one.nested <- coef.HUB.one[,which.min(cvRisk.lambda.one.nested)]
  cvRisk.huber.optimal.one.nested <- cvRisk.HUBER.one[,which.min(cvRisk.lambda.one.nested)]
  
  time_train = proc.time() - time_train_start
  
  #########################################################################################################
  # step 4: now fit all algorithms in library on entire data set (X) and predict on newX
  print("Generate predictions")
  
  m <- dim(newX)[1L]
  predY <- matrix(NA, nrow = m, ncol = k)
  
  .screenFun <- function(fun, list) {
    screen_fn = get(fun, envir = env)
    testScreen <- try(do.call(screen_fn, list))
    if (inherits(testScreen, "try-error")) {
      warning(paste("replacing failed screening algorithm,", fun, ", with All() in full data", "\n "))
      out <- rep(TRUE, ncol(list$X))
    } else {
      out <- testScreen
    }
    return(out)
  }
  
  time_predict_start = proc.time()
  
  whichScreen <- sapply(library$screenAlgorithm, FUN = .screenFun, list = list(Y = Y, X = X, family = family, id = id, obsWeights = obsWeights), simplify = FALSE)
  whichScreen <- do.call(rbind, whichScreen)
  
  
  .predFun <- function(index, lib, Y, dataX, newX, whichScreen, family, id, obsWeights, verbose, control, libraryNames) {
    pred_fn = get(lib$predAlgorithm[index], envir = env)
    testAlg <- try(do.call(pred_fn, list(Y = Y,
                                         X = subset(dataX,
                                                    select = whichScreen[lib$rowScreen[index], ], drop=FALSE),
                                         newX = subset(newX, select = whichScreen[lib$rowScreen[index], ], drop=FALSE),
                                         family = family, id = id, obsWeights = obsWeights)))
    # testAlg <- try(do.call(lib$predAlgorithm[index], list(Y = Y, X = dataX[, whichScreen[lib$rowScreen[index], drop = FALSE]], newX = newX[, whichScreen[lib$rowScreen[index], drop = FALSE]], family = family, id = id, obsWeights = obsWeights)))
    if (inherits(testAlg, "try-error")) {
      warning(paste("Error in algorithm", lib$predAlgorithm[index], " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n" ))
      out <- rep.int(NA, times = nrow(newX))
    } else {
      out <- testAlg$pred
      if (control$saveFitLibrary) {
        eval(bquote(fitLibrary[[.(index)]] <- .(testAlg$fit)), envir = fitLibEnv)
      }
    }
    if (verbose) {
      message(paste("full", libraryNames[index]))
    }
    invisible(out)
  }
  
  
  predY <- do.call('cbind', lapply(seq(k), FUN = .predFun,
                                   lib = library$library, Y = Y, dataX = X,
                                   newX = newX, whichScreen = whichScreen,
                                   family = family, id = id,
                                   obsWeights = obsWeights, verbose = verbose,
                                   control = control,
                                   libraryNames = libraryNames))
  
  # check for errors
  errorsInLibrary <- apply(predY, 2, function(algorithm) anyNA(algorithm))
  if (sum(errorsInLibrary) > 0) {
    if (sum(coef[as.logical(errorsInLibrary)]) > 0) {
      warning(paste0("Re-running estimation of coefficients removing failed algorithm(s)\n",
                     "Original coefficients are: \n", paste(coef, collapse = ", "), "\n"))
      z[, as.logical(errorsInLibrary)] <- 0
      if (all(z == 0)) {
        stop("All algorithms dropped from library")
      }
      
      getCoef.LS <- method$`Square Loss`$computeCoef(Z = z, Y = Y, libraryNames = wholelibrary,
                                       obsWeights = obsWeights, control = control,
                                       verbose = verbose,
                                       errorsInLibrary = errorsInLibrary)
      coef.LS <- getCoef.LS$coef
      names(coef.LS) <- wholelibrary
      
      for (i in 1:n.lambda){
        lam <- lambda[i]
        getCoef.HUBER <- method$`Huber Loss`$computeCoef(Z=z,Y=Y,libraryNames=wholelibrary,lambda = lam,
                                                       verbose=verbose)
        # Set a default in case the method does not return the optimizer result.
        if (!("optimizer" %in% names(getCoef.HUBER))) {
          getCoef.HUBER["optimizer"] <- NA
        }
        coef.HUBtemp <- getCoef.HUBER$coef
        cvRisk.HUBtemp <- getCoef.HUBER$cvRisk
        # Discrete Learner -- Huber Loss
        discrete.HUBtemp <- wholelibrary[which.min(getCoef.HUBER$cvRisk)]
        discrete.HUB <- c(discrete.HUB,discrete.HUBtemp)
        coef.HUB <- cbind(coef.HUB,coef.HUBtemp)
        lamcand <- c(lamcand,lam)
        cvRisk.HUBER <- cbind(cvRisk.HUBER,cvRisk.HUBtemp)
      }
      colnames(coef.HUB) <- lamname
      rownames(coef.HUB) <- wholelibrary
      discrete.HUBER <- data.frame('Lambda'=lamname,'Optimal Learner'=discrete.HUB)
      colnames(cvRisk.HUBER) <- lamname
      
    } else {
      warning("Coefficients already 0 for all failed algorithm(s)")
    }
  }
  
  ######################## Ensemble predictions ###########################
  # one-stage super learner
  # standard & huber-partial CV & huber-nested CV
  SL.pred.MSE.one <- 
    method$`Square Loss`$computePred(predY = predY, coef = coef.LS.one, control=control)
  SL.pred.Huber.partial.one <- 
    method$`Huber Loss`$computePred(predY = predY, coef = coef.huber.optimal.one.partial, control=control)
  SL.pred.Huber.nested.one <- 
    method$`Huber Loss`$computePred(predY = predY, coef = coef.huber.optimal.one.nested, control=control)
  
  getPred.SL.one <- list('Standard (MSE)' = SL.pred.MSE.one,
                         'Huber-partial CV' = SL.pred.Huber.partial.one,
                         'Huber-nested CV' = SL.pred.Huber.nested.one)
  
  # one-stage discrete learner
  # standard & huber-partial CV & huber-nested CV
  discrete.pred.MSE.one <- 
    method$`Square Loss`$computePred(predY = predY, coef = coef.LS.discrete.one, control=control)
  discrete.pred.Huber.partial.one <- 
    method$`Huber Loss`$computePred(predY = predY, coef = coef.huber.optimal.one.discrete.partial, control=control)
  discrete.pred.Huber.nested.one <- 
    method$`Huber Loss`$computePred(predY = predY, coef = coef.huber.optimal.one.discrete.nested, control=control)
  
  getPred.discrete.one <- list('Standard (MSE)' = discrete.pred.MSE.one,
                               'Huber-partial CV' = discrete.pred.Huber.partial.one,
                               'Huber-nested CV' = discrete.pred.Huber.nested.one)
  
  time_predict = proc.time() - time_predict_start
  
  # Add names of algorithms to the predictions.
  colnames(predY) <- libraryNames
  
  # Clean up when errors in library.
  if(sum(errorsInCVLibrary) > 0) {
    getCoef.LS$cvRisk[as.logical(errorsInCVLibrary)] <- NA
  }
  
  # Finish timing the full SuperLearner execution.
  time_end = proc.time()
  
  # Compile execution times.
  times = list(everything = time_end - time_start,
               train = time_train,
               predict = time_predict)
  
  ## discrete learner
  # one-stage
  # standard & huber-partial CV & huber-nested CV
  discrete.learner.one <- list('Standard (MSE)' = discrete.LS.one,
                               'Huber-partial CV' = 
                                 libraryNames[which.min(cvRisk.huber.optimal.one.discrete.partial)],
                               'Huber-nested CV' = 
                                 libraryNames[which.min(cvRisk.huber.optimal.one.discrete.nested)])
  
  ## super learner weights
  # one-stage
  # standard & huber-partial CV & huber-nested CV
  coef.one <- list('Standard (MSE)' = coef.LS.one,
                   'Huber-partial CV' = coef.huber.optimal.one.partial,
                   'Huber-nested CV' = coef.huber.optimal.one.nested)
  
  ## cvRisk
  # one-stage discrete learner
  # standard & huber-partial CV & huber-nested CV
  cvRisk.one.discrete <- list('Standard (MSE)' = onestage.fit$cvRisk,
                              'Huber-partial CV' = cvRisk.huber.optimal.one.discrete.partial,
                              'Huber-nested CV' = cvRisk.huber.optimal.one.discrete.nested)
  
  # one-stage super learner
  # standard & huber-partial CV & huber-nested CV
  cvRisk.one.SL <- list('Standard (MSE)' = onestage.fit$cvRisk,
                        'Huber-partial CV' = cvRisk.huber.optimal.one.partial,
                        'Huber-nested CV' = cvRisk.huber.optimal.one.nested)
  
  
  # predictions
  getPred.all <- list("one-stage discrete learner" = getPred.discrete.one,
                      "one-stage super learner" = getPred.SL.one
  )
  
  # cvRisk
  cvRisk.all <- list("one-stage discrete learner" = cvRisk.one.discrete,
                     "one-stage super learner" = cvRisk.one.SL
  )
  
  # optimal lambda
  optimal.lambda <- list("one-stage huber discrete learner" = list("partial-CV" = huber.optimal.one.discrete.partial,
                                                                   "nested-CV" = huber.optimal.one.discrete.nested),
                         "one-stage super learner" = list("partial-CV" = huber.optimal.one.partial,
                                                          "nested-CV" = huber.optimal.one.nested)
  )
  
  # Put everything together in a list.
  out <- list(
    call = call,
    libraryNames = libraryNames,
    SL.library = library,
    SL.predict = getPred.all,
    coef = coef.one,
    library.predict = predY,
    Z = Z,
    cvRisk = cvRisk.all,
    optimal.learner = discrete.learner.one,
    huber.optimal.lambda = optimal.lambda,
    family = family,
    fitLibrary = fitLibrary,
    cvfitLibrary = cvfitLibrary,
    varNames = varNames,
    validRows = folds,
    number0 = num.0,
    method = method,
    whichScreen = whichScreen,
    control = control,
    cvControl = cvControl,
    errorsInCVLibrary = errorsInCVLibrary,
    errorsInLibrary = errorsInLibrary,
    env = env,
    times = times
  )
  class(out) <- c("SuperLearner")
  out
}
