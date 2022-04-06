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
             'Standard (MSE): two-stage super learner' = 
               cbind(Risk = x$cvRisk$`two-stage super learner`$`Standard (MSE)`,
                     Coef = x$coef$`two-stage`$`Standard (MSE)`),
             'Huber-partial CV: one-stage super learner' = 
               cbind(Risk = x$cvRisk$`one-stage super learner`$`Huber-partial CV`,
                     Coef = x$coef$`one-stage`$`Huber-partial CV`),
             'Huber-partial CV: two-stage super learner' = 
               cbind(Risk = x$cvRisk$`two-stage super learner`$`Huber-partial CV`,
                     Coef = x$coef$`two-stage`$`Huber-partial CV`),
             'Huber-nested CV: one-stage super learner' = 
               cbind(Risk = x$cvRisk$`one-stage super learner`$`Huber-nested CV`,
                     Coef = x$coef$`one-stage`$`Huber-nested CV`),
             'Huber-nested CV: two-stage super learner' = 
               cbind(Risk = x$cvRisk$`two-stage super learner`$`Huber-nested CV`,
                     Coef = x$coef$`two-stage`$`Huber-nested CV`))
  )
}







## One-Stage Huber loss-based Super Learner
##############################################################################################################

onestage.HuberSL <- function(Y, X, newX = NULL, library.2stage, library.1stage, lambda,
                             family.1=binomial, family.2=gaussian, family.single=gaussian,
                             id=NULL, verbose=FALSE, control = list(), 
                             cvControl = list(), obsWeights = NULL, env = parent.frame()){
  
  # Begin timing how long two-stage SuperLearner takes to execute
  time_start = proc.time()
  
  # Get details of estimation algorithm for the algorithm weights (coefficients)
  method <- list('Square Loss' = method.CC_LS.scale(),
                 'Huber Loss' = method.CC_HUBER())
  
  # get defaults for controls and make sure in correct format
  control <- do.call('SuperLearner.control', control)
  # change the logical for saveCVFitLibrary to TRUE (we are gonna use that)
  control$saveCVFitLibrary <- TRUE
  
  cvControl <- do.call('SuperLearner.CV.control', cvControl)
  
  # put together the library
  library.stage1 <- library.2stage$stage1
  library.stage2 <- library.2stage$stage2
  library.stage_1 <- SuperLearner:::.createLibrary(library.stage1)
  library.stage_2 <- SuperLearner:::.createLibrary(library.stage2)
  library.stage_single <- SuperLearner:::.createLibrary(library.1stage)
  SuperLearner:::.check.SL.library(library = c(unique(library.stage_1$library$predAlgorithm),
                                               library.stage_1$screenAlgorithm))
  SuperLearner:::.check.SL.library(library = c(unique(library.stage_2$library$predAlgorithm),
                                               library.stage_2$screenAlgorithm))
  SuperLearner:::.check.SL.library(library = c(unique(library.stage_single$library$predAlgorithm),
                                               library.stage_single$screenAlgorithm))
  call <- match.call(expand.dots = TRUE)
  
  # should we be checking X and newX for data.frame?
  # data.frame not required, but most of the built-in wrappers assume a data.frame
  if(!inherits(X, 'data.frame')) message('X is not a data frame. Check the algorithms in SL.library to make sure they are compatible with non data.frame inputs')
  varNames <- colnames(X)
  N <- dim(X)[1L]
  p <- dim(X)[2L]
  k.1 <- nrow(library.stage_1$library)
  k.2 <- nrow(library.stage_2$library)
  k.single <- nrow(library.stage_single$library)
  k.2stage <- k.1*k.2
  k.all <- k.1*k.2+k.single
  kScreen.1 <- length(library.stage_1$screenAlgorithm)
  kScreen.2 <- length(library.stage_2$screenAlgorithm)
  kScreen.single <- length(library.stage_single$screenAlgorithm)
  kscreen.2stage <- kScreen.1*kScreen.2
  kscreen.all <- kScreen.1*kScreen.2+kScreen.single
  n.lambda <- length(lambda)
  
  # family can be either character or function, so these lines put everything together
  #family for stage 1
  if(is.character(family.1))
    family.1 <- get(family.1, mode="function", envir=parent.frame())
  if(is.function(family.1))
    family.1 <- family.1()
  if (is.null(family.1$family)) {
    print(family.1)
    stop("'family' not recognized")
  }
  # family for stage 2
  if(is.character(family.2))
    family.2 <- get(family.2, mode="function", envir=parent.frame())
  if(is.function(family.2))
    family.2 <- family.2()
  if (is.null(family.2$family)) {
    print(family.2)
    stop("'family' not recognized")
  }
  # family for single stage
  if(is.character(family.single))
    family.single <- get(family.single, mode="function", envir=parent.frame())
  if(is.function(family.single))
    family.single <- family.single()
  if (is.null(family.single$family)) {
    print(family.single)
    stop("'family' not recognized")
  }
  
  # check if the model use method.AUC
  if (family.1$family != "binomial" & isTRUE("cvAUC" %in% method$require)){
    stop("'method.AUC' is designed for the 'binomial' family only")
  }
  if (family.2$family != "binomial" & isTRUE("cvAUC" %in% method$require)){
    stop("'method.AUC' is designed for the 'binomial' family only")
  }
  if (family.single$family != "binomial" & isTRUE("cvAUC" %in% method$require)){
    stop("'method.AUC' is designed for the 'binomial' family only")
  }
  
  # chekc whether screen algorithm compatible with number of columns
  if(p < 2 & !identical(library.stage_1$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }
  if(p < 2 & !identical(library.stage_2$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }
  if(p < 2 & !identical(library.stage_single$screenAlgorithm, "All")) {
    warning('Screening algorithms specified in combination with single-column X.')
  }
  
  # generate library names
  # stage 1
  libname.stage.1 <- NULL
  lib.stage1 <- library.stage_1$library$predAlgorithm
  lib.stage1.screen <- library.stage_1$screenAlgorithm[library.stage_1$library$rowScreen]
  repname <- function(x) {
    name <- rep(x,k.2)
    return(name)
  }
  libname.stage.1 <- unlist(lapply(lib.stage1,repname), use.names=FALSE)
  libname.stage.1.screen <- unlist(lapply(lib.stage1.screen,repname), use.names=FALSE)
  # stage 2
  libname.stage.2 <- NULL
  lib.stage2 <- library.stage_2$library$predAlgorithm
  lib.stage2.screen <- library.stage_2$screenAlgorithm[library.stage_2$library$rowScreen]
  libname.stage.2 <- rep(lib.stage2,k.1)
  libname.stage.2.screen <- rep(lib.stage2.screen,k.1)
  # single stage
  libname.stage.single <- library.stage_single$library$predAlgorithm
  libname.stage.single.screen <- library.stage_single$screenAlgorithm[library.stage_single$library$rowScreen]
  
  twostage.library <- paste("S1:",paste(libname.stage.1,libname.stage.1.screen,sep="_"),
                            "+ S2:",paste(libname.stage.2,libname.stage.2.screen,sep="_"))
  wholelibrary <- c(twostage.library,
                    paste("Single:",paste(libname.stage.single,libname.stage.single.screen,sep="_")))
  
  # add family for two stages and single stage
  family <- list(stage1 = family.1,stage2 = family.2,
                 stage.single = family.single)
  
  # add library for two stages
  lib <- list(twostage=data.frame("predAlgorithm"=paste("S1:",libname.stage.1,
                                                        "+ S2:",libname.stage.2),
                                  "rowScreen.Stage.1"=rep(library.stage_1$library$rowScreen,each=k.2),
                                  "rowScreen.Stage.2"=rep(library.stage_2$library$rowScreen,k.1)),
              singlestage=library.stage_single$library)
  library <- list("library"=lib,
                  "screenAlgorithm"=list(stage.1 = library.stage_1$screenAlgorithm,
                                         stage.2 = library.stage_2$screenAlgorithm,
                                         stage.single = library.stage_single$screenAlgorithm))
  
  # if newX is missing, use X
  if(is.null(newX)) {
    newX <- X
  }
  
  # Various chekcs for data structure
  if(!identical(colnames(X), colnames(newX))) {
    stop("The variable names and order in newX must be identical to the variable names and order in X")
  }
  if (sum(is.na(X)) > 0 | sum(is.na(newX)) > 0 | sum(is.na(Y)) > 0) {
    stop("missing data is currently not supported. Check Y, X, and newX for missing values")
  }
  if (!is.numeric(Y)) {
    stop("the outcome Y must be a numeric vector")
  }
  
  # errors records if an algorithm stops either in the CV step and/or in full data
  errorsInCVLibrary <- rep(0, k.all)
  errorsInLibrary <- rep(0, k.all)
  
  ########################################################################################################
  # Step 0: make valid rows
  # ensure each folds have approximately equal number of obs with y=0
  V <- cvControl$V
  ord <- order(Y)
  cvfold <- rep(c(1:V,V:1),N)[1:N]
  folds <- split(ord, factor(cvfold))
  folds <- lapply(folds,sort,decreasing=FALSE)
  # check
  tab <- rep(NA,V)
  for (i in 1:V) {
    tab[i] <- sum(Y[folds[[i]]]==0)
  }
  num.0 <- data.frame("fold"=paste0("fold ",c(1:cvControl$V)),"number.of.0"=tab)
  
  cvControl$validRows = folds
  
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
  
  #########################################################################################################
  # Step 1: fit superlearner for modeling prob of y=0
  print("Fit stage-1: P(Y=0)")
  
  time_train_start = proc.time()
  
  # list all the algorithms considered
  # save cross-validated fits (V) in the control option
  step1.fit <- SuperLearner(Y=as.numeric(Y==0),X=X,family=family.1,
                            SL.library=library.stage1,
                            verbose=verbose,
                            method=method.CC_nloglik,
                            control=list(saveCVFitLibrary=T),
                            cvControl=cvControl)
  
  # get the cross-validated predicted values for each algorithm in SL.library
  # P(Y=0|X)
  z1 <- step1.fit$Z
  
  # get the cross-validated fits (V fits for V training set) for each algorithm
  stage1.cvFitLibrary <- step1.fit$cvFitLibrary
  
  
  ##########################################################################################################
  # step 2: fit model for E[Y|Y>0,X]
  print("Fit stage-2: E[Y|Y>0,X]")
  
  # create function for the cross-validation step at stage 2:
  .crossValFUN <- function(valid, Y, dataX, id, obsWeights, library, family,
                           kScreen, k, p, libraryNames, saveCVFitLibrary) {
    tempLearn <- dataX[-valid, , drop = FALSE]
    tempOutcome <- Y[-valid]
    tempValid <- dataX[valid, , drop = FALSE]
    tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
    tempId <- id[-valid]
    tempObsWeights <- obsWeights[-valid]
    
    # create subset with only obs y>0
    pid <- as.numeric(row.names(tempLearn))
    dat.p <- cbind(pid,tempLearn,tempOutcome)
    tempOutcome.p <- tempOutcome[tempOutcome>0]
    tempLearn.p <- dat.p[dat.p$tempOutcome>0,-c(1,ncol(dat.p))]
    tempId.p <- dat.p[dat.p$tempOutcome>0,1]
    tempObsWeights.p <- obsWeights[tempId.p]
    
    # should this be converted to a lapply also?
    for(s in seq(kScreen)) {
      screen_fn = get(library$screenAlgorithm[s], envir = env)
      testScreen <- try(do.call(screen_fn,
                                list(Y = tempOutcome.p,
                                     X = tempLearn.p,
                                     family = family,
                                     id = tempId.p,
                                     obsWeights = tempObsWeights.p)))
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
      testAlg <- try(do.call(pred_fn,
                             list(Y = tempOutcome.p,
                                  X = subset(tempLearn.p,
                                             select = tempWhichScreen[library$library$rowScreen[s], ],
                                             drop=FALSE),
                                  newX = subset(tempValid,
                                                select = tempWhichScreen[library$library$rowScreen[s], ],
                                                drop=FALSE),
                                  family = family,
                                  id = tempId.p,
                                  obsWeights = tempObsWeights.p)))
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
    invisible(list(out = out, model_out = model_out, pred_id = valid))
  }
  
  # the lapply performs the cross-validation steps to create Z for stage 2
  # additional steps to put things in the correct order
  # rbind unlists the output from lapply
  # need to unlist folds to put the rows back in the correct order
  
  crossValFUN_out <- lapply(cvControl$validRows, FUN = .crossValFUN,
                            Y = Y, dataX = X, id = id,
                            obsWeights = obsWeights, family = family.2,
                            library = library.stage_2, kScreen = kScreen.2,
                            k = k.2, p = p, libraryNames = library.stage_2$library$predAlgorithm,
                            saveCVFitLibrary = control$saveCVFitLibrary)
  
  # create matrix to store results
  z2 <- matrix(NA,nrow = N,ncol=k.2)
  z2[unlist(cvControl$validRows, use.names = FALSE), ] <- do.call('rbind', lapply(crossValFUN_out, "[[", "out"))
  
  if(control$saveCVFitLibrary){
    stage2.cvFitLibrary <- lapply(crossValFUN_out, "[[", "model_out")
  }else{
    stage2.cvFitLibrary <- NULL
  }
  
  # z1 for E[P(Y=0|X)]
  # z2 for E[Y|Y>0,X]
  # multiply (1-z1)*z2 to generate z
  z.two_stage <- NULL
  for (i in 1:k.1){
    for (j in 1:k.2){
      temp <- (1 - z1[,i]) * z2[,j]
      z.two_stage <- cbind(z.two_stage,temp)
    }
  }
  
  ########################################################################################################
  # step 3: fit the whole model using one stage option (rather than two stages)
  print("Fit single-stage: E[Y|X]")
  
  # list all the algorithms considered
  # save cross-validated fits (V) in the control option
  onestage.fit <- SuperLearner(Y=Y, X=X, family=family.single,
                               SL.library=library.1stage,
                               verbose=verbose,
                               method=method.CC_LS.scale,
                               control=list(saveCVFitLibrary=T),
                               cvControl=cvControl)
  
  # get the cross-validated predicted values for each algorithm in SL.library
  z.single <- onestage.fit$Z
  
  # get the cross-validated fits (V fits for V training set) for each algorithm
  single.stage.cvFitLibrary <- onestage.fit$cvFitLibrary
  
  # combine 2 stages output z with 1 stage prediction output z
  z <- cbind(z.two_stage, z.single)
  
  # Check for errors. If any algorithms had errors, replace entire column with
  # 0 even if error is only in one fold.
  errorsInCVLibrary <- apply(z, 2, function(x) anyNA(x))
  if (sum(errorsInCVLibrary) > 0) {
    z[, as.logical(errorsInCVLibrary)] <- 0
  }
  if (all(z == 0)) {
    stop("All algorithms dropped from library")
  }
  
  ########################################################################################################
  # step 4: calculate super learner weights for square loss (MSE)
  print("Calculate super learner weights - MSE")
  
  # using loss functions: Squared Loss
  # get optimum weights for each algorithm
  
  # Two-stage Super Learner
  getCoef.LS.two <- method$`Square Loss`$computeCoef(Z=z,Y=Y,libraryNames=wholelibrary,
                                                     verbose=verbose)
  coef.LS.two <- getCoef.LS.two$coef
  names(coef.LS.two) <- wholelibrary
  
  # Two stage Discrete Learner
  discrete.LS.two <- wholelibrary[which.min(getCoef.LS.two$cvRisk)]
  coef.LS.discrete.two <- as.numeric(wholelibrary == discrete.LS.two)
  names(coef.LS.discrete.two) <- wholelibrary
  
  # One-stage Super Learner
  coef.LS.one <- onestage.fit$coef
  
  # One-stage Discrete Learner
  discrete.LS.one <- onestage.fit$libraryNames[which.min(onestage.fit$cvRisk)]
  coef.LS.discrete.one <- as.numeric(onestage.fit$libraryNames == discrete.LS.one)
  names(coef.LS.discrete.one) <- onestage.fit$libraryNames
  
  # Set a default in case the method does not return the optimizer result.
  if (!("optimizer" %in% names(getCoef.LS.two))) {
    getCoef.LS.two["optimizer"] <- NA
  }
  
  
  ########################################################################################################
  # step 5: calculate super learner weights for Huber loss
  print("Calculate super learner weights for Huber loss")
  
  # Calculate super learner weights for all possible lambda
  
  # Huber Two-stage - Different Lambda
  discrete.HUB <- NULL
  lamcand <- NULL
  coef.HUB <- NULL
  coef.discrete.HUB <- NULL
  cvRisk.HUBER <- NULL
  
  for (i in 1:n.lambda){
    lam <- lambda[i]
    # Super Learner
    getCoef.HUBER <- method$`Huber Loss`$computeCoef(Z=z,Y=Y,libraryNames=wholelibrary,lambda = lam,
                                                     verbose=verbose)
    # Set a default in case the method does not return the optimizer result.
    if (!("optimizer" %in% names(getCoef.HUBER))) {
      getCoef.HUBER["optimizer"] <- NA
    }
    coef.HUBtemp <- getCoef.HUBER$coef
    cvRisk.HUBtemp <- getCoef.HUBER$cvRisk
    # Discrete Learner
    discrete.HUBtemp <- wholelibrary[which.min(getCoef.HUBER$cvRisk)]
    coef.discrete.HUBtemp <- as.numeric(wholelibrary == discrete.HUBtemp)
    discrete.HUB <- c(discrete.HUB,discrete.HUBtemp)
    coef.HUB <- cbind(coef.HUB,coef.HUBtemp)
    coef.discrete.HUB <- cbind(coef.discrete.HUB,coef.discrete.HUBtemp)
    lamcand <- c(lamcand,lam)
    cvRisk.HUBER <- cbind(cvRisk.HUBER,cvRisk.HUBtemp)
  }
  lamname <- paste('lambda =',10000*lamcand)
  colnames(coef.HUB) <- lamname
  rownames(coef.HUB) <- wholelibrary
  colnames(coef.discrete.HUB) <- lamname
  rownames(coef.discrete.HUB) <- wholelibrary
  discrete.HUBER <- data.frame('Lambda'=lamname,'Optimal Learner'=discrete.HUB)
  colnames(cvRisk.HUBER) <- lamname
  
  
  # Huber One-stage - Different Lambda
  libname <- onestage.fit$libraryNames
  discrete.HUB.one <- NULL
  lamcand.one <- NULL
  coef.HUB.one <- NULL
  coef.HUB.one.discrete <- NULL
  cvRisk.HUBER.one <- NULL
  
  for (i in 1:n.lambda){
    lam <- lambda[i]
    # Super Learner
    getCoef.HUBER <- method$`Huber Loss`$computeCoef(Z=z.single,Y=Y,libraryNames=libname,lambda = lam,
                                                     verbose=verbose)
    # Set a default in case the method does not return the optimizer result.
    if (!("optimizer" %in% names(getCoef.HUBER))) {
      getCoef.HUBER["optimizer"] <- NA
    }
    coef.HUBtemp <- getCoef.HUBER$coef
    cvRisk.HUBtemp <- getCoef.HUBER$cvRisk
    # Discrete Learner
    discrete.HUBtemp <- libname[which.min(getCoef.HUBER$cvRisk)]
    coef.discrete.HUBtemp <- as.numeric(libname == discrete.HUBtemp)
    discrete.HUB.one <- c(discrete.HUB.one,discrete.HUBtemp)
    coef.HUB.one <- cbind(coef.HUB.one,coef.HUBtemp)
    coef.HUB.one.discrete <- cbind(coef.HUB.one.discrete,coef.discrete.HUBtemp)
    lamcand.one <- c(lamcand.one,lam)
    cvRisk.HUBER.one <- cbind(cvRisk.HUBER.one,cvRisk.HUBtemp)
  }
  colnames(coef.HUB.one) <- lamname
  rownames(coef.HUB.one) <- libname
  colnames(coef.HUB.one.discrete) <- lamname
  rownames(coef.HUB.one.discrete) <- libname
  discrete.HUBER.one <- data.frame('Lambda'=lamname,'Optimal Learner'=discrete.HUB.one)
  colnames(cvRisk.HUBER.one) <- lamname
  
  
  ########################################################################################################
  ## Huber Loss: Partial-CV 
  ########################################################################################################
  print("Huber partial-CV")
  
  ## Split cv prediction matrix Z into V folds (same split) -- use same v-fold cv in selecting optimal lambda
  ## generate empty matrix to store hold-out predictions
  z.lambda.two.discrete.partial <- matrix(NA, nrow = N, ncol = n.lambda)
  z.lambda.two.partial <- matrix(NA, nrow = N, ncol = n.lambda)
  z.lambda.one.discrete.partial <- matrix(NA, nrow = N, ncol = n.lambda)
  z.lambda.one.partial <- matrix(NA, nrow = N, ncol = n.lambda)
  
  # Only one CV
  for(i in 1:V){
    tempLearn <- z[-folds[[i]], , drop = FALSE] # training set (two)
    tempLearn.single <- z.single[-folds[[i]], , drop = FALSE] # training set (one)
    tempOutcome <- Y[-folds[[i]]] # training outcome Y
    tempValid <- z[folds[[i]], , drop = FALSE] # individial predictions on validation set (two)
    tempValid.single <- z.single[folds[[i]], , drop = FALSE] # individial predictions on validation set (one)
    
    # Huber Loss - Different Lambda (just fit on the training set)
    coef.HUB.1 <- NULL
    coef.discrete.HUB.1 <- NULL
    libname <- onestage.fit$libraryNames
    coef.HUB.one.2 <- NULL
    coef.HUB.one.discrete.2 <- NULL
    # different lambda
    for (j in 1:n.lambda){
      lam.1 <- lambda[j]
      # two-stage
      getCoef.HUBER.1 <- method$`Huber Loss`$computeCoef(Z=tempLearn,Y=tempOutcome,libraryNames=wholelibrary,
                                                         lambda = lam.1,verbose=verbose)
      # one-stage
      getCoef.HUBER.2 <- method$`Huber Loss`$computeCoef(Z=tempLearn.single,Y=tempOutcome,libraryNames=libname,
                                                         lambda = lam.1,verbose=verbose)
      # Set a default in case the method does not return the optimizer result.
      if (!("optimizer" %in% names(getCoef.HUBER.1))) {
        getCoef.HUBER.1["optimizer"] <- NA
      }
      if (!("optimizer" %in% names(getCoef.HUBER.2))) {
        getCoef.HUBER.2["optimizer"] <- NA
      }
      
      coef.HUBtemp.1 <- getCoef.HUBER.1$coef
      coef.HUBtemp.2 <- getCoef.HUBER.2$coef
      
      # Discrete Learner -- Huber Loss
      # two-stage
      discrete.HUBtemp.1 <- wholelibrary[which.min(getCoef.HUBER.1$cvRisk)]
      coef.discrete.HUBtemp.1 <- as.numeric(wholelibrary == discrete.HUBtemp.1)
      # one-stage
      discrete.HUBtemp.2 <- libname[which.min(getCoef.HUBER.2$cvRisk)]
      coef.discrete.HUBtemp.2 <- as.numeric(libname == discrete.HUBtemp.2)
      
      # Super Learner -- Huber Loss
      # two-stage
      coef.HUB.1 <- cbind(coef.HUB.1,coef.HUBtemp.1)
      coef.discrete.HUB.1 <- cbind(coef.discrete.HUB.1,coef.discrete.HUBtemp.1)
      # one-stage
      coef.HUB.one.2 <- cbind(coef.HUB.one.2,coef.HUBtemp.2)
      coef.HUB.one.discrete.2 <- cbind(coef.HUB.one.discrete.2,coef.discrete.HUBtemp.2)
    }
    
    # two-stage
    colnames(coef.HUB.1) <- lamname
    rownames(coef.HUB.1) <- wholelibrary
    colnames(coef.discrete.HUB.1) <- lamname
    rownames(coef.discrete.HUB.1) <- wholelibrary
    # one-stage
    colnames(coef.HUB.one.2) <- lamname
    rownames(coef.HUB.one.2) <- libname
    colnames(coef.HUB.one.discrete.2) <- lamname
    rownames(coef.HUB.one.discrete.2) <- libname
    
    # Get predicitons on the validation set (from the model just fit on the training set)
    # two-stage discrete learner
    z.lambda.two.discrete.partial[folds[[i]],] = tempValid %*% coef.discrete.HUB.1
    # two-stage super learner
    z.lambda.two.partial[folds[[i]],] = tempValid %*% coef.HUB.1
    # one-stage discrete learner
    z.lambda.one.discrete.partial[folds[[i]],] = tempValid.single %*% coef.HUB.one.discrete.2
    # one-stage super learner
    z.lambda.one.partial[folds[[i]],] = tempValid.single %*% coef.HUB.one.2
  }
  
  
  # two-stage discrete learner
  cvRisk.lambda.two.discrete.partial <- apply(z.lambda.two.discrete.partial, 2, function(x) mean(obsWeights*(x-Y)^2))
  huber.optimal.two.discrete.partial <- lamname[which.min(cvRisk.lambda.two.discrete.partial)]
  coef.huber.optimal.two.discrete.partial <- coef.discrete.HUB[,which.min(cvRisk.lambda.two.discrete.partial)]
  cvRisk.huber.optimal.two.discrete.partial <- cvRisk.HUBER[,which.min(cvRisk.lambda.two.discrete.partial)]
  
  # two-stage super learner
  cvRisk.lambda.two.partial <- apply(z.lambda.two.partial, 2, function(x) mean(obsWeights*(x-Y)^2))
  huber.optimal.two.partial <- lamname[which.min(cvRisk.lambda.two.partial)]
  coef.huber.optimal.two.partial <- coef.HUB[,which.min(cvRisk.lambda.two.partial)]
  cvRisk.huber.optimal.two.partial <- cvRisk.HUBER[,which.min(cvRisk.lambda.two.partial)]
  
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
  z.lambda.two.discrete.nested <- matrix(NA, nrow = N, ncol = n.lambda)
  z.lambda.two.nested <- matrix(NA, nrow = N, ncol = n.lambda)
  z.lambda.one.discrete.nested <- matrix(NA, nrow = N, ncol = n.lambda)
  z.lambda.one.nested <- matrix(NA, nrow = N, ncol = n.lambda)
  
  ## Outer CV (fold = V): use same v-fold cv as square loss (MSE)
  for(i in 1:V){
    nest.X <- X[-folds[[i]], , drop = FALSE] # training set X
    nest.Outcome <- Y[-folds[[i]]] # training outcome Y
    nest.Valid <- z[folds[[i]], , drop = FALSE] # individual predictions on validation set (two)
    nest.Valid.single <- z.single[folds[[i]], , drop = FALSE] # individual predictions on validation set (one)
    
    ## Inner CV (fold = D): create a new cross-validation
    # D-fold cross-validation on training set
    D <- 5 # number of folds
    tempN <- length(nest.Outcome)
    tempid <- seq(1,tempN,1)
    temp.Weights <- rep(1,tempN)
    cvControl.2 <- do.call('SuperLearner.CV.control', list(V=D))
    # tempord <- order(nest.Outcome)
    # tempcvfold <- rep(c(1:V,V:1),N)[1:tempN]
    # tempfolds <- split(tempord, factor(tempcvfold))
    # tempfolds <- lapply(tempfolds,sort,decreasing=FALSE)
    cvControl.2$validRows = CVFolds(N = tempN, id = tempid, Y = nest.Outcome, cvControl = cvControl.2)
    
    # Step 1: fit superlearner for modeling Prob(Y=0|X)
    step1.fit.temp <- SuperLearner(Y=as.numeric(nest.Outcome==0),X=nest.X,family=family.1,
                                   SL.library=library.stage1,verbose=verbose,
                                   method=method.CC_nloglik,
                                   control=list(saveCVFitLibrary=F),
                                   cvControl=cvControl.2)
    
    # get the cross-validated predicted values for each algorithm in SL.library
    # P(Y=0|X)
    z1.temp <- step1.fit.temp$Z
    
    # step 2: fit model for E[Y|Y>0,X]
    crossValFUN_out.temp <- lapply(cvControl.2$validRows, FUN = .crossValFUN,
                                   Y = nest.Outcome, dataX = nest.X, id = tempid,
                                   obsWeights = NULL, family = family.2,
                                   library = library.stage_2, kScreen = kScreen.2,
                                   k = k.2, p = p, libraryNames = library.stage_2$library$predAlgorithm,
                                   saveCVFitLibrary = FALSE)
    
    # create matrix to store results
    z2.temp <- matrix(NA,nrow = tempN,ncol=k.2)
    z2.temp[unlist(cvControl.2$validRows, use.names = FALSE), ] <- do.call('rbind', lapply(crossValFUN_out.temp, "[[", "out"))
    
    # z1.temp for P(Y=0|X)
    # z2.temp for E[Y|Y>0,X]
    # multiply (1-z1.temp)*z2.temp to generate z.temp
    z.two_stage.temp <- NULL
    for (a in 1:k.1){
      for (b in 1:k.2){
        temp <- (1 - z1.temp[,a]) * z2.temp[,b]
        z.two_stage.temp <- cbind(z.two_stage.temp,temp)
      }
    }
    
    # step 3: fit the whole model using one stage option (rather than two stages)
    onestage.fit.temp <- SuperLearner(Y=nest.Outcome,X=nest.X,family=family.single,
                                      SL.library=library.1stage,verbose=verbose,
                                      method=method.CC_LS.scale,
                                      control=list(saveCVFitLibrary=F),
                                      cvControl=cvControl.2)
    
    # get the cross-validated predicted values for each algorithm in SL.library
    z.single.temp <- onestage.fit.temp$Z
    
    # combine 2 stages output z with 1 stage prediction output z
    z.temp <- cbind(z.two_stage.temp,z.single.temp)
    
    # Huber Loss - Different Lambda (just fit on the training set)
    coef.HUB.nest.two <- NULL
    coef.discrete.HUB.nest.two <- NULL
    libname <- onestage.fit$libraryNames
    coef.HUB.nest.one <- NULL
    coef.discrete.HUB.nest.one <- NULL
    
    # different lambda
    for (j in 1:n.lambda){
      lam.1 <- lambda[j]
      
      # two-stage
      getCoef.HUBER.nest.two <- method$`Huber Loss`$computeCoef(Z=z.temp,Y=tempOutcome,libraryNames=wholelibrary,
                                                                lambda = lam.1,verbose=verbose)
      # one-stage
      getCoef.HUBER.nest.one <- method$`Huber Loss`$computeCoef(Z=z.single.temp,Y=tempOutcome,libraryNames=libname,
                                                                lambda = lam.1,verbose=verbose)
      
      # Set a default in case the method does not return the optimizer result.
      if (!("optimizer" %in% names(getCoef.HUBER.nest.two))) {
        getCoef.HUBER.nest.two["optimizer"] <- NA
      }
      if (!("optimizer" %in% names(getCoef.HUBER.nest.one))) {
        getCoef.HUBER.nest.one["optimizer"] <- NA
      }
      
      # coefficients
      coef.HUBtemp.nest.two <- getCoef.HUBER.nest.two$coef
      coef.HUBtemp.nest.one <- getCoef.HUBER.nest.one$coef
      
      # Discrete Learner
      # two stage
      discrete.HUBtemp.nest.two <- wholelibrary[which.min(getCoef.HUBER.nest.two$cvRisk)]
      coef.discrete.HUBtemp.nest.two <- as.numeric(wholelibrary == discrete.HUBtemp.nest.two)
      # one stage
      discrete.HUBtemp.nest.one <- libname[which.min(getCoef.HUBER.nest.one$cvRisk)]
      coef.discrete.HUBtemp.nest.one <- as.numeric(libname == discrete.HUBtemp.nest.one)
      
      # Super Learner
      coef.HUB.nest.two <- cbind(coef.HUB.nest.two, coef.HUBtemp.nest.two)
      coef.discrete.HUB.nest.two <- cbind(coef.discrete.HUB.nest.two, coef.discrete.HUBtemp.nest.two)
      coef.HUB.nest.one <- cbind(coef.HUB.nest.one, coef.HUBtemp.nest.one)
      coef.discrete.HUB.nest.one <- cbind(coef.discrete.HUB.nest.one, coef.discrete.HUBtemp.nest.one)
    }
    
    # Get predicitons for the validation set (from the model just fit on the training set)
    # two-stage discrete learner
    z.lambda.two.discrete.nested[folds[[i]],] = nest.Valid %*% coef.discrete.HUB.nest.two
    # two-stage super learner
    z.lambda.two.nested[folds[[i]],] = nest.Valid %*% coef.HUB.nest.two
    # one-stage discrete learner
    z.lambda.one.discrete.nested[folds[[i]],] = nest.Valid.single %*% coef.discrete.HUB.nest.one
    # one-stage super learner
    z.lambda.one.nested[folds[[i]],] = nest.Valid.single %*% coef.HUB.nest.one
  }
  
  
  # two-stage discrete learner
  cvRisk.lambda.two.discrete.nested <- apply(z.lambda.two.discrete.nested, 2, function(x) mean(obsWeights*(x-Y)^2))
  huber.optimal.two.discrete.nested <- lamname[which.min(cvRisk.lambda.two.discrete.nested)]
  coef.huber.optimal.two.discrete.nested <- coef.discrete.HUB[,which.min(cvRisk.lambda.two.discrete.nested)]
  cvRisk.huber.optimal.two.discrete.nested <- cvRisk.HUBER[,which.min(cvRisk.lambda.two.discrete.nested)]
  
  # two-stage super learner
  cvRisk.lambda.two.nested <- apply(z.lambda.two.nested, 2, function(x) mean(obsWeights*(x-Y)^2))
  huber.optimal.two.nested <- lamname[which.min(cvRisk.lambda.two.nested)]
  coef.huber.optimal.two.nested <- coef.HUB[,which.min(cvRisk.lambda.two.nested)]
  cvRisk.huber.optimal.two.nested <- cvRisk.HUBER[,which.min(cvRisk.lambda.two.nested)]
  
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
  # step 6: now fit all algorithms in library on entire data set (X) and predict on newX
  print("Generate predictions")
  
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
  
  # stage 1
  whichScreen.stage1 <- sapply(library$screenAlgorithm$stage.1, FUN = .screenFun,
                               list = list(Y = Y, X = X, family = family, id = id, obsWeights = NULL),
                               simplify = FALSE)
  whichScreen.stage1 <- do.call(rbind, whichScreen.stage1)
  # stage 2
  whichScreen.stage2 <- sapply(library$screenAlgorithm$stage.2, FUN = .screenFun,
                               list = list(Y = Y, X = X, family = family, id = id, obsWeights = NULL),
                               simplify = FALSE)
  whichScreen.stage2 <- do.call(rbind, whichScreen.stage2)
  # single stage
  whichScreen.stage.single <- sapply(library$screenAlgorithm$stage.single, FUN = .screenFun,
                                     list = list(Y = Y, X = X, family = family, id = id, obsWeights = NULL),
                                     simplify = FALSE)
  whichScreen.stage.single <- do.call(rbind, whichScreen.stage.single)
  # combine together
  whichScreen <- list(stage1 = whichScreen.stage1,
                      stage2 = whichScreen.stage2,
                      single.stage = whichScreen.stage.single)
  
  # Prediction for each algorithm
  .predFun <- function(index, lib, Y, dataX, newX, whichScreen, family, id, obsWeights,
                       verbose, control, libraryNames) {
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
  
  ######################### stage 1 #############################
  # put fitLibrary at stage 1 in it's own environment to locate later
  fitLibEnv <- new.env()
  assign('fitLibrary', vector('list', length = k.1), envir = fitLibEnv)
  assign('libraryNames', library.stage_1$library$predAlgorithm, envir = fitLibEnv)
  evalq(names(fitLibrary) <- library.stage_1$library$predAlgorithm, envir = fitLibEnv)
  # get prediction for stage 1
  predY.stage1 <- do.call('cbind', lapply(seq(k.1), FUN = .predFun,
                                          lib = library.stage_1$library, 
                                          Y = as.numeric(Y==0), 
                                          dataX = X,
                                          newX = newX, 
                                          whichScreen = whichScreen$stage1,
                                          family = family.1, id = id,
                                          obsWeights = obsWeights, verbose = verbose,
                                          control = control,
                                          libraryNames = library.stage_1$library$predAlgorithm))
  
  # save fit library for stage 1
  stage1.fitlib <- get("fitLibrary",envir = fitLibEnv)
  
  ######################### stage 2 #############################
  # put fitLibrary at stage 2 in it's own environment to locate later
  fitLibEnv <- new.env()
  assign('fitLibrary', vector('list', length = k.2), envir = fitLibEnv)
  assign('libraryNames', library.stage_2$library$predAlgorithm, envir = fitLibEnv)
  evalq(names(fitLibrary) <- library.stage_2$library$predAlgorithm, envir = fitLibEnv)
  
  # get prediction for stage 2
  # create subset with only obs y>0
  pid <- c(1:N)
  dat.p <- cbind(pid,X,Y)
  Y.p <- Y[Y>0]
  X.p <- dat.p[dat.p$Y>0,-c(1,ncol(dat.p))]
  p.id <- dat.p[dat.p$Y>0,1]
  p.obsWeights <- obsWeights[p.id]
  
  predY.stage2 <- do.call('cbind', lapply(seq(k.2), FUN = .predFun,
                                          lib = library.stage_2$library, Y = Y.p, dataX = X.p,
                                          newX = newX, whichScreen = whichScreen$stage2,
                                          family = family.2, id = p.id,
                                          obsWeights = p.obsWeights, verbose = verbose,
                                          control = control,
                                          libraryNames = library.stage_2$library$predAlgorithm))
  
  # save fit library for stage 2
  stage2.fitlib <- get("fitLibrary",envir = fitLibEnv)
  
  
  
  ######################### single stage #############################
  # put fitLibrary at single in it's own environment to locate later
  fitLibEnv <- new.env()
  assign('fitLibrary', vector('list', length = k.single), envir = fitLibEnv)
  assign('libraryNames', library.stage_single$library$predAlgorithm, envir = fitLibEnv)
  evalq(names(fitLibrary) <- library.stage_single$library$predAlgorithm, envir = fitLibEnv)
  
  # save fit library for single stage
  predY.stage.single <- do.call('cbind', lapply(seq(k.single), FUN = .predFun,
                                                lib = library.stage_single$library, Y = Y, dataX = X,
                                                newX = newX, whichScreen = whichScreen$single.stage,
                                                family = family.single, id = id,
                                                obsWeights = obsWeights, verbose = verbose,
                                                control = control,
                                                libraryNames = library.stage_single$library$predAlgorithm))
  
  # save fit library for single stage
  stage.single.fitlib <- get("fitLibrary",envir = fitLibEnv)
  
  ######################## All individual predictions ###########################
  # get prediction for 2-stage model
  predY.twostage <- NULL
  for (i in 1:k.1){
    for (j in 1:k.2){
      pred <- (1 - predY.stage1[,i]) * predY.stage2[,j]
      predY.twostage <- cbind(predY.twostage, pred)
    }
  }
  
  # combine with prediction from singe-stage model
  predY <- cbind(predY.twostage, predY.stage.single)
  
  #generate Fitlibrary
  fitLibrary = list("stage.1"=stage1.fitlib,
                    "stage.2"=stage2.fitlib,
                    "stage.sinlge"=stage.single.fitlib)
  
  #generate cross-validation Fitlibrary
  cvfitLibrary <- list("stage.1"=stage1.cvFitLibrary,
                       "stage.2"=stage2.cvFitLibrary,
                       "stage.single"=single.stage.cvFitLibrary)
  
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
      
      lancand <- NULL
      coef.HUB <- NULL
      discrete.HUB <- NULL
      cvRisk.HUBER <- NULL
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
        # Discrete Learner
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
  # two-stage super learner
  # standard & huber-partial CV & huber-nested CV
  SL.pred.MSE.two <- 
    method$`Square Loss`$computePred(predY = predY, coef = coef.LS.two, control=control)
  SL.pred.Huber.partial.two <- 
    method$`Huber Loss`$computePred(predY = predY, coef = coef.huber.optimal.two.partial, control=control)
  SL.pred.Huber.nested.two <- 
    method$`Huber Loss`$computePred(predY = predY, coef = coef.huber.optimal.two.nested, control=control)
  
  getPred.SL.two <- list('Standard (MSE)' = SL.pred.MSE.two,
                         'Huber-partial CV' = SL.pred.Huber.partial.two,
                         'Huber-nested CV' = SL.pred.Huber.nested.two)
  
  # two-stage discrete learner
  # standard & huber-partial CV & huber-nested CV
  discrete.pred.MSE.two <- 
    method$`Square Loss`$computePred(predY = predY, coef = coef.LS.discrete.two, control=control)
  discrete.pred.Huber.partial.two <- 
    method$`Huber Loss`$computePred(predY = predY, coef = coef.huber.optimal.two.discrete.partial, control=control)
  discrete.pred.Huber.nested.two <- 
    method$`Huber Loss`$computePred(predY = predY, coef = coef.huber.optimal.two.discrete.nested, control=control)
  
  getPred.discrete.two <- list('Standard (MSE)' = discrete.pred.MSE.two,
                               'Huber-partial CV' = discrete.pred.Huber.partial.two,
                               'Huber-nested CV' = discrete.pred.Huber.nested.two)
  
  # one-stage super learner
  # standard & huber-partial CV & huber-nested CV
  SL.pred.MSE.one <- 
    method$`Square Loss`$computePred(predY = predY.stage.single, coef = coef.LS.one, control=control)
  SL.pred.Huber.partial.one <- 
    method$`Huber Loss`$computePred(predY = predY.stage.single, coef = coef.huber.optimal.one.partial, control=control)
  SL.pred.Huber.nested.one <- 
    method$`Huber Loss`$computePred(predY = predY.stage.single, coef = coef.huber.optimal.one.nested, control=control)
  
  getPred.SL.one <- list('Standard (MSE)' = SL.pred.MSE.one,
                         'Huber-partial CV' = SL.pred.Huber.partial.one,
                         'Huber-nested CV' = SL.pred.Huber.nested.one)
  
  # one-stage discrete learner
  # standard & huber-partial CV & huber-nested CV
  discrete.pred.MSE.one <- 
    method$`Square Loss`$computePred(predY = predY.stage.single, coef = coef.LS.discrete.one, control=control)
  discrete.pred.Huber.partial.one <- 
    method$`Huber Loss`$computePred(predY = predY.stage.single, coef = coef.huber.optimal.one.discrete.partial, control=control)
  discrete.pred.Huber.nested.one <- 
    method$`Huber Loss`$computePred(predY = predY.stage.single, coef = coef.huber.optimal.one.discrete.nested, control=control)
  
  getPred.discrete.one <- list('Standard (MSE)' = discrete.pred.MSE.one,
                               'Huber-partial CV' = discrete.pred.Huber.partial.one,
                               'Huber-nested CV' = discrete.pred.Huber.nested.one)
  
  time_predict = proc.time() - time_predict_start
  
  # Add names of algorithms to the predictions.
  colnames(predY) <- wholelibrary
  
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
  
  # number of algorithms used in each stage
  library.num <- list(stage1 = k.1,
                      stage2 = k.2,
                      stage.single = k.single)
  
  # original library for each stage's algorithm
  orig.library <- list(stage1 = library.stage_1,
                       stage2 = library.stage_2,
                       stage.single = library.stage_single)
  
  
  ## discrete learner
  # two-stage
  # standard & huber-partial CV & huber-nested CV
  discrete.learner.two <- list('Standard (MSE)' = discrete.LS.two,
                               'Huber-partial CV' = 
                                 wholelibrary[which.min(cvRisk.huber.optimal.two.discrete.partial)],
                               'Huber-nested CV' = 
                                 wholelibrary[which.min(cvRisk.huber.optimal.two.discrete.nested)])
  
  # one-stage
  # standard & huber-partial CV & huber-nested CV
  discrete.learner.one <- list('Standard (MSE)' = discrete.LS.one,
                               'Huber-partial CV' = 
                                 libname[which.min(cvRisk.huber.optimal.one.discrete.partial)],
                               'Huber-nested CV' = 
                                 libname[which.min(cvRisk.huber.optimal.one.discrete.nested)])
  
  
  ## super learner weights
  # two-stage
  # standard & huber-partial CV & huber-nested CV
  coef.two <- list('Standard (MSE)' = coef.LS.two,
                   'Huber-partial CV' = coef.huber.optimal.two.partial,
                   'Huber-nested CV' = coef.huber.optimal.two.nested)
  
  # one-stage
  # standard & huber-partial CV & huber-nested CV
  coef.one <- list('Standard (MSE)' = coef.LS.one,
                   'Huber-partial CV' = coef.huber.optimal.one.partial,
                   'Huber-nested CV' = coef.huber.optimal.one.nested)
  
  ## cvRisk
  # two-stage discrete learner
  # standard & huber-partial CV & huber-nested CV
  cvRisk.two.discrete <- list('Standard (MSE)' = getCoef.LS.two$cvRisk,
                              'Huber-partial CV' = cvRisk.huber.optimal.two.discrete.partial,
                              'Huber-nested CV' = cvRisk.huber.optimal.two.discrete.nested)
  
  # two-stage super learner
  # standard & huber-partial CV & huber-nested CV
  cvRisk.two.SL <- list('Standard (MSE)' = getCoef.LS.two$cvRisk,
                        'Huber-partial CV' = cvRisk.huber.optimal.two.partial,
                        'Huber-nested CV' = cvRisk.huber.optimal.two.nested)
  
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
  
  # create list to include all results
  library.name <- list("one-stage" = libname,
                       "two-stage" = wholelibrary)
  
  # predictions
  getPred.all <- list("one-stage discrete learner" = getPred.discrete.one,
                      "one-stage super learner" = getPred.SL.one,
                      "two-stage discrete learner" = getPred.discrete.two,
                      "two-stage super learner" = getPred.SL.two
  )
  
  # super learner weights
  coef.all <- list("one-stage" = coef.one,
                   "two-stage" = coef.two)
  
  # discrete learner
  discrete.learner.all <- list("one-stage" = discrete.learner.one,
                               "two-stage" = discrete.learner.two)  
  
  # cvRisk
  cvRisk.all <- list("one-stage discrete learner" = cvRisk.one.discrete,
                     "one-stage super learner" = cvRisk.one.SL,
                     "two-stage discrete learner" = cvRisk.two.discrete,
                     "two-stage super learner" = cvRisk.two.SL
  )
  
  # optimal lambda
  optimal.lambda <- list("one-stage huber discrete learner" = list("partial-CV" = huber.optimal.one.discrete.partial,
                                                                   "nested-CV" = huber.optimal.one.discrete.nested),
                         "one-stage super learner" = list("partial-CV" = huber.optimal.one.partial,
                                                          "nested-CV" = huber.optimal.one.nested),
                         "two-stage discrete learner" = list("partial-CV" = huber.optimal.two.discrete.partial,
                                                             "nested-CV" = huber.optimal.two.discrete.nested),
                         "two-stage super learner" = list("partial-CV" = huber.optimal.two.partial,
                                                          "nested-CV" = huber.optimal.two.nested)
  )
  
  # Put everything together in a list.
  out <- list(
    call = call,
    libraryNames = library.name,
    library.Num = library.num,
    orig.library = orig.library,
    SL.library = library,
    SL.predict = getPred.all,
    coef = coef.all,
    library.predict = predY,
    Z = z,
    cvRisk = cvRisk.all,
    optimal.learner = discrete.learner.all,
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
