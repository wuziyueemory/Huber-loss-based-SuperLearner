
# Individual algorithms used in super learner library

#=========================================================#

# (1) GLM (Gamma distribution + Identity-link)

#=========================================================#

SL.gammaIdentityGLM <- function(Y, X, newX, family, obsWeights,...){
  if(family$family=="gaussian"){
    fit.glm <- glm(Y ~ ., data=X, family=Gamma(link='identity'), 
                   weights=obsWeights,
                   control=list(maxit=1000), start=c(mean(Y),rep(0,ncol(X))))
    pred <- predict(fit.glm, newdata=newX, type="response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm"
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.gammaIdentityGLM not written for binomial family")
  }
}


#=========================================================#

# (2) log-OLS: OLS on ln(y) + smear retransformation
# GLM with Gaussian family and id link on log(Y) + Duan (1983) correction

#=========================================================#

SL.logOLS.smear <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    logY <- log(Y)
    fit.logGLM <- glm(logY ~ ., data=X, family=family, weights=obsWeights)
    mu <- predict(fit.logGLM, type="response", newdata=X)
    resid <- logY - mu
    pred <- exp(predict(fit.logGLM, type="response",newdata=newX))*mean(exp(resid))
    fit <- list(object=fit.logGLM, mean(exp(resid)))
    class(fit) <- "SL.logOLS.smear"
  }else{
    stop("SL.logGLM.smear not written for binomial family")
  }
  out <- list(fit=fit, pred=pred)
  return(out)
}

# predict function for SL.logOLS.smear
predict.SL.logOLS.smear <- function(object, newdata, ...){
  mu <- predict(object$object, newdata=newdata, type="response")
  correction <- object[[2]]
  return(exp(mu)*correction) 
}



#=========================================================#

# (3) GLM (Gamma distribution + log-link)

#=========================================================#

SL.gammaLogGLM <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    fit.glm <- glm(Y ~ ., data=X, family=Gamma(link='log'), weights=obsWeights,
                   control=list(maxit=100))
    pred <- predict(fit.glm, newdata=newX, type="response")
    fit <- list(object = fit.glm)
    class(fit) <- "SL.glm" # can use predict.SL.glm
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.logGLM not written for binomial family")
  }
}

