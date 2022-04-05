# estimators
#==================================================#

# (1)	Zero-Inflated Poisson Model

#==================================================#

SL.zip <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    # round outcome Y to be interger
    Y.int <- round(Y)
    suppressWarnings(
      fit.zip <- zeroinfl(Y.int ~ . | ., data=X,weights = obsWeights)
    )
    pred <- predict(fit.zip, newdata=newX, type="response")
    fit <- list(object = fit.zip)
    class(fit) <- "SL.glm" # can use predict.SL.glm
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.zip not written for binomial family")
  }
}



#==================================================#

# (2)	Tweedie Model

#==================================================#
# ? or use zcpglm (zero-inflated version)

SL.tweedie <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    # using optimizer bobyqa
    suppressWarnings(
      fit.tweedie <-  cpglm(Y~.,data=X,optimizer = "bobyqa")
    )
    pred <- predict(fit.tweedie, newdata=newX, type="response")
    fit <- list(object = fit.tweedie)
    class(fit) <- "SL.tweedie" 
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.tweedie not written for binomial family")
  }
}

predict.SL.tweedie <- function(object, newdata, ...) {
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  pred <- predict(object = object$object, newdata = newdata, type = "response")
  pred
}




#==================================================#

# (3)	Tobit Model

#==================================================#

SL.tobit <- function(Y, X, newX, family, obsWeights, ...){
  if(family$family=="gaussian"){
    suppressWarnings(
      fit.tobit <- vglm(Y ~., tobit(Lower = 0,type.fitted = "censored"),data=X,maxit=100)
    )
    pred <- predict(fit.tobit, newdata=newX, type="response")
    # in case generate negative prediction
    pred[pred<0]=0
    fit <- list(object = fit.tobit)
    class(fit) <- "SL.tobit"
    out <- list(pred=pred, fit=fit)
    return(out)
  }else{
    stop("SL.tobit not written for binomial family")
  }
}

predict.SL.tobit <- function(object, newdata, ...) {
  # newdata must be a dataframe, not a matrix.
  if (is.matrix(newdata)) {
    newdata = as.data.frame(newdata)
  }
  pred <- predict(object = object$object, newdata = newdata, type = "response")
  # in case generate negative prediction
  pred[pred<0]=0
  pred
}
