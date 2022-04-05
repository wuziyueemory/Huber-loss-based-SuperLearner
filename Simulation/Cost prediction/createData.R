# data generating process

createData <- function(n,skew){
  # generate covariates (10)
  # w1 & w6 have bernoulli distribution
  # w2 & w7 have Uniform distribution
  # w3 & w8 have a normal distribution
  # w4 & w9 have a gamma distribution
  # w5 & w10 have a Poisson distribution
  w1 <- rbinom(n,1,0.5)
  w2 <- runif(n,0,1)
  w3 <- rnorm(n,0,1)
  w4 <- rgamma(n,1,1)
  w5 <- rpois(n,1)
  w6 <- rbinom(n,1,0.2)
  w7 <- runif(n,-1,1)
  w8 <- rnorm(n,0,3)
  w9 <- rgamma(n,0.5,1)
  w10 <- rpois(n,2)
  y <- rep(NA,n)
  
  # zero proportion - 35%
  # probability of y=0
  prob <- plogis(0.6 + 0.1*(w1 + w2 - w3 + w4 - w5 + w1*w2 - w2*w3 + w3*w4 - w4*w5))
  
  g <- rbinom(n,1,prob)
  # assign g=0 costs
  ind <- g==0
  y[ind] <- 0
  
  # skewness
  # assign g=1 costs
  ind <- g==1
  main <- (w1 + w2 + w3 + w4 + w5 + w1*w2 + w2*w3 + w3*w4 + w4*w5)
  if (skew=="low") { # low skewness
    y[ind] <- rgamma(sum(ind),shape = 10*abs(main[ind]), scale = 1.5)
    y <- rescale(y,to=c(0,1000000))
  } else if (skew=="medium") { # medium skewness
    if (n==250) { # n=250
      y[ind] <- rgamma(sum(ind),shape = 10*abs(main[ind]), scale = 1.5)
      id <- which(y > quantile(y, probs = c(0.75)))
      y[id] <- y[id] + rgamma(length(id),shape = 1.13*(main[id])^2, scale = 1.5)
      y <- rescale(y,to=c(0,1000000))
    } else 
      y[ind] <- rgamma(sum(ind),shape = 10*abs(main[ind]), scale = 1.5)
    id <- which(y > quantile(y, probs = c(0.75)))
    y[id] <- y[id] + rgamma(length(id),shape = 0.71*(main[id])^2, scale = 1.5)
    y <- rescale(y,to=c(0,1000000))
  } else if (skew=="high") { # high skewness
    if (n==250) { # n=250
      y[ind] <- rgamma(sum(ind),shape = 10*abs(main[ind]), scale = 1.5)
      id <- which(y > quantile(y, probs = c(0.75)))
      y[id] <- y[id] + rgamma(length(id),shape = 38*(main[id])^2, scale = 1.5)
      y <- rescale(y,to=c(0,1000000))
    } else 
      y[ind] <- rgamma(sum(ind),shape = 10*abs(main[ind]), scale = 1.5)
    id <- which(y > quantile(y, probs = c(0.75)))
    y[id] <- y[id] + rgamma(length(id),shape = 2.9*(main[id])^2, scale = 1.5)
    y <- rescale(y,to=c(0,1000000))
  }
  # generate a dataframe
  return(data.frame(pid=1:n,w1=w1,w2=w2,w3=w3,w4=w4,w5=w5,w6=w6,
                    w7=w7,w8=w8,w9=w9,w10=w10,y=y))
}
