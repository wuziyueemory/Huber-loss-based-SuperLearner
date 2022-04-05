##################################################################################################

## data generating process 
# training data used for ATE estimation 

createData <- function(n, skew){
  # generate covariates (10)
  # A is the treatment variable with bernoulli distribution
  # w2 & w7 have Uniform distribution
  # w3 & w8 have a normal distribution
  # w4 & w9 have a gamma distribution
  # w5 & w10 have a Poisson distribution
  # w6 have bernoulli distribution
  
  A <- rbinom(n,1,0.5)
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
  
  # main term
  main <- A + w2 + w3 + A*w4 + A*w5 + w2*w3 + w4*w5
  if (skew=="low") { # low skewness
    y <- 8300 * rTweedie(mu = 15 + abs(main), p = 1.5, phi = 5)
  } else if (skew=="medium") { # medium skewness
    y <- 1000 * rTweedie(mu = 3.5*main^2, p = 1.2, phi = 3.5)
  } else if (skew=="high") { # high skewness
    y <- 1200 * rTweedie(mu = (main)^2, p = 1.932, phi = 10)
  }
  # generate a dataframe
  return(data.frame(pid=1:n,A=A,w2=w2,w3=w3,w4=w4,w5=w5,w6=w6,
                    w7=w7,w8=w8,w9=w9,w10=w10,y=y))
}






# data used for calculating true ATE (A=1)

createData.1 <- function(n, skew){
  # generate covariates (10)
  # A is the treatment variable with bernoulli distribution
  # w2 & w7 have Uniform distribution
  # w3 & w8 have a normal distribution
  # w4 & w9 have a gamma distribution
  # w5 & w10 have a Poisson distribution
  # w6 have bernoulli distribution
  
  A <- rep(1, n) # A set to 1
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
  
  # main term
  main <- A + w2 + w3 + A*w4 + A*w5 + w2*w3 + w4*w5
  if (skew=="low") { # low skewness
    y <- 8300 * rTweedie(mu = 15 + abs(main), p = 1.5, phi = 5)
  } else if (skew=="medium") { # medium skewness
    y <- 1000 * rTweedie(mu = 3.5*main^2, p = 1.2, phi = 3.5)
  } else if (skew=="high") { # high skewness
    y <- 1200 * rTweedie(mu = (main)^2, p = 1.932, phi = 10)
  }
  # generate a dataframe
  return(data.frame(pid=1:n,A=A,w2=w2,w3=w3,w4=w4,w5=w5,w6=w6,
                    w7=w7,w8=w8,w9=w9,w10=w10,y=y))
}





# data used for calculating true ATE (A=0)

createData.0 <- function(n, skew){
  # generate covariates (10)
  # A is the treatment variable with bernoulli distribution
  # w2 & w7 have Uniform distribution
  # w3 & w8 have a normal distribution
  # w4 & w9 have a gamma distribution
  # w5 & w10 have a Poisson distribution
  # w6 have bernoulli distribution
  
  A <- rep(0, n) # A set to 0
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
  
  # main term
  main <- A + w2 + w3 + A*w4 + A*w5 + w2*w3 + w4*w5
  if (skew=="low") { # low skewness
    y <- 8300 * rTweedie(mu = 15 + abs(main), p = 1.5, phi = 5)
  } else if (skew=="medium") { # medium skewness
    y <- 1000 * rTweedie(mu = 3.5*main^2, p = 1.2, phi = 3.5)
  } else if (skew=="high") { # high skewness
    y <- 1200 * rTweedie(mu = (main)^2, p = 1.932, phi = 10)
  }
  # generate a dataframe
  return(data.frame(pid=1:n,A=A,w2=w2,w3=w3,w4=w4,w5=w5,w6=w6,
                    w7=w7,w8=w8,w9=w9,w10=w10,y=y))
}






