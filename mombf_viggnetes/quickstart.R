#initialization
library(mombf)
set.seed(1234)

# define data
x <- matrix(rnorm(100*3),nrow=100,ncol=3)
theta <- matrix(c(1,1,0),ncol=1)
y <- x %*% theta + rnorm(100)

# define model
#Default MOM prior on parameters
priorCoef <- momprior(tau=0.348)
#Beta-Binomial prior for model space
priorDelta <- modelbbprior(1,1)

#Model selection
fit1 <- modelSelection(y ~ x[,1]+x[,2]+x[,3], priorCoef=priorCoef, priorDelta=priorDelta)

#Posterior model probabilities
postProb(fit1)
