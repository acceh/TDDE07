setwd("~/Programming/TDDE07/Lab 3")
data <-read.table("eBayNumberOfBidderData.dat", header=TRUE)

# a)

library(glmnet)

model <- glm(nBids ~ . - Const,family=poisson, data=data)

#The covariates that affects the B more is MinBidShare, Sealed and VerifyID.


# b)
library(mvtnorm)
X <- as.matrix(data[,-1])
y <- as.vector(data$nBids)
XTX_inv <- solve(t(X)%*%X)
mu = rep(0,ncol(X))
initVal <- rmvnorm(1,mean=rep(0,nrow(XTX_inv)),sigma=(100 * XTX_inv))


LogPostLogistic <- function(betaVect, y, X, mu, Sigma) {
	nPara <- length(betaVect)
	print("------")
	print(betaVect)
	print(y)
	print(X)
	print(mu)
	print(Sigma)
	print("------")
	
	linPred <- X %*% betaVect
	#print(betaVect)
	# evaluating the log-likelihood
	#logLik <-  sum(linPred * y - log(1 + exp(linPred)))
	
	logLik <-  sum(y*log(linPred)-linPred-log(factorial(y)))
	
	print(logLik)
	
	if (abs(logLik) == Inf || is.nan(logLik)) 
		logLik = -20000
	# Likelihood is not finite, stear the optimizer away from here!
	
	# evaluating the prior
	logPrior <- dmvnorm(betaVect, mu, Sigma, log = TRUE)
	
	#print(logLik + logPrior)
	#print(betaVect)
	# add the log prior and log-likelihood together to get log posterior
	return(logLik + logPrior)
}


OptimResults <-
	optim(
		initVal,
		LogPostLogistic,
		gr = NULL,
		y,
		X,
		mu,
		100 * XTX_inv,
		method = c("BFGS"),
		control = list(fnscale = -1),
		hessian = TRUE
	)

inverse_hessian <- solve(OptimResults$hessian)
