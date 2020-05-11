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
B_prior <- rmvnorm(1,mean=rep(0,nrow(XTX_inv)),sigma=(100 * XTX_inv))



LogPostLogistic <- function(betaVect, y, X, mu, Sigma) {
	nPara <- length(betaVect)
	
	linPred <- X %*% betaVect
	
	print(betaVect)
	# evaluating the log-likelihood
	logLik <- sum(linPred * y - log(1 + exp(linPred)))
	
	if (abs(logLik) == Inf)
		logLik = -20000
	# Likelihood is not finite, stear the optimizer away from here!
	
	# evaluating the prior
	logPrior <- dmvnorm(betaVect, matrix(0, nPara, 1), Sigma, log = TRUE)
	
	print(logLik + logPrior)
	print(betaVect)
	# add the log prior and log-likelihood together to get log posterior
	return(logLik + logPrior)
}


OptimResults <-
	optim(
		B_prior,
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
