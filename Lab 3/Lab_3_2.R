setwd("~/Programming/TDDE07/Lab 3")
data <- read.table("eBayNumberOfBidderData.dat", header = TRUE)

# a)

library(glmnet)

model <- glm(nBids ~ . - Const, family = poisson, data = data)

#The covariates that affects the B more is MinBidShare, Sealed and VerifyID.


# b)
library(mvtnorm)
X <- as.matrix(data[, -1])
y <- as.vector(data$nBids)
XTX_inv <- solve(t(X) %*% X)
mu = rep(0, ncol(X))
initBeta <-
	as.vector(rmvnorm(1, mean = rep(0, nrow(XTX_inv)), sigma = (100 * XTX_inv)))


LogPostPoisson <- function(betaVect, y, X, mu, Sigma) {
	linPred <-  X %*% betaVect
	lambda <- exp(linPred)
	
	logLik <- sum(-log(factorial(y)) + y * linPred - lambda)
	
	if (abs(logLik) == Inf)
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
		initBeta,
		LogPostPoisson,
		gr = NULL,
		y,
		X,
		mu,
		100 * XTX_inv,
		method = c("BFGS"),
		control = list(fnscale = -1),
		hessian = TRUE
	)


#Results 2 b)
inverse_hessian <- solve(-OptimResults$hessian)
B_tilde <- OptimResults$par


# c)


sigmaInput <- inverse_hessian
initalThetaInput <- c(rep(0, nrow(sigmaInput)))
c_Input <- 0.5
noIterationsInput <- 5000



rwmSampler <-
	function(logPostFunc,
					 sigma,
					 initalTheta,
					 c,
					 noIterations,
					 ...) {
		outSteps <- matrix(nrow = noIterations, ncol = length(initalTheta))
		outSteps[1, ] <- initalTheta
		
		for (i in c(2:noIterations)) {
			nextTheta <- as.vector(rmvnorm(1, outSteps[i - 1, ], c * sigma))
			r <-
				exp(logPostFunc(nextTheta, ...) - logPostFunc(outSteps[i - 1, ], ...))
			
			
			if (r > runif(1)) {
				outSteps[i, ] <- nextTheta
			} else{
				outSteps[i, ] <- outSteps[i - 1, ]
			}
		}
		return(outSteps)
	}



steps <-
	rwmSampler(
		LogPostPoisson,
		sigmaInput,
		initalThetaInput,
		c_Input,
		noIterationsInput,
		y,
		X,
		mu,
		100 * XTX_inv
	)


for (i in c(1:9)) {
	hist(steps[, i], main = paste("posterior beta num", i))
	
}

for (i in c(1:9)) {
	plot(steps[, i],
			 type = "l",
			 main = paste("Steps for posterior beta num", i))
	
}



#d)

b_tilde <- c()
features <- c(1, 1, 1, 1, 0, 0, 0, 1, 0.5)
for (i in c(1:9)) {
	b_tilde[i] <- median(steps[, i])
}
lambda_tilde <- exp(t(features) %*% b_tilde)

predicted_draws <- rpois(10000, lambda_tilde)

hist(predicted_draws, main = "Predicitve distribution")

prob_for_zero_bidders <-
	sum(ifelse(predicted_draws == 0, 1, 0)) / (length(predicted_draws))
