setwd("~/Programming/TDDE07/Lab 2")
data <-
	read.delim("TempLinkoping.txt",
						 header = TRUE,
						 sep = "\t",
						 dec = ".")

library(LaplacesDemon)
library(mvtnorm)
set.seed(12345)
my_0 <- matrix(c(-10, 130, -130), nrow = 3, ncol = 1)

omega_0 <- 3 * diag(3)

v_0 <- 10
sigma_sq_0 <- 2

number_of_draws <- 100

sigma_sq <- rinvchisq(number_of_draws, v_0, sigma_sq_0)

omega_0.inv <- solve(omega_0)


betas <- c()


for (sigma in sigma_sq) {
	beta <- rmvnorm(1, my_0, sigma * omega_0.inv)
	#tmp <- cbind(sigma, beta)
	betas <- rbind(betas, beta)
}


linear_matrix <- c()

times <- seq(0, 1, 0.01)
for (rownumber in 1:nrow(betas)) {
	tmprow <- c()
	for (time in times) {
		tmp <- betas[rownumber, 1] + betas[rownumber, 2] * time + betas[rownumber, 3] *
			time ^ 2
		tmprow <- cbind(tmprow, tmp)
	}
	linear_matrix <- rbind(linear_matrix, tmprow)
}

plot(times,linear_matrix[1,], type='l', ylim=c(-20,30))
for (i in 2:nrow(linear_matrix)) {
	lines(times, linear_matrix[i,])
	
}


################################### Ass 2 #####################################################
library("mvtnorm") # This command reads the mvtnorm package into R's memory. NOW we can use dmvnorm function.

data<-read.table("WomenWork.dat",header=TRUE)  # Spam data from Hastie et al.

y <- as.vector(data$Work)
X <- as.matrix(data[,2:9])


tau <- 10
mu <- rep(0, 8)
sigma <- tau^2 * diag(8)

initVal <- rmvnorm(1, mu, sigma)



LogPostLogistic <- function(betaVect,y,X,mu,Sigma){
	nPara <- length(betaVect);
	linPred <- X%*%betaVect;
	print(betaVect)
	# evaluating the log-likelihood                                    
	logLik <- sum( linPred*y -log(1 + exp(linPred)));
	if (abs(logLik) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
	
	# evaluating the prior
	logPrior <- dmvnorm(betaVect, matrix(0,nPara,1), Sigma, log=TRUE);
	print(logLik + logPrior)
	print(betaVect)
	# add the log prior and log-likelihood together to get log posterior
	return(logLik + logPrior)
}


OptimResults<-optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,sigma,method=c("BFGS"),control=list(fnscale=-1),hessian=TRUE)

inverse_hessian <- solve(OptimResults$hessian)

postModeNSmallChild <- OptimResults$par[7]
stdPostNSmallChild <- sqrt(diag(-inverse_hessian)[7])



betaGrid <-
	seq(
		postModeNSmallChild - 4 * stdPostNSmallChild,
		postModeNSmallChild + 4 * stdPostNSmallChild,
		length = 1000
	)
ciBetaNsmallChild <-
	quantile(rnorm(10000, mean = postModeNSmallChild, sd = stdPostNSmallChild),
						 c(0.025, 0.975))
ggplot() + geom_line(aes(
	x = betaGrid,
	y = dnorm(x = betaGrid, mean = postModeNSmallChild, sd = stdPostNSmallChild)
)) + geom_vline(xintercept = ciBetaNsmallChild[1]) + geom_vline(xintercept = ciBetaNsmallChild[2])


glmModel <- glm(Work ~ 0 + ., data = data, family = binomial)

glmSmallChildVar <- glmModel[["coefficients"]][["NSmallChild"]]

#b)

beta_posterior <- OptimResults$par

features <- c(1, 10.0, 8,(8/10)^2 , 10, 40, 1, 1)

prob_working_women <- (exp(t(features)%*%t(beta_posterior))) / (1+exp(t(features)%*%t(beta_posterior)))

func_working_women_prob <- function(feature, betaposterior) {
	tmp <- (exp(t(feature) %*% betaposterior)) / (1 + exp(t(feature) %*% betaposterior))
	return(tmp)
}

stdbeta <- diag(-inverse_hessian)

sim_beta <- rmvnorm(10000, mean=beta_posterior, sigma=-solve(OptimResults$hessian))

tmp <- c()
for(i in 1:nrow(sim_beta)) {

	tmp <- rbind(tmp,(func_working_women_prob(features, sim_beta[i,])))
}

barplot(tmp)

tmp_list <- ifelse(tmp>0.5, 1, 0)



