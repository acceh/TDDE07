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
