library(rstan)

#a)
T <- 200
mu <- 10
sigma_2 <- 2

phi <- seq(-1,1,0.05)

AR_func <- function(mu, phi, x_prev, e_mu, e_sigma_sq) {
	err <- rnorm(1, e_mu, sqrt(e_sigma_sq))
	print(err)
	#err <- 0
	return(mu + phi * (x_prev - mu) + err)
}


calced_val <- matrix(nrow = T,ncol=length(phi))
xt <- mu
#xt <- rnorm(1,mu,)
for (i in 1:length(phi)) {
	for (t in 1:T) {
		xt <- AR_func(mu,phi[i],xt,0,sigma_2)
		calced_val[t,i] <- xt
	}
	temp_mean <- mean(calced_val[,i])
	plot(c(1:T),calced_val[,i],type="b",main=paste("phi =",phi[i],"mean =",temp_mean),ylim=c(4,16))
	Sys.sleep(0.5)
}


i <- 2
plot(c(1:T),calced_val[,i],type="b",main=paste("phi =",phi[i]))


#b)
x_t <- mu
y_t <- mu
x <- c()
y <- c()
for (t in 1:T) {
	x_t <- AR_func(mu,0.3,x_t,0,sigma_2)
	y_t <- AR_func(mu,0.95,y_t,0,sigma_2)
	x <- rbind(x,x_t)
	y <- rbind(y,y_t)
}


StanModel <- 
	'data {
	int<lower=0> N;
	vector[N] x;
	}
  parameters {
  real mu;
  real phi;
  real<lower=0> sigma_2;
  }
  model {
  for (n in 2:N)
    x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma);
  }'


xfit <- stan(model_code=StanModel, data=list(x=x, N=T))
print(xfit)



y = stan(model_code=ARStanModel, data=list(x=x.95, N=T))





options(mc.cores = parallel::detectCores())











y = c(4,5,6,4,0,2,5,3,8,6,10,8) 
N = length(y)
StanModel = ' data {
int<lower=0> N; // Number of observations
int<lower=0> y[N]; // Number of flowers }
parameters { 
real mu;
real<lower=0> sigma2; }
model {
mu ~ normal(0,100); // Normal with mean 0, st.dev. 100
sigma2 ~ scaled_inv_chi_square(1,2); // Scaled-inv-chi2 with nu 1, sigma 2 for(i in 1:N)
y[i] ~ normal(mu,sqrt(sigma2)); }'

data = list(N=N, y=y)
burnin = 1000
niter = 2000
fit = stan(model_code=StanModel,data=data,
					 warmup=burnin,iter=niter,chains=4) # Print the fitted model
print(fit,digits_summary=3) # Extract posterior samples
postDraws <- extract(fit)
# Do traceplots of the first chain
par(mfrow = c(1,1)) plot(postDraws$mu[1:(niter-burnin)],type="l",ylab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit)
# Bivariate posterior plots
pairs(fit)












data <- read.delim("campy.dat", header = TRUE, sep = "\n")[, 1]

