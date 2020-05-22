library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
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


# b)

## i)

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

## Creates a stan model for given dist and function
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
  x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma_2);
}'


xfit <- stan(model_code=StanModel, data=list(x=c(x), N=T))
sum_xfit <- summary(xfit)

print(sum_xfit$summary[,"n_eff"][-4])
print(sum_xfit$summary[,"mean"][-4])
print(sum_xfit$summary[,"2.5%"][-4])
print(sum_xfit$summary[,"97.5%"][-4])


yfit <- stan(model_code=StanModel, data=list(x=c(y), N=T))
sum_yfit <- summary(yfit)
print(sum_yfit$summary[,"n_eff"][-4])
print(sum_yfit$summary[,"mean"][-4])
print(sum_yfit$summary[,"2.5%"][-4])
print(sum_yfit$summary[,"97.5%"][-4])

## ii)

params_xfit_mu <- extract(xfit)$mu
params_xfit_phi <- extract(xfit)$phi

plot(params_xfit_mu, params_xfit_phi)

# The plot above shows the joint posterior of $\phi$ and $\mu$ for data based on $\phi=0.3$. As one can see there is a high density in the middle of the cluster of the sampled values.
# The middle of the cluster is the mean of each of the paramaters, which can be seen above in i).

params_yfit_mu <- extract(yfit)$mu
params_yfit_phi <- extract(yfit)$phi

plot(params_yfit_mu, params_yfit_phi)

# The plot above shows the joint posterior of $\phi$ and $\mu$ for data based on $\phi=0.95$. Here one can see that the sampled values of $\phi$ is more uniformly distributed between about 1 and 0.9
# while the majority of the sampled values for $\mu$ are all almost the same value, except a small portion of large outliers.


# c)

data <- read.delim("campy.dat", header = TRUE, sep = "\n")[, 1]


PoisStanModel <- 
'data {
int<lower=0> N;
int c[N];
}
parameters {
real mu;
real phi;
real<lower=0> sigma_2;
vector[N] x;
}
model {


for (n in 2:N)
  x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma_2);

for (n in 1:N)
  c[n] ~ poisson(exp(x[n]));
}'

poisfit <- stan(model_code=PoisStanModel, data=list(c=c(data), N=length(c(data))))
sum_poisfit <- summary(poisfit)$summary[-c(1,2,3,144),]

poisfit_upper <- sum_poisfit[,"97.5%"]
poisfit_lower <- sum_poisfit[,"2.5%"]
poisfit_mean <- sum_poisfit[,"mean"]

plot(data)
lines(exp(poisfit_upper), col="red")
lines(exp(poisfit_lower), col="red")
lines(exp(poisfit_mean), col="blue")

# d)


PoisStanModel2 <- 
  'data {
int<lower=0> N;
int c[N];

}
parameters {
real mu;
real phi;
real<lower=0> sigma_2;
vector[N] x;
}
model {

sigma_2 ~scaled_inv_chi_square(N, 0.1);
for (n in 2:N)
  x[n] ~ normal(mu + phi * (x[n-1] - mu), sigma_2);

for (n in 1:N)
  c[n] ~ poisson(exp(x[n]));
}'

poisfit2 <- stan(model_code=PoisStanModel2, data=list(c=c(data), N=length(c(data))))
sum_poisfit2 <- summary(poisfit2)$summary[-c(1,2,3,144),]



poisfit_upper2 <- sum_poisfit2[,"97.5%"]
poisfit_lower2 <- sum_poisfit2[,"2.5%"]
poisfit_mean2 <- sum_poisfit2[,"mean"]

plot(data)
lines(exp(poisfit_upper2), col="red")
lines(exp(poisfit_lower2), col="red")
lines(exp(poisfit_mean2), col="blue")

# The posterior has changed as it is smoother and outliers doesn't affect as much as it did before.


