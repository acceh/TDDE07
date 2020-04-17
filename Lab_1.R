library(ggplot2)
set.seed(12345)
#1 a)
a0 = 2
b0 = 2
s = 5
f = 15
x_vector = seq(1, 10000, 1)
set.seed(12345)
distribution =rbeta(x_vector, a0+s, b0+f)

y_calced_mean = c(1:length(x_vector))
for (i  in x_vector) {
	y_calced_mean[i] = mean(distribution[0:i])
}

qplot(x=x_vector,y=y_calced_mean, geom="line")+
	labs(x="Number of draws", y="Mean", title="Mean as a function of number of draws")


y_calced_sd = c(1:length(x_vector))
for (i  in x_vector) {
	y_calced_sd[i] = sd(distribution[0:i])
}

qplot(x=x_vector,y=y_calced_sd, geom="line")+
	labs(x="Number of draws", y="Standard deviation", title="Standard deviation as a function of number of draws")

#1 b)

classified = ifelse(distribution > 0.3, 1, 0)
classified_table = table(classified)
counted_classified = classified_table[names(classified_table)==1]

classified_probabilty = counted_classified/length(distribution)
classified_probabilty

true_probability_value = 1- pbeta(0.3, a0+s, b0+f)

#1 c)

log_odds = log(distribution/(1-distribution))
hist(log_odds)

density(log_odds)

#2 a)

library(LaplacesDemon)
observations <- c(44, 25, 45, 52, 30, 63, 19, 50, 34, 69)
n <- length(observations)
my <- 3.7
number_of_draws <- 10000
set.seed(12345)


tausquared <-
	sum((log(observations) - my) ^ 2) / length(observations)

sigma_sq <- rinvchisq(number_of_draws, n, tausquared)

interval = seq(min(sigma_sq), max(sigma_sq), 0.001)
invchisq = dinvchisq(interval, n, tausquared)


plot <-
	ggplot() + geom_histogram(aes(x = sigma_sq, y = ..density..), fill = "gray", binwidth = 0.02) + geom_line(aes(x =interval, y = invchisq), position = "stack") + labs(x="Sigma Squared", y = "Density", title="Draws from posterior against the theoretical posterior distributon.")

plot

#2 b)

gini <- 2*pnorm(sqrt(sigma_sq)/sqrt(2))

plot_gini <- ggplot() + geom_histogram(aes(x=gini, y=..density..), binwidth = 0.001)
plot_gini

#2 c)

library(HDInterval)

interval_pin <- p.interval(gini, HPD=FALSE, prob=0.9)

gini_density <- density(gini)

interval_hdi <- hdi(gini_density, credMass=0.9)



plot_CI <- plot_gini + geom_vline(aes(xintercept=c(interval_pin[1,1], interval_pin[1,2]),color="ETI")) + 
	geom_vline(aes(xintercept=c(interval_hdi["lower"], interval_hdi["upper"]),color="HDI")) + scale_color_manual(name = '',
																																																							 values = c(
																																																							 	'ETI' = 'red',
																																																							 	'HDI' = 'green'
																																																							 ))
plot_CI

#3 a)

directions <-
	c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

mu <- 2.39

p_y_given_k_mu <- c()

k_v = seq(0.001,10,0.01)

posterior <- c()
likelihood <- c()
prior <- c()

for (o in 1:length(k_v)) {
	probs <- c()
	k <- k_v[o]
	for (i in 1:length(directions)) {
		y <- directions[i]
		
		probs[i] <- exp(k * cos(y - mu)) / (2 * pi * besselI(k, 0)) #likelihood
	}
	likelihood[o] <- prod(probs)
	posterior[o] <- prod(probs) * dexp(k,1) #product of likelihood for y1...yn * prior for k
	prior[o] <- dexp(k,1)
}

plotposterior <-
	ggplot() + geom_line((aes(
		x = k_v,
		y = posterior / sum(posterior),
		color = "Posterior"
	))) + geom_line((
		aes(
			x = k_v,
			y = likelihood / sum(likelihood),
			color = "Likelihood"
		)
	))  + geom_line((aes(
		x = k_v,
		y = prior / sum(prior),
		color = "Prior"
	))) + labs(x = "k", y = "Density", title = "Posterior distribution of k") + scale_color_manual(
		name = '',
		values = c(
			'Prior' = 'red',
			'Posterior' = 'green',
			'Likelihood' = 'blue'
		)
	)
plotposterior

#3 b)

plotkmode <- plotposterior + geom_vline(aes(xintercept=k_v[which.max(posterior)]))
plotkmode

print(k_v[which.max(posterior)])