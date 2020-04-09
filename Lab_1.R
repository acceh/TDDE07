library(manipulate)

BetaPlot <- function(a,b){
	xGrid <- seq(0.001, 0.999, by=0.001)
	prior = dbeta(xGrid, a, b)
	maxDensity <- max(prior) # Use to make the y-axis high enough
	plot(xGrid, prior, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "theta", 
			 ylab = 'Density', main = 'Beta(a,b) density')
}

manipulate(
	BetaPlot(a,b),
	a = slider(1, 10, step=1, initial = 2, label = "The hyperparameter a in Beta(a,b) density"),
	b = slider(1, 10, step=1, initial = 2, label = "The hyperparameter b in Beta(a,b) density")
)

####################################################################
## Plotting the prior-to-posterior mapping for the Bernoulli model.
####################################################################




BetaPriorPostPlot <- function(a,b,n,p){
	xGrid <- seq(0.001, 0.999, by=0.001)
	normalizedLikelihood = dbeta(xGrid, n*p+1, n*(1-p)+1)
	prior = dbeta(xGrid, a, b)
	posterior = dbeta(xGrid, a+n*p, b+n*(1-p))
	maxDensity <- max(normalizedLikelihood, prior, posterior) # Use to make the y-axis high enough
	plot(xGrid, normalizedLikelihood, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = "theta", 
			 ylab = 'Density', main = 'Bernoulli model - Beta(a,b) prior')
	lines(xGrid, posterior, lwd = 3, col = "red")
	lines(xGrid, prior, lwd = 3, col = "green")
	legend(x = 0.01, y = maxDensity*0.95, legend = c("Likelihood (normalized)", "Prior", "Posterior"), col = c("blue","green","red"), lwd = c(3,3,3), cex = 0.7)
}

manipulate(
	BetaPriorPostPlot(a,b,n,p),
	a = slider(1, 100, step=1, initial = 2, label = "The hyperparameter a in Beta(a,b) prior"),
	b = slider(1, 100, step=1, initial = 2, label = "The hyperparameter b in Beta(a,b) prior"),
	n = slider(1, 1000, step=1, initial = 20, label = "The number of trials, n"),
	p = slider(0, 1, step=0.01, initial = 0.4, label = "Success proportion in the sample")
)

BetaPriorPostPlot(1,1,20,0.25)


# Equal tail interval
# 
G <- gini
equal_tail_interval = quantile(G, probs=c(0.05, 0.95))

# Highest Posterior Density
gd = density(G)
# Order x and y by y from high to low
ordered_x = gd$x[order(-gd$y)]
ordered_y = gd$y[order(-gd$y)]

# Iterate until 95% of prob. is captured
prob_mass = 0
total_mass = sum(gd$y)
for(i in 1:length(gd$y)){
	prob_mass = prob_mass + ordered_y[i]
	if(prob_mass / total_mass >= 0.90){
		break
	}
}



# Calculate the interval
# 
 grid_w = 6
grid_h = 5
a = min(ordered_x[1:i])
b = max(ordered_x[1:i])
highest_posterior_density = c(a, b)

plot(gd, col="red", lwd=2, main="2.3 Credibility intervals")
lines(equal_tail_interval, rep(0.12, 2), col="black", lwd=3)
lines(highest_posterior_density, rep(0.02, 2), col="gray", lwd=3)
legend("topright", 
			 legend = c("ETI","HPD"),
			 fill = c("black", "gray"),
			 inset = 0.02)
