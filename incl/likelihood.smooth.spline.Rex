# Define f(x)
f <- expression(0.1*x^4 + 1*x^3 + 2*x^2 + x + 10*sin(2*x))

# Simulate data from this function in the range [a,b]
a <- -2; b <- 5
x <- seq(a, b, length=3000)
y <- eval(f)

# Add some noise to the data
y <- y + rnorm(length(y), 0, 10)

# Plot the function and its second derivative
plot(x,y, type="l", lwd=4)

# Fit a cubic smoothing spline and plot it
g <- smooth.spline(x,y, df=16)
lines(g, col="yellow", lwd=2, lty=2)

# Calculating the (log) likelihood of the fitted spline
l <- likelihood(g)

cat("Log likelihood with unique x values:\n")
print(l)

# Note that this is not the same as the log likelihood of the
# data on the fitted spline iff the x values are non-unique
x[1:5] <- x[1]  # Non-unique x values
g <- smooth.spline(x,y, df=16)
l <- likelihood(g)

cat("\nLog likelihood of the *spline* data set:\n");
print(l)

# In cases with non unique x values one has to proceed as
# below if one want to get the log likelihood for the original
# data.
l <- likelihood(g, x=x, y=y)
cat("\nLog likelihood of the *original* data set:\n");
print(l)





