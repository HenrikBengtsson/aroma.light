library("aroma.light")

data(cars)
attach(cars)
plot(speed, dist, main = "data(cars)  &  robust smoothing splines")

# Fit a smoothing spline using L_2 norm
cars.spl <- smooth.spline(speed, dist)
lines(cars.spl, col = "blue")

# Fit a smoothing spline using L_1 norm
cars.rspl <- robustSmoothSpline(speed, dist)
lines(cars.rspl, col = "red")

# Fit a smoothing spline using L_2 norm with 10 degrees of freedom
lines(smooth.spline(speed, dist, df=10), lty=2, col = "blue")

# Fit a smoothing spline using L_1 norm with 10 degrees of freedom
lines(robustSmoothSpline(speed, dist, df=10), lty=2, col = "red")

legend(5,120, c(
    paste("smooth.spline [C.V.] => df =",round(cars.spl$df,1)),
    paste("robustSmoothSpline [C.V.] => df =",round(cars.rspl$df,1)),
    "standard with s( * , df = 10)", "robust with s( * , df = 10)"
  ), col = c("blue","red","blue","red"), lty = c(1,1,2,2), bg='bisque')
