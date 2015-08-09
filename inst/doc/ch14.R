## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(boot)  # boot

data("CH06FI05", package = "ALSM")
data("CH14TA01", package = "ALSM")
data("CH14TA02", package = "ALSM")
data("CH14TA03", package = "ALSM")
data("APPENC11", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH14TA01
names(.data) <- c("x", "y", "fitted")

## ------------------------------------------------------------------------
fit <- glm(y ~ x, family = binomial(link = "logit"), .data)
cbind(
  'Experience'   = .data$x,
  'Task Success' = .data$y,
  'Fitted Value' = .data$fitted,
  'Residual'     = resid(fit))

summary(fit)  # Estimates and Standard Deviations
exp(coef(fit)[2])  # Estimated Odds Ratio
exp(15 * coef(fit)[2])  # Estimated Odds Ratio for 15 Months Training (p 567)

# Plot Logistic Regression and Loess Fit
xx <- with(.data, seq(min(x), max(x), len = 200))
plot(y ~ x, .data, pch = 19, col = "gray40", xlab = "Experience (X)", ylab = "Fitted Value")
lines(xx, predict(loess(y ~ x, .data), data.frame(x = xx)), lty = 2, col = 'blue')
lines(xx, predict(fit, data.frame(x = xx), type = "resp"), lwd = 2)
title("Scatter Plot with Loess (blue) and Logistic Mean Response Functions")

## ------------------------------------------------------------------------
plot(y ~ x, .data, pch = 19, col = "gray40", xlab = "Experience (X)", ylab = "Fitted Value")
title("Logistic (Gray), Probit (Blue), and \nComplementary Log-Log (Red) Fitted Models")

fit <- glm(y ~ x, data = .data, family = binomial(link = "probit"))
lines(xx, predict(fit, data.frame(x = xx), type = "resp"), col = 'blue', lwd = 2, lty = 2)

fit <- glm(y ~ x, data = .data, family = binomial(link = "cloglog"))
lines(xx, predict(fit, data.frame(x = xx), type = "resp"), col = 'red', lwd = 2, lty = 2)

fit <- glm(y ~ x, family = binomial(link = "logit"), .data)
lines(xx, predict(fit, data.frame(x = xx), type = "resp"), col = 'gray40', lwd = 2)

## ------------------------------------------------------------------------
.data <- CH14TA02
names(.data) <- c("x", "n", "y", "p")

## ------------------------------------------------------------------------
# Recreate the data from the summary table
fit <- glm(y ~ x, family = binomial(link = 'logit'), 
           data = data.frame(y = c(rep(1, sum(.data$y)), rep(0, 1000 - sum(.data$y))),
                             x = c(rep(.data$x, .data$y), rep(.data$x, 200 - .data$y))))
exp(coef(fit)) 

plot(p ~ x, .data, pch = 19, col = "gray40", ylim = c(0, 1), xlim = c(0, 40),
     xlab = "Price Reduction ($)", ylab = "Proportion Redeemed")
lines(seq(0, 40, len = 200), predict(fit, data.frame(x = seq(0, 40, len = 200)), type = "resp"))
title("Proportion of Coupons Redeemed and \nFitted Logistic Response Function")

## ------------------------------------------------------------------------
.data <- CH06FI05
names(.data) <- c("x1", "x2", "y")
.data <- transform(.data, y = ifelse(y < median(y), 0, 1))  # Create dichotomous response

## ------------------------------------------------------------------------
glm3d <- function(model, res, ...) {
  ivs <- labels(terms(model))
  x <- model.frame(model)[, ivs[1]]     # 1st independent variable
  y <- model.frame(model)[, ivs[2]]     # 2nd independent variable  
  xr <- seq(min(x), max(x), len = res)  # equidistant sequence from range of 1st iv
  yr <- seq(min(y), max(y), len = res)  # equidistant sequence from range of 2nd iv
  
  # Create a list object of these sequences to be turned into a grid
  newdata <- list(xr, yr)
  names(newdata) <- ivs
  newdata <- do.call('expand.grid', newdata)
  
  zr <- matrix(predict(model, newdata, type = 'resp'), res, res)
  persp(xr, yr, zr, ...)
}
fit <- glm(y ~ x1 + x2, binomial(link = 'logit'), .data)
glm3d(fit, res=45, tick = "detailed", shade = 0.5, expand = 2/3, theta = 300, phi = 30, 
      xlab = "x1", ylab = "x2", zlab = "y")
glm3d(fit, res=45, tick = "detailed", shade = 0.5, expand = 2/3, theta = 330, phi = 0, 
      xlab = "x1", ylab = "x2", zlab = "y")

## ------------------------------------------------------------------------
.data <- CH14TA03
names(.data) <- c('id', 'x1', 'x2', 'x3', 'x4', 'y')

## ------------------------------------------------------------------------
fit <- glm(y ~ x1 + x2 + x3 + x4, data = .data, family = binomial(link = 'logit'))
cbind(.data, 'fitted' = fitted(fit))

cbind(summary(fit)$coefficients, 'OR' = exp(coef(fit)))
vcov(fit)

## ------------------------------------------------------------------------
.data <- APPENC11
.data <- .data[2:3]
names(.data) <- c('Y', 'X')
.data <- transform(.data, X = log(X))  # log scale X
.data <- transform(.data, x = scale(X, T, F), xx = scale(X, T, F)^2)  # Center log X and square it

## ----poly-plot, fig.width=12, fig.height=6-------------------------------
par(mfrow = c(1,2))
xr <- with(.data, seq(min(X), max(X), len = 200))  # Plotting range values from log of X

fit <- glm(Y ~ X, binomial(link = 'logit'), .data)
plot(Y ~ X, .data, pch = 19, col = "gray")
title("(a) First-Order Fit")
lines(xr, predict(loess(Y ~ X, .data), data.frame(X = xr)), lty = 2, lwd = 2, col = 'blue')
lines(xr, predict(fit, data.frame(X = xr), type = "resp"), lwd = 2)

fit <- glm(Y ~ x + xx, binomial(link = "logit"), .data)
plot(Y ~ X, .data, pch = 19, col = "gray")
title("(b) Second-Order Fit")
lines(xr, predict(loess(Y ~ X, .data), data.frame(X = xr)), lty = 2, lwd = 2, col = 'blue')
lines(xr, predict(fit, data.frame(x = scale(xr, T, F), xx = scale(xr, T, F)^2), type = "resp"), lwd = 2)

summary(fit)  # TABLE 14.5

## ------------------------------------------------------------------------
.data <- CH14TA01
names(.data) <- c("x", "y", "fitted")

## ------------------------------------------------------------------------
fit <- glm(y ~ x, data = .data, family = binomial)
summary(fit)  # the results are already supplied

## ----wald-test, eval=FALSE-----------------------------------------------
#  
#  cat('z* =', coef(fit)[2] / sqrt(diag(vcov(fit))[2]))
#  cat('z(0.95) =', qnorm(0.95))
#  cat('p-value =', 1-pnorm(coef(fit)[2] / sqrt(diag(vcov(fit))[2])))  # Input 'z*'

## ----wal-test-show, echo=FALSE-------------------------------------------
cat('z* =', coef(fit)[2] / sqrt(diag(vcov(fit))[2]))
cat('z(0.95) =', qnorm(0.95))
cat('p-value =', 1-pnorm(coef(fit)[2] / sqrt(diag(vcov(fit))[2])))  # Input 'z*'

## ----int-est, eval=FALSE-------------------------------------------------
#  confint.default(fit)[2, ]                       # this result is based on asymptotic normality
#  exp(5 * coef(fit)[2])                           # Point estimate for 5 months of training
#  exp(5 * confint.default(fit)[2, ])              # Confidence limits for 5 months of training
#  (exp(5 * confint.default(fit)[2, ]) - 1) * 100  # Percentage change form

## ----int-est-show, echo=FALSE--------------------------------------------
confint.default(fit)[2, ]
exp(5 * coef(fit)[2])
exp(5 * confint.default(fit)[2, ])
(exp(5 * confint.default(fit)[2, ]) - 1) * 100

## ----boot, message=FALSE, warning=FALSE----------------------------------
boot.beta <- function(data, indices) {coef(glm(y ~ x, binomial, data = data[indices, ]))[2]}
z <- boot(.data, boot.beta, R = 1000)
boot.ci(z, type = c("perc", "bca"))

## ------------------------------------------------------------------------
.data <- CH14TA03
names(.data) <- c('id', 'x1', 'x2', 'x3', 'x4', 'y')

