## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(boot)  # boot
library(HH)    # hov
library(nnet)  # nnet

data("CH13TA01", package = "ALSM")
data("CH13TA04", package = "ALSM")
data("APPENC09", package = "ALSM")


## ----example-plots, fig.width=12, fig.height=6---------------------------
par(mfrow = c(1, 2))
curve(100 - 50*exp(-2*x), xlim = c(0, 4), ylim = c(50, 100), xlab = "", ylab = "")
title("(a) Exponential Model (13.8)")

curve(10/(1 + 20*exp(-2*x)), xlim = c(0, 4), ylim = c(0, 10), xlab = "", ylab = "")
title("(b) Logistic Model (13.10)")

## ------------------------------------------------------------------------
.data <- CH13TA01
names(.data) <- c("index", "days")
fit <- nls(index ~ a * exp(b * days), .data, start = c(a = 56.6646, b = -0.03797))

## ------------------------------------------------------------------------
cbind('Patient'           = seq(15),
      'Days Hospitalized' = .data$days,
      'Prognostic Index'  = .data$index)

## ------------------------------------------------------------------------
plot(index ~ days, .data, pch = 19)
curve(coef(fit)[1] * exp(coef(fit)[2] * x), add = TRUE)

## ------------------------------------------------------------------------
nlm <- function(formula, data, g, FUN = function(x, a, b) a * exp(b * x), TOL = 1e-5)
{
  YX  <- model.frame(formula, data)
  x   <- YX[, -1]
  y   <- YX[, 1]
  n   <- length(y)
  SSE <- NULL
  i   <- 1
  G   <- matrix(nrow = 50, ncol = 2)

  repeat
  {
      G[i, ] = g
      Y   = as.matrix(y - FUN(x, g[1], g[2]))               # (13.25a)
      D   = cbind(exp(g[2] * x), g[1] * x * exp(g[2] * x))  # (13.25b)
      b   = coef(lm(Y ~ D - 1))                             # (13.25c)
      SSE = c(SSE, sum(Y^2))
      g <- g + b
      is.done = ((length(SSE) > 1) && (SSE[i-1] - SSE[i] < TOL))
      i = i + 1
      if (is.done) break
  }

  G <- G[!is.na(G[,1]), ]
  i = i-1
  
  tab <- list()
  tab$Estimates = cbind('g' = G, 'SSE' = SSE)
  tab$LS        = cbind(G[i, ], rbind(SSE[i-1] / (n-2), SSE[i-1]))
  tab$se        = (SSE[i] / (n-2)) * solve(t(D) %*% D)
  tab$fitted    = G[i, 1] * exp(G[i, 2] * x)
  tab$Y         = Y
  tab$D         = D
  
  return(tab)
}

(coefs <- coef(lm(log(index) ~ days, .data)))  # Linear Transform Fit
(coefs <- c(exp(coefs[1]), coefs[2]))       # Exp. Transformation
nlm(index ~ days, .data, g = coefs)

## ----diag-plots, message=FALSE-------------------------------------------
plot(resid(fit) ~ fitted(fit), xlab = "Fitted Values", ylab = "Residual")
title("(a) Residual Plot against Fitted")
abline(0, 0)

qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
qqline(resid(fit))
title("(b) Normal Probability Plot")
abline(0, 0)

hov(resid(fit) ~ cut(.data$days, 2))  # BF test. P-value = 0.64
deviance(fit) / df.residual(fit)   # (13.31)
vcov(fit)                          # (13.32b)

## ----boot, message=FALSE, warning=FALSE----------------------------------
boot.coef <- function(data, indices, g)
{
  x   <- data[indices, 2]
  y   <- data[indices, 1]
  mod <- nls(y ~ a * exp(b*x), .data, start = c(a = g[[1]], b = g[[2]]))
  
  return(coef(mod))
}

# Warnings are produced
(z <- boot(.data, boot.coef, R = 1000, g = coef(lm(log(index) ~ days, .data))))

# NLS coefficients and Bootstrap CIs
coef(fit)
boot.ci(z, index = 1)
boot.ci(z, index = 2)


g0 = expression(g[0]^symbol("\052"))
g1 = expression(g[1]^symbol("\052"))
hist(z$t[, 1], 40, freq = FALSE, xlab = g0, main = "")
title(expression("(a) Histogram of Bootstrap Estimates " * g[0]^symbol("\052")))

hist(z$t[, 2], 40, freq = FALSE, xlab = g1, main = "")
title(expression("(b) Histogram of Bootstrap Estimates " * g[1]^symbol("\052")))

# Simultaneous Interval Estimation (p 532)
s = sqrt(diag(vcov(fit)))
(coef(fit) + 2.16 * c(lwr = -s, upr = s))[c(1,3, 2,4)]

## ------------------------------------------------------------------------
.data <- CH13TA04
names(.data) <- c("location", "week", "efficiency")

## ------------------------------------------------------------------------
xtabs(efficiency ~ week + location, .data)  

## ------------------------------------------------------------------------
boot.coef <- function(data, indices, g) 
{
  y   <- data[indices, 1]
  x1  <- data[indices, 2]
  x2  <- data[indices, 3]
  mod <- nls(y ~ a  + b*x1 + c*exp(d*x2), start = g)
  
  return(coef(mod))
}

g   <- c(a = 1.025, b = -0.0459, c = -0.5, d = -0.122)
fit <- nls(efficiency ~ a + b*location + c*exp(d*week), .data, start = g)
z   <- boot(.data[c(3, 1:2)], boot.coef, R = 1000, g = g)

cbind(
  'g'   = g,
  '(1)' = summary(fit)$coef[, 1],
  '(2)' = summary(fit)$coef[, 2],
  '(3)' = apply(z$t, 2, mean),
  '(4)' = apply(z$t, 2, sd)) 

# Residual Plots That Were Not Shown
par(mfrow = c(2, 2))
plot(resid(fit) ~ fitted(fit))
title("Residuals aganist Fitted")

plot(resid(fit) ~ location, .data)
title("Residuals against Location")

plot(resid(fit) ~ week, .data)
title("Residuals against Week")

qqnorm(resid(fit), main = "")
qqline(resid(fit))
title("Normal Probability Plot")

## ------------------------------------------------------------------------
p <- coef(fit)
plot(efficiency ~ week, .data, col = factor(location), pch = 19)

curve(p[1] + p[2]*(0) + p[3] * exp(p[4] * x), to = 90, add = T)  # Location = 0
curve(p[1] + p[2]*(1) + p[3] * exp(p[4] * x), to = 90, add = T)  # Location = 1

## ----boot-hist, fig.width=12, fig.height=6-------------------------------
par(mfcol = c(2,4))
for(i in seq(4)) {
  hist(z$t[, i], 50, freq = FALSE, main = "", xlab = paste("g", i, sep = ""))
  title(paste("(", letters[i], ")", sep = ""))
  
  plot(density(z$t[, i]), xlab = paste("g", i, sep = ""), main = "")
  title(paste("(", letters[i], ")", sep = ""))
}  

## ----log-plot, fig.width=12, fig.height=4--------------------------------
f <- function(x, a, b) {(1 + exp(-a - b*x))^(-1)}

par(mfrow = c(1, 3), lwd = 2)
curve(f(x, a =  0, b = 0.1), -10, 10, ylim = c(0, 1))
curve(f(x, a =  0, b =   1), -10, 10, add = TRUE, lty = 5)
curve(f(x, a =  0, b =  10), -10, 10, add = TRUE, lty = 3)
title("(a)")

curve(f(x, a =  0, b = -0.1), -10, 10, ylim = c(0, 1))
curve(f(x, a =  0, b =   -1), -10, 10, add = TRUE, lty = 5)
curve(f(x, a =  0, b =  -10), -10, 10, add = TRUE, lty = 3)
title("(b)")

curve(f(x, a =  5, b = 1), -10, 10, ylim = c(0, 1))
curve(f(x, a =  0, b = 1), -10, 10, add = TRUE, lty = 5)
curve(f(x, a = -5, b = 1), -10, 10, add = TRUE, lty = 3)
title("(c)")

## ------------------------------------------------------------------------
.data <- APPENC09
names(.data) <- c("id", "cost", "age", "gender", "intervention", "drugs", 
               "visits", "complications", "comorbidities", "duration")

## ---- results=FALSE, eval=FALSE------------------------------------------
#  set.seed(666)
#  indices  <- sample(seq(nrow(.data)), 400)
#  .data.train <- .data[indices, ]
#  .data.test  <- .data[-indices, ]
#  fit <- nnet(log(cost) ~ intervention + drugs + comorbidities + complications, .data.train, size = 5, maxit = 50)

