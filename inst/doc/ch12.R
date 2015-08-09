## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(lmtest)

data("CH12TA02", package = "ALSM")

## ------------------------------------------------------------------------
error <- function(u, estart = 3.0, plot = FALSE, ...)
{
  e  <- vector(mode = "numeric", length = length(u))  
  y  <- vector(mode = "numeric", length = length(u))
  Xt <- 0:10

  for(i in seq(u)) 
    if (i != 1) e[i] <- e[i-1] + u[i] else e[i] <- estart
  
  y <- 2 + 0.5 * Xt + e  # Define the response from systematic Xt and error term u

  if (plot) 
    plot(e, pch = 19, ...)

  return(data.frame(u = u, e = e, y = y, x = 0:10))
}
 
u <- c(0, 0.5, -0.7, 0.3, 0, -2.3, -1.9, 0.2, -0.3, 0.2, -0.1)
.data <- error(u, xlab = "Time", plot = TRUE, ylim = c(-3, 4))
abline(0, 0)
print(.data)

## ------------------------------------------------------------------------
plot(y ~ x, .data, ylim = c(0, 10), main = "(a) True Regression Line e0 = 3", pch = 19)
abline(2, 0.5)

plot(y ~ x, .data, ylim = c(0, 10), main = "(b) Fitted Regression Line e0 = 3", pch = 19)
abline(lm(y ~ x, .data))

# Create new data set
.data <- error(rnorm(11), estart = -0.2)
plot(y ~ x, .data, ylim = c(0, 10), main = "(c) Fitted Regression Line e0 = -0.2", pch = 19)
abline(lm(y ~ x, .data))

## ------------------------------------------------------------------------
.data <- CH12TA02
names(.data) <- c("Y", "X")
fit  <- lm(Y ~ X, .data)
t    <- seq(2, 20)

## ------------------------------------------------------------------------
plot(resid(fit), xlab = "Time", ylim = c(-0.5, .5), pch = 19)
abline(0,0)

## ----dwtest, message=FALSE-----------------------------------------------
tab <- as.table(cbind(
  '(1)' = .data$Y,
  '(2)' = .data$X,
  '(3)' = resid(fit),
  '(4)' = c(NA,  resid(fit)[t] - resid(fit)[t-1]),
  '(5)' = c(NA, (resid(fit)[t] - resid(fit)[t-1])^2),
  '(6)' = resid(fit)^2))

round(tab, 4)

tab <- (resid(fit)[t] - resid(fit)[t-1])^2 / 
        sum(resid(fit)^2)
sum(tab)  # D = 0.735 (p 488)

dwtest(fit)  # alternative hypotheses defaults to rho > 0

fit.ordinary <- fit  # Save for end of chapter

## ------------------------------------------------------------------------
# An NA value is prefixed to each vector to fill out the t+1 = n observations
as.table(cbind(
  '(1)' = resid(fit),
  '(2)' = c(NA, resid(fit)[t-1]),
  '(3)' = c(NA, resid(fit)[t-1] * resid(fit)[t]),
  '(4)' = c(NA, resid(fit)[t-1]^2))) 

c("rho" = sum(resid(fit)[t-1] * resid(fit)[t]) / sum(resid(fit)[t-1]^2))

## ----cochrane-orcutt, message=FALSE--------------------------------------
cochrane.orcutt <- function(model) 
{
  x      <- model.matrix(model)[, -1]
  y      <- model.response(model.frame(model))
  e      <- resid(model)
  n      <- length(e)
  t      <- 2:n
  r      <- sum(e[t-1] * e[t]) / sum(e[t-1]^2)     # (12.22)    Step 1
  y      <- y[t] - r * y[t-1]                      # (12.18a)
  x      <- x[t] - r * x[t-1]                      # (12.18b)
  model  <- lm(y ~ x)                              # (12.19)    Step 2
  model$rho <- r
  return(model)
}

fit <- lm(Y ~ X, .data)
fit <- cochrane.orcutt(fit)
cbind(.data, 
      rbind(NA, model.frame(fit)))

# Generated in function above and stored with model to be used below.
c('rho' = fit$rho)  
summary(fit)
anova(fit)

# Durbin-Watson test on Cochrane-Orcutt model after 1 iteration
dwtest(fit)  # DW = 1.6502, fail to reject H0

# (12.23) Convert transformed coefficients back to original
cat("Y = ", coef(fit)[1] / (1 - fit$rho), " + ", coef(fit)[2], "X\n", sep = "")

fit.cochrane <- fit  # Save for end of chapter

## ------------------------------------------------------------------------
rho <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98)

hildreth.lu <- function(rho, model)
{
  x <- model.matrix(model)[, -1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[t] - rho * y[t-1]
  x <- x[t] - rho * x[t-1]
  
  return(lm(y ~ x))
}

fit <- lm(Y ~ X, .data)

# Minimum SSE = 0.07167118 for rho = 0.96. See plot below.
tab <- data.frame('rho' = rho, 
             'SSE' = sapply(rho, function(r) {deviance(hildreth.lu(r, fit))}))
round(tab, 4)

plot(SSE ~ rho, tab, type = 'l')
abline(v = tab[tab$SSE == min(tab$SSE), 'rho'], lty = 3)

fit <- hildreth.lu(0.96, fit)
summary(fit)
dwtest(fit)

# (12.28) Convert transformed coefficients back to original
cat("Y = ", coef(fit)[1] / (1 - 0.96), " + ",  coef(fit)[2], "X", sep = "")

fit.hildreth <- fit  # Save for end of chapter

## ------------------------------------------------------------------------
# One could also just do lm(diff(Y) ~ 0 + diff(X), .data), but the coefficients would be oddly named diff(X).
fit <- lm(y ~ x - 1, data.frame(y = diff(.data$Y), x = diff(.data$X)))

cbind(
  '(1)' = .data$Y, 
  '(2)' = .data$X, 
  '(3)' = c(NA, diff(.data$Y)),
  '(4)' = c(NA, diff(.data$X))) 

summary(fit)
anova(fit)
dwtest(fit)

# (12.33) Convert transformed coefficient back to original
cat("Y = ", mean(.data$Y) - coef(fit)*mean(.data$X), " + ", coef(fit), "X\n", sep = "")

fit.difference <- fit  # Save for end of chapter

## ------------------------------------------------------------------------
# First column of slope coefficients
c('Cochrane'   = c(coef(fit.cochrane)[2]),
  'Hildreth'   = c(coef(fit.hildreth)[2]), 
  'Difference' = c(coef(fit.difference)), 
  'Ordinary'   = c(coef(fit.ordinary)[2]))

# Second column of sloope coefficient standard errors
c('Cochrane'   = round(summary(fit.cochrane)[['coefficients']][2, 2], 4),
  'Hildreth'   = round(summary(fit.hildreth)[['coefficients']][2, 2], 4),
  'Difference' = round(summary(fit.difference)[['coefficients']][1, 2], 4),
  'Ordinary'   = round(summary(fit.ordinary)[['coefficients']][2, 2], 4))

