## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(MASS)     # rlm
library(car)      # spm, vif
library(bestglm)  # bestglm, regsubsets
library(leaps)    # regsubsets
library(rpart)    # rpart
library(boot)     # boot

data("CH07TA01", package = "ALSM")
data("CH11TA01", package = "ALSM")
data("CH11TA04", package = "ALSM")
data("CH11TA07", package = "ALSM")
data("CH11TA08", package = "ALSM")
data("CH11TA09", package = "ALSM")
data("CH11TA10", package = "ALSM")
data("CH11TA11", package = "ALSM")
data("APPENC04", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH11TA01
names(.data) <- c("X", "Y")

## ------------------------------------------------------------------------
fit <- lm(Y ~ X, .data)
w <- fitted(lm(abs(resid(fit)) ~ X, .data))^(-2)
fit <- update(fit, weights = w)
 
cbind(
  "Subject" = seq(nrow(.data)),
  "Age"     = .data$X,
  "DBP"     = .data$Y,
  "e"       = round(resid(lm(Y ~ X, .data)), 2),
  "|e|"     = round(abs(resid(lm(Y ~ X, .data))), 2),
  "s"       = round(sqrt(1 / w), 4),
  "w"       = round(fit$weights, 5))

## ----diag-plot, fig.height=3---------------------------------------------
fit <- lm(Y ~ X, .data)
par(mfrow = c(1, 3), pch = 19)

plot(Y ~ X, .data, xlab = "Age", ylab = "Blood Pressure")
title("(a) Scatter Plot")
abline(fit)

plot(resid(fit) ~ X, .data, xlab = "Age", ylab = "Residual")
title("(b) Residual Plot Against X")
abline(0, 0)

plot(abs(resid(fit)) ~ X, .data, xlab = "Age", ylab = "Absolute Residual")
title("(c) Absolute Residual Plot Against X")
abline(lm(abs(resid(fit)) ~ X, .data))

## ------------------------------------------------------------------------
.data <- CH07TA01
names(.data) <- c("x1", "x2", "x3", "y")
fit <- lm(y ~ x1 + x2 + x3, .data)

## ------------------------------------------------------------------------
c <- c(0, .002, .004, .006, .008, .01, .02, .03, .04, .05, .1, .5, 1)
fun1 <- function(c, rxx, rxy) {solve(rxx + c * diag(ncol(rxx)), rxy)}  # (11.34)
fun2 <- function(c, rxx) 
{
  bias = solve(rxx + c*diag(ncol(rxx)))
  diag(bias %*% rxx %*% bias)                                          # (11.36)
}

data.frame(lambda = c, 
           t(sapply(c, fun1, rxx = cor(.data)[-4, -4], rxy = cor(.data)[-4, 4])),
           VIF = t(sapply(c, fun2, rxx = cor(.data)[-4, -4])))

## ------------------------------------------------------------------------
lambda <- seq(0.001, 1, by = 0.001)
f <- function(c, x) {solve(cor(x)[-4, -4] + c * diag(3), cor(x)[-4, 4])}
tab <- data.frame(c = lambda, 
                  t(sapply(lambda, f, x = .data)))

par(lwd = 2)
plot(x1 ~ c, tab, type = "l", log = "x", xaxt = "n", lty = 5,
     xlim = c(0.001, 1), ylim = c(-1, 3), xlab = "", ylab = "")
axis(1, c(0.001, 0.01, .10, 1.00))
lines(x2 ~ c, tab, lty = 3)
lines(x3 ~ c, tab, lty = 1)
abline(0,0)
abline(v = 0.02, lwd = 1)  # Authors choice of c
text(0.02, -1, pos=4, labels = "0.02")

## ------------------------------------------------------------------------
.data <- CH11TA04
names(.data) <- c("State", "Y", "X1", "X2", "X3", "X4", "X5")

## ----many-fits, fig.width=10, fig.height=14------------------------------
library(MASS)
par(mfrow = c(3, 2), pch = 19)
xlab <- expression(X[2])

# Linear Regression Fit
fit  <- lm(Y ~ X2, .data)
outliers <- subset(.data, State %in% c("Guam", "D.C.", "Virgin_Islands"))

plot(Y ~ X2, .data, xlab = xlab, ylab = "Y")                          # (a)
lines(loess.smooth(.data$X2, .data$Y), lty = 2)
title("(a) Lowess and Linear Regression Fits")
abline(fit)
with(outliers, text(Y ~ X2, labels = State, pos = 4))

plot(resid(fit) ~ X2, .data, xlab = xlab, ylab = "Residual")          # (b)
title("(b) Residuals from Linear Regression")
abline(0, 0)

# OLS Quadratic Fit
fit <- lm(Y ~ x + xsq, transform(.data,
  x   = scale(X2, T, F),
  xsq = scale(X2, T, F)^2))
x <- data.frame(X2 = .data$X2, Y = fitted(fit))

plot(Y ~ X2, .data, xlab = expression(X[2]), ylab = "Y")              # (c)
title("(c) OLS Quadratic Fit")
lines(Y ~ X2, x[order(x$X2), ])

plot(cooks.distance(fit), xlab = "Index", ylab = "D", type = "o")  # (d)
title("(d) Cook's Distance--OLS Quadratic Fit")

# Robust Quadratic Fit
fit <- rlm(formula(fit), model.frame(fit))
x <- data.frame(X2 = .data$X2, Y = fitted(fit))

plot(Y ~ X2, .data, xlab = expression(X[2]), ylab = "Y")              # (e)
title("(e) Robust Quadratic Fit")
lines(Y ~ X2, x[order(x$X2), ])

plot(fit$w, xlab = "Index", ylab = "W", type = "o")                # (f)
title("(f) Robust Weights")

## ------------------------------------------------------------------------
with(.data, cbind("State"    = State,
               "MATHPROF" = Y,
               "PARENTS"  = X1,
               "HOMELIB"  = X2,
               "READING"  = X3,
               "TVWATCH"  = X4,
               "ABSENCES" = X5)) 

## ------------------------------------------------------------------------
fit <- lm(Y ~ x + xsq, transform(.data,
                                 x   = scale(X2, T, F),
                                 xsq = scale(X2, T, F)^2))

IRLS <- function(model, i = 1) 
{
  e   <- resid(model)
  MAD <- median(abs(e - median(e))) / 0.6745  # (11.46)
  u   <- e / MAD  # (11.47)
  w   <- apply(data.frame(u), 1, function(x) if (abs(x) <= 1.345) 1 else 1.345/abs(x))  # (11.44)

  model <- update(model, weights = w)
  if (i > 1) return(IRLS(model, i-1)) else return(model)  # Recursive return definition
}

tab <- cbind(
  "e0" = resid(fit),
  "u0" = resid(fit) / (median(abs(resid(fit) - median(resid(fit)))) / 0.6745),
  "w1" = IRLS(fit, 1)$weights,
  "e1" = resid( IRLS(fit, 1) ),
  "w2" = IRLS(fit, 2)$weights,
  "e2" = resid( IRLS(fit, 2) ),
  "w7" = IRLS(fit, 7)$weights,
  "e7" = resid( IRLS(fit, 7) ))

round(tab, 4)

## ----spm-plot, message=FALSE, warning=FALSE, fig.width=12, fig.height=14----
library(car)
cor(.data[-1])
spm(.data[-1], span = 0.75, diagonal = "none", reg.line = FALSE, spread = FALSE,
  var.labels = c("MATHPROF", "PARENTS", "HOMELIB", 
                 "READING", "TVWATCH", "ABSENCES"))

## ------------------------------------------------------------------------
fit <- lm(Y ~ X1 + X2 + X3 + X4 + X5, .data)
data.frame("State" = levels(.data$State),
           "h"     = round(hatvalues(fit), 2),
           "t"     = round(rstudent(fit), 2),
           "D"     = round(cooks.distance(fit), 2))

## ----not-presented, fig.width=7, fig.height=10---------------------------
par(mfrow = c(3, 2))
for (n in seq(5))
{
    x <- model.frame(fit)[, n+1]
    predictor <- paste("X", n, sep = "")
    plot(resid(fit) ~ x, ylab = "Residual", main = predictor)
    abline(0, 0)
}

## ----best-subsets, message=FALSE-----------------------------------------
library(bestglm)
bestglm(model.frame(fit)[c(2:6, 1)], IC = "AIC")$Subsets  # Author's choice is 3 on each search
bestglm(model.frame(fit)[c(2:6, 1)], IC = "BIC")$Subsets

summary(rlm(Y ~ X2 + X3 + X4, .data))  # (11.53) 
summary(lm (Y ~ X2 + X3 + X4, .data))  # (11.54)

with(summary(regsubsets(formula(fit), .data, nbest=2)),
     cbind(which[, -1], rss, rsq, adjr2, cp, bic))

## ------------------------------------------------------------------------
.data <- CH11TA07
names(.data) <- c("X1", "X2", "Y")

## ------------------------------------------------------------------------
fit <- lm(Y ~ X1 + X2, .data)  # (6.1)
lo <- function(model, newdata = data.frame(X1 = 30, X2 = 3), q = 0.5) 
{
  x<- transform(model.frame(model),
    d = sqrt(scale(X1, newdata$X1, sd(X1))^2 + scale(X2, newdata$X2, sd(X2))^2))
  dq <- with(x, sort(d)[nrow(x) * q])
  x <- transform(x, w = ifelse(d > dq, 0, (1-(d/dq)^3)^3))
  with(x, update(model, weights = w))
}

transform(.data, 
  d = sqrt(scale(X1, 30, sd(X1))^2 + scale(X2, 3, sd(X2))^2),
  w = lo(fit)$weights)

predict(lo(fit), data.frame(X1 = 30, X2 = 3))  # 4.664607

# Compare
fit <- loess(Y ~ X1 + X2, .data, span = 0.5, degree = 1)
predict(fit, data.frame(X1 = 30, X2 = 3))  # returns 4.2894 with normalize = FALSE

## ------------------------------------------------------------------------
lm3d <- function(model, res, ...) {
  ivs <- labels(terms(model))
  x <- model.frame(model)[, ivs[1]]     # 1st independent variable
  y <- model.frame(model)[, ivs[2]]     # 2nd independent variable  
  xr <- seq(min(x), max(x), len = res)  # equidistant sequence from range of 1st iv
  yr <- seq(min(y), max(y), len = res)  # equidistant sequence from range of 2nd iv

  # Create a list object of these sequences to be turned into a grid
  newdata <- list(xr, yr)
  names(newdata) <- ivs
  newdata <- do.call('expand.grid', newdata)

  zr <- matrix(predict(model, newdata), res, res)
  persp(xr, yr, zr, ...)
  
  invisible(list(x = xr, y = yr, z = zr))  # included in this example for `contour`
}

fit <- loess(Y ~ X1 + X2, .data, span = 0.5, degree = 1)
lm3d(fit, res = 30, theta = 270, phi = 0, ticktype = "detailed", expand = 2/3, shade = 0.5)

# Store the dimensional sequences for the plots below
.data <- lm3d(fit, res = 30, theta = 300, phi = 30, ticktype = "detailed", expand = 2/3, shade = 0.5)
with(.data, contour(x, y, z))

## ----coplot, fig.width=12, fig.height=4----------------------------------
par(mfrow = c(1, 3))
plot(.data$x, predict(fit, data.frame(X1=.data$x, X2=3)), type='l', ylim=c(0, 360), ylab = 'y', xlab = 'X1')
title("X2 = 3")

plot(.data$x, predict(fit, data.frame(X1=.data$x, X2=6)), type='l', ylim=c(0, 360), ylab = 'y', xlab = 'X1')
title("X2 = 6")

# Very different scale than book
plot(.data$x, predict(fit, data.frame(X1=.data$x, X2=9)), type='l', ylim=c(-100, 360), ylab = 'y', xlab = 'X1')
title("X2 = 9")

## ------------------------------------------------------------------------
.data <- CH11TA08
names(.data) <-c("Y", "X")

## ------------------------------------------------------------------------
.data <- transform(.data, 
                R = cut(X, c(8, 9, 10, 13, 14, 25), include.lowest = TRUE, right = FALSE))
with(.data, tapply(Y, R, mean))

## ------------------------------------------------------------------------
library(rpart)

fit <- rpart(Y ~ X, .data, control = rpart.control(minsplit = 2, cp = 0.013))

# FIGURE 11.9a
plot(Y ~ X, .data, xlab = "Age", ylab = "Steroid Level")
lines(predict(fit, data.frame(X = sort(X))) ~ sort(X), .data, type = 's')  # Step plot

# FIGURE 11.9b 
plot(resid(fit) ~ predict(fit), xlab = "Predicted", ylab = "Residual")

# FIGURE 11.9c
plot(fit)
text(fit, use.n = TRUE, xpd = TRUE)

## ------------------------------------------------------------------------
rplot <- function(formula, data, ...)
{
  .data <- model.frame(formula, data)
  model <- rpart(Y ~ X, .data, ...)
  plot(formula, .data, xlab = "AGE", ylab = "Steroid Level")
  lines(predict(model, data.frame(X = sort(X))) ~ sort(X), .data, type = 's')
}

par(mfrow = c(2, 2))
rplot(Y ~ X, .data, minsplit = 2, cp = 0.095)
rplot(Y ~ X, .data, minsplit = 2, cp = 0.04)
rplot(Y ~ X, .data, minsplit = 2, cp = 0.03)
rplot(Y ~ X, .data, minsplit = 2, cp = 0.015)

## ------------------------------------------------------------------------
.data <- APPENC04
names(.data) <- c("id", "gpa", "rank", "score", "year")

## ------------------------------------------------------------------------
library(rpart)
set.seed(38383)  # I want more control over what model gets generated here.
indices <- sample(nrow(.data), 400)
.data.test <- .data[indices, ]
.data.valid <- .data[-indices, ]
fit <- rpart(gpa ~ score + rank, .data.test, cp = 0.015, model = TRUE)

mspr <- function(model, validation) 
{
  x <- model.frame(formula(model), validation)
  x <- sum((x[, 1] - predict(model, x[-1]))^2) / nrow(x)  # (9.20)
  
  return(x)
}

x <- NULL
for (cp in seq(0.01, 0.04, by = 0.001))
  x <- append(x, mspr(rpart(gpa ~ rank + score, .data.test, cp = cp), .data.valid))

x <- NULL
for (n in seq(30)) {
  for (cp in seq(0.01, 0.03, by = 0.001)) {
    x <- append(x, mspr(rpart(gpa ~ rank + score, .data.test, maxdepth = n, cp = cp), .data.valid))
  }
}
x <- data.frame(depth = rep(seq(30), each = 21), cp = seq(0.01, 0.03, by = 0.001), mspr = x)

# FIGURE 11.12 (a)



# FIGURE 11.12 (b)
lm3d(fit, 30, theta = 300, phi = 30, ticktype = "detailed", expand = 2/3, shade = 0.5, xlab = "Score", ylab = "Rank", zlab = "GPA")
lm3d(fit, 30, theta = 35, phi = 30, ticktype = "detailed", expand = 2/3, shade = 0.5, xlab = "Score", ylab = "Rank", zlab = "GPA")



# FIGURE 11.12 (c)
plot(fit, compress = TRUE, uniform = TRUE)
text(fit, use.n = TRUE, xpd = TRUE, digits = 4)
title("(c)")

# FIGURE 11.12 (d)
plot(predict(fit), resid(fit), xlab = "Predicted GPA", ylab = "Residual")
title("(d)")

## ------------------------------------------------------------------------
.data <- CH11TA09
names(.data) <- c("X", "Y")
fit <- lm(Y ~ X, .data)

## ------------------------------------------------------------------------
estar <- sample(resid(fit), size = nrow(.data), replace = TRUE)
cbind(.data,
      'Yh'    = fitted(fit),
      'ei'    = resid(fit),
      'estar' = estar,
      'ystar' = fitted(fit) + estar)

## ------------------------------------------------------------------------
bootstrap <- function(model, times = 1000, alpha = 0.05) 
{
  b     <- coef(model)[[2]]
  n     <- nrow(model.frame(model))
  coefs <- vector(mode = "numeric", length = times)
  for (i in seq(times))
  {
    estar <- sample(resid(model), size = n, replace = TRUE)
    fstar <- fitted(model) + estar
    bootmodel  <- lm(fstar ~ ., data = model.frame(model))
    coefs[i] <- coef(bootmodel)[[2]]
  }

  p <- quantile(coefs, probs = c(alpha/2, 1-alpha/2))

  statistics <- cbind(
    "mean"      = mean(coefs),
    "sd"        = sd(coefs),
    "b(a/2)"    = p[[1]],
    "b(1-a/2)"  = p[[2]])

  confint <- cbind(
    'd1'  = b - p[[1]],
    'd2'  = p[[2]] - b,
    'lwr' = 2*b - p[[2]],
    'upr' = 2*b - p[[1]])

  return (list(coefs = coefs, statistics = statistics, confint = confint))
}

z <- bootstrap(fit)
z$statistics
z$confint
hist(z$coefs, breaks = 40, freq = F, main = "", xlab = "Bootstrap b1*")

## ------------------------------------------------------------------------
.data <- CH11TA10
names(.data) <- c("X", "Y")
fit  <- lm(Y ~ X, .data)

## ------------------------------------------------------------------------
xstar <- sample(.data$X, size = nrow(.data), replace = TRUE)
ystar <- sample(.data$Y, size = nrow(.data), replace = TRUE)
estar <- abs(resid(lm(ystar ~ xstar) ))
cbind(.data,
  'X*' = xstar,
  'Y*' = ystar,
  'e*' = round( resid( lm(ystar ~ xstar) ), 2),
  's*' = round( fitted( lm(estar ~ xstar) ), 2),
  'w*' = round( 1 / (fitted( lm(estar ~ xstar) ))^2, 4)) 

## ----bootstrap, message=FALSE--------------------------------------------
library(boot)
bootstrap <- function(model, times = 1000, alpha = 0.05)
{
  b     <- coef(model)[[2]]
  n     <- nrow(model.frame(model))
  coefs <- vector(mode = "numeric", length = times)
  for(i in seq(times)) 
  {
    indices  <- sample(1:n, size = n, replace = TRUE)
    xstar    <- model.frame(model)[indices, 2]
    ystar    <- model.frame(model)[indices, 1]
    mod      <- lm(ystar ~ xstar)
    mod      <- update(mod, weights = fitted(lm(abs(resid(mod)) ~ xstar))^(-2))
    coefs[i] <- coef(mod)[[2]]
  }

  p <- quantile(coefs, probs = c(alpha/2, 1-alpha/2))

  statistics <- cbind(
    "mean"     = mean(coefs),
    "sd"       = sd(coefs),
    "b(a/2)"   = p[[1]],
    "b(1-a/2)" = p[[2]])

  confint <- cbind(
    'd1'  = b - p[[1]],
    'd2'  = p[[2]] - b,
    'lwr' = 2*b - p[[2]],
    'upr' = 2*b - p[[1]])

  return (list(coefs = coefs, statistics = statistics, confint = confint))
}

boot.coef <- function(data, indices) 
{
  x   <- data[indices, 2]
  y   <- data[indices, 1]
  mod <- lm(y ~ x)
  mod <- update(mod, weights = fitted(lm(abs(resid(mod)) ~ x))^(-2))
  coef(mod)[[2]]
}

z <- bootstrap(fit)
z$statistics
z$confint
hist(z$coefs, breaks = 40, freq = F, main = "", xlab = "Bootstrap b1*")

# Using the boot package
z <- boot(data = model.frame(fit), statistic = boot.coef, R = 1000)
print(z)
boot.ci(z)
hist(z$t, breaks = 40, freq = F, main = "", xlab = "Bootstrap b1*")

## ------------------------------------------------------------------------
.data <- CH11TA11
names(.data) <- c("AADT", "CTYPOP", "LANES", "WIDTH", 
               "CONTROL", "CLASS", "TRUCK", "LOCALE")
.data <- transform(.data,
                CONTROL = factor(CONTROL),
                CLASS   = factor(CLASS),
                TRUCK   = factor(TRUCK),
                LOCALE  = factor(LOCALE))
fit <- lm(AADT ~ ., .data)

## ------------------------------------------------------------------------
{
  print(head(.data))
  cat("...\n")
  print(tail(.data))
}

## ----spm-mndot, warning=FALSE, message=FALSE, fig.width=12, fig.height=14----
library(car)
spm(.data, diagonal = "none", span = 0.75, reg.line = F, by.groups = TRUE)

## ----all-regs, message=FALSE---------------------------------------------
library(leaps)
library(car)

# Some Mentioned Diagnostics
plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual", pch = 19)
vif(fit)
max(cooks.distance(fit), na.rm = TRUE)

# Best Subsets Analysis
best <- regsubsets(formula(fit), .data, force.in = c(1, 2), nbest=5, nvmax = 7)
best$xnames <- paste("X", 0:13, sep = "")
with(summary(best), cbind(
    which[,-c(1:3)],
    Cp     = round(cp, 4),
    Rsq    = round(rsq, 4),
    AdjRsq = round(adjr2, 4))) 

## ----best-subs2----------------------------------------------------------
levels(.data$CLASS)[2] <- "4"
.data <- transform(.data, 
  CTYPOP2 = scale(CTYPOP, T, F)^2,
  LANES2  = scale(LANES, T, F)^2)

pool <- AADT ~ (CTYPOP + LANES + CONTROL + CLASS)^2 + CTYPOP2 + LANES2
best <- regsubsets(pool, .data, force.in = c(1, 2), nbest = 5)

best$xnames <- paste("X", 0:16, sep = "")
with(summary(best), cbind(
    which[,-c(1:3)],
    Cp     = round(cp, 4),
    Rsq    = round(rsq, 4),
    AdjRsq = round(adjr2, 4)))

## ----student-plot, fig.width=8, fig.height=4-----------------------------
par(mfrow = c(1, 2), pch = 19)

fit <- lm(AADT ~ CTYPOP + LANES + CONTROL + CLASS, .data)
plot(rstudent(fit) ~ fitted(fit), xlab = "Fitted Values", ylab = "Studentized Residual")
title("(a) First-Order OLS Model")

fit <- lm(AADT ~ CTYPOP + LANES + LANES2 + CONTROL + CTYPOP:CONTROL, .data)
plot(rstudent(fit) ~ fitted(fit), xlab = "Fitted Values", ylab = "Studentized Residual")
title("(b) Second-order OLS Model")

## ------------------------------------------------------------------------
fit <- lm(AADT ~ CTYPOP + LANES + LANES2 + CONTROL + CTYPOP:CONTROL, .data)
summary(fit)  # Compare before iterations
w   <- fitted(lm(abs(resid(fit)) ~ CTYPOP + LANES, .data))^(-2)
fit <- update(fit, weights = w)                           # 1st iteration
w   <- fitted(lm(abs(resid(fit)) ~ CTYPOP + LANES, .data))^(-2)
fit <- update(fit, weights = w)                           # 2nd iteration
w   <- fitted(lm(abs(resid(fit)) ~ CTYPOP + LANES, .data))^(-2)
fit <- update(fit, weights = w)                           # 3nd iteration
w   <- fitted(lm(abs(resid(fit)) ~ CTYPOP + LANES, .data))^(-2)
fit <- update(fit, weights = w)                           # 4nd iteration

summary(fit)  # Compare after iterations

## ----resid-plot, fig.width=8, fig.height=4-------------------------------
par(mfrow = c(1, 2), pch = 19)
plot(rstudent(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Student. Residual")
title("(a) Residual Plot against Fitted")

qqnorm(rstudent(fit), xlab = "Expected", ylab = "Student. Residual", main = "")
qqline(rstudent(fit))
title("(b) Normal Probability Plot")

