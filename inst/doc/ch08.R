## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(car)    # Used by alr3
library(alr3)   # pureErrorAnova
library(Rcmdr)  # scatter3d

data("CH08TA01", package = "ALSM")
data("CH07TA01", package = "ALSM")
data("CH08TA02", package = "ALSM")
data("CH08TA05", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH08TA01
names(.data) <- c("Y", "X1", "X2")
.data <- transform(.data, 
  x1   = scale(X1, scale = 0.4),
  x2   = scale(X2, scale = 10),
  x1sq = scale(X1, scale = 0.4)^2,
  x2sq = scale(X2, scale = 10)^2,
  x1x2 = scale(X1,, 0.4) * scale(X2,, 10))

fit <- lm(Y ~ ., .data)

## ------------------------------------------------------------------------
model.frame(.data)

## ------------------------------------------------------------------------
par(mfrow = c(2, 2), pch = 19)
plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual")
title("(a) Residual Plot Against Fitted Values")

plot(resid(fit) ~ x1, .data, xlab = "x1", ylab = "Residual")
title("(b) Residual Plot against x1")

plot(resid(fit) ~ x2, .data, xlab = "x2", ylab = "Residual")
title("(c) Residual Plot against x2")

qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
qqline(resid(fit))
title("(d) Normal Probability Plot")

## ----poly-out, message=FALSE---------------------------------------------
fit <- lm(Y ~ ., .data)                                 # Restate the model
summary(fit)
pureErrorAnova(fit)                                     # alr3 function
anova(lm(Y ~ 1, .data), fit)                            # Explicit F-Test


# Calculate F*     (p 302)
pureErrorAnova(fit)[" Lack of fit", "Mean Sq"] / pureErrorAnova(fit)[" Pure Error",  "Mean Sq"]  # (6.68b)

# Correlation (p 301)
with(.data, rbind(
  "X and Xsq" = round(c("1" = cor(X1, X1^2), "2" = cor(X2, X2^2)), 3),
  "x and xsq" = round(c("1" = cor(x1, x1^2), "2" = cor(x2, x2^2)), 3)))

anova(lm(Y ~ x1 + x2, .data), fit)                      # partial F-test
fit <- lm(Y ~ x1 + x2, .data)                           # refit model (p 304)

# Diagnostic Plots
par(mfrow = c(2, 2), pch = 19)                       # Same as FIGURE 8.5
plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual")
title("(a) Residual Plot Against Fitted Values")

plot(resid(fit) ~ x1, .data, xlab = "x1", ylab = "Residual")
title("(b) Residual Plot against x1")

plot(resid(fit) ~ x2, .data, xlab = "x2", ylab = "Residual")
title("(c) Residual Plot against x2")

qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
title("(d) Normal Probability Plot")
qqline(resid(fit))

confint(lm(Y ~ X1 + X2, .data))                         # Bonferroni procedure on original coefs. (p 305)

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
}

lm3d(lm(Y ~ X1 + X2, .data), res = 30, ticktype = "detailed", shade = 0.5, expand = 2/3,
     xlab = "Charge Rate", ylab = "Temp", zlab = "Cycles",
     theta = 320, phi = 30)

## ------------------------------------------------------------------------
.data <- CH07TA01
names(.data) <- c("X1", "X2", "X3", "Y")

# Add centered variables
.data <- transform(.data,
  x1 = scale(X1, T, F),
  x2 = scale(X2, T, F), 
  x3 = scale(X3, T, F))

fit <- lm(Y ~ (X1 + X2 + X3)^2, .data)

## ------------------------------------------------------------------------
cor(model.matrix(fit)[, -1])         # Interactions, high correlation (p 312)

fit <- lm(Y ~ (x1 + x2 + x3)^2, .data)  # fit centered values instead
cor(model.matrix(fit)[, -1])         # less correlation, but still present
summary(fit)
anova(fit)
anova(lm(Y ~ x1 + x2 + x3, .data), fit) # partial F-test: F = 0.53 < qf(.95, 3,13)

## ------------------------------------------------------------------------
.data <- CH08TA02
names(.data) <- c("Y", "X1", "X2")
.data <- transform(.data, X2 = factor(X2, labels = c("Mutual", "Stock")))
fit <- lm(Y ~ X1 + X2, .data)

## ------------------------------------------------------------------------
plot(Y ~ X1, .data, pch = 19, ylim = c(0, 40), xlim = c(0, 310),
     xlab = "Size of Firm", ylab = "Months Elapsed")

abline(coef(fit)[[1]], coef(fit)[[2]],)  # Mutual Firm
abline(coef(fit)[[1]] + coef(fit)[[3]], coef(fit)[[2]], col = "red")  # Stock  Firm

## ------------------------------------------------------------------------
data.frame(
    "Months"      = .data$Y,
    "Firm Size"   = .data$X1,
    "Firm Type"   = .data$X2,
    "Indicator"   = model.matrix(fit)[, 3],
    "Interaction" = .data$X1 * (unclass(.data$X2) - 1)
)

## ------------------------------------------------------------------------
summary(fit)
anova(fit)
confint(fit)

## ------------------------------------------------------------------------
fit <- lm(Y ~ X1 + X2, .data)
plot(Y ~ X1, .data, col = X2, pch = 19, ylim = c(0, 40), xlim = c(0, 310),
     xlab = "Size of Firm", ylab = "Months Elapsed")
abline(coef(fit)[[1]], coef(fit)[[2]])  # Mutual Firm
abline(coef(fit)[[1]] + coef(fit)[[3]],  coef(fit)[[2]], col = "red")  # Stock  Firm

## ------------------------------------------------------------------------
fit <- lm(Y ~ X2/X1 - 1, .data)  # Using nested model for ease of presentation
plot(Y ~ X1, .data, col = X2, pch = 19, ylim = c(0, 40), xlim = c(0, 310),
     xlab = "Size of Firm", ylab = "Months Elapsed")
abline(coef(fit)[[1]], coef(fit)[[3]])               # Mutual Firm
abline(coef(fit)[[2]], coef(fit)[[4]], col = "red")  # Stock  Firm

## ------------------------------------------------------------------------
fit <- lm(Y ~ X1 * X2, .data)  # Going back to intended model
summary(fit)
anova(fit)
anova(lm(Y ~ X1 + X2, .data), fit)

## ------------------------------------------------------------------------
.data <- CH08TA05
names(.data) <- c("Y", "X1", "X2")
.data <- transform(.data, X2 = factor(X2, labels = c("Line1", "Line2")))

## ------------------------------------------------------------------------
data.frame(
  "Scrap"      = .data$Y,
  "Line_Speed" = .data$X1,
  "Indicator"  = .data$X2)

## ------------------------------------------------------------------------
fit <- lm(Y ~ X2/X1 - 1, .data)  # Using nested model for ease of presentation
plot(Y ~ X1, .data, col = X2, ylim = c(100, 500), xlim = c(100, 350), pch = 19,
     xlab = "Line Speed", ylab = "Amount of Scrap")
abline(coef(fit)[[1]], coef(fit)[[3]])               # Production Line 1
abline(coef(fit)[[2]], coef(fit)[[4]], col = "red")  # Production Line 2

## ------------------------------------------------------------------------
plot(fitted(fit), resid(fit), col = .data$X2, pch = 19,
     xlim = c(100, 500), ylim = c(-40,40),
     ylab = "Residual", xlab = "Fitted")
title("(a) Production Line 1")
abline(0, 0)

plot(jitter(as.numeric(.data$X2), 0.1), resid(fit), pch = 19, 
     xlab = "Production Line", ylab = "Residual", sub = "jittered added to X2")

qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
qqline(resid(fit))
title("Normal Probability Plot")

## ------------------------------------------------------------------------
fit <- lm(Y ~ X1 * X2, .data)  # No longer using the nested model
summary(fit)
anova(fit)

