## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(SuppDists)  # qmaxFratio
library(HH)         # hov
library(MASS)       # boxcox

data("CH17TA02", package = "ALSM")
data("CH18TA02", package = "ALSM")
data("CH18TA05", package = "ALSM")
data("CH18TA07", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH17TA02
names(.data) <- c("y", "x1", "x2")
fit <- lm(y ~ factor(x1)-1, .data)

## ------------------------------------------------------------------------
xtabs(resid(fit) ~ x2 + x1, .data)

## ------------------------------------------------------------------------
plot(resid(fit) ~ fitted(fit), ylab = "Residual", xlab = expression(hat(Y)), pch = 19)
abline(0,0)
title(expression(paste("(a) Residual against ", hat(Y))))

stripchart(split(resid(fit), .data$x1), method = "stack",  pch = 19)
abline(h = seq(2, 4)-0.1)
title("(b) Aligned Residual Dot Plot")

qqnorm(resid(fit), xlab = "Exp Val", ylab = "Residual", pch = 19, main = "")
qqline(resid(fit))
title("(c) Normal Probability Plot")

## ------------------------------------------------------------------------
.data <- CH18TA02
names(.data) <- c("y", "x1", "x2")

## ----fratio, message=FALSE-----------------------------------------------
tab <- xtabs(y ~ x2 + x1, .data)
round(addmargins(tab, 1, FUN = list(list(mean, median, s = var))), 3)

vars <- with(.data, tapply(y, x1, var))  # (18.10)
cat("H* = ", max(vars) / min(vars))   # (18.8)
qmaxFratio(0.95, 7, 5)                # (18.9)

## ------------------------------------------------------------------------
stripchart(round(y) ~ x1, .data, pch = 19, method = "stack", xlim = c(0, 30), xlab = "Pull Strength", ylab = "Type")
abline(h = seq(2,5) - 0.1)

## ---- message=FALSE------------------------------------------------------
hov(y ~ factor(x1), .data)
data.frame(with(.data, tapply(y, x1, function(x) list(abs(x - median(x))))))

## ------------------------------------------------------------------------
w = 1 / vars          # Factor level variances (18.10) defined earlier
w = rep(w, each = 8)  # (18.14) -- Repeat each weight by # of x2 factor levels
fit <- lm(y ~ x1 -1, transform(.data, x1 = factor(x1)), weights = w)  # (18.17)

data.frame(
  'i' = .data$x1,
  'j' = .data$x2,
  'Y' = .data$y,
  model.matrix(fit),
  'Weights' = w,
  'Reduced.Model' = 1)

summary(fit)
anova(fit)
summary(lm(y ~ 1, .data, weights = w))  # Reduced Model (18.19)
anova(lm(y ~ 1, .data, weights = w))
anova(lm(y ~ 1, .data, weights = w), fit)  # Confirm that x1 is significant to keep

## ------------------------------------------------------------------------
.data <- CH18TA05
names(.data) <- c("y", "x1", "x2")

## ------------------------------------------------------------------------
addmargins(xtabs(y ~ x2 + x1, .data), 1, list(list(mean, var)))
addmargins(xtabs(rank(y) ~ x2 + x1, .data), 1, list(list(mean, var)))

cbind(
  '(1)' = with(.data, tapply(y, x1, function(x) var(x) / mean(x))),
  '(2)' = with(.data, tapply(y, x1, function(x) sd(x)  / mean(x))),
  '(3)' = with(.data, tapply(y, x1, function(x) sd(x)  / mean(x)^2)))

## ----boxcox, message=FALSE-----------------------------------------------
fit <- lm(y ~ factor(x1) - 1, .data)
boxcox(fit)

## ----prob-plots, fig.width=12, fig.height=6------------------------------
with(.data, tapply(log(y), x1, var))   # More stable variances (p 792)
with(.data, tapply(log(y), x1, mean))  # Transformed means (p 793)

par(mfrow = c(1, 2), pch = 16)
qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
qqline(resid(fit))
title("(a) Original Data")

fit <- update(fit, log(.) ~ .)
qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
qqline(resid(fit))
title("(b) Transformed Data")

  ## Some statistics (p 793)
anova(lm(log(y) ~ 1, .data), fit)  # anova(reduced, full)
qf(0.90, 2, 12)                 # Compare with above F*
summary(fit)                    # Authors typo mean(location2) = 2.792
with(.data, pairwise.t.test(log(y), x1, p.adj = "bonf", conf.level = 0.9))
TukeyHSD(aov(log(y) ~ factor(x1), .data), conf.level = 0.90)

## ------------------------------------------------------------------------
fit <- lm(rank(y) ~ factor(x1), .data)
addmargins(xtabs(rank(y) ~ x2 + x1, .data), 1, list(list(mean, var)))  # Same results as earlier
anova(fit)  # Compare F-value with value below
cat("F* =", qf(0.90, 2, 12))
kruskal.test(y ~ x1, .data)  # Compare the chi-squared value with value below
cat("chi* =", qchisq(0.90, 2))


# Multiple Pairwise Comparisons
B = qnorm(1 - 0.1 / 6) * sqrt((15 * (15 + 1)) / 12 * (1/5 + 1/5))
m = with(.data, tapply(rank(y), x1, mean))
round(rbind(
  "1-2" = c('diff' = m[[1]] - m[[2]], (m[1] - m[2]) + c('lwr' = -B, 'upr' =  B)),
  "3-2" = c('diff' = m[[3]] - m[[2]], (m[3] - m[2]) + c('lwr' = -B, 'upr' =  B)),
  "3-1" = c('diff' = m[[3]] - m[[1]], (m[3] - m[1]) + c('lwr' = -B, 'upr' =  B))), 1)

## ------------------------------------------------------------------------
.data <- CH18TA07
names(.data) <-  c("y", "x")
.data <- transform(.data, x = factor(x, labels = c("Low", "Medium", "High")))

## ------------------------------------------------------------------------
with(.data, split(y, x))

## ----diag-plot, message=FALSE--------------------------------------------
fit <- lm(y ~ x - 1, .data)

stripchart(y ~ x, .data, method = "jitter", pch = 19, xlab = "Survival Time")
title("(a) Dot Plots of Survival Times")

stripchart(rstudent(fit) ~ x, .data, method = "jitter", pch = 19, xlab = "Studentized Residual")
title("(b) Dot Plots of Studentized Residuals")
  
qqnorm(rstudent(fit), pch = 19, xlab = "Expected Value", ylab = "Studentized Residual", main = "")
qqline(rstudent(fit))
title("(c) Normal Probability Plot\n of Studentized Residuals")

hov(y ~ x, .data)  # (p 799) test for constancy of variance
boxcox(fit)

## ------------------------------------------------------------------------
fit <- lm(log(y) ~ x, .data)
summary(fit)

stripchart(rstudent(fit) ~ x, .data, method = "jitter", pch = 19, xlab = "Studentized Residual")
title("(a) Dot Plots of Studentized Residuals")
  
qqnorm(rstudent(fit), pch = 19, xlab = "Expected Value", ylab = "Studentized Residual", main = "")
qqline(rstudent(fit))
title("(b) Normal Probability Plot\n of Studentized Residuals")

