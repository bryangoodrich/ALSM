## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(car)  # avPlots, residualPlots

data("CH10TA01", package = "ALSM")
data("CH07TA01", package = "ALSM")
data("CH09TA01", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH10TA01
names(.data) <- c("X1", "X2", "Y")
fit <- lm(Y ~ ., .data)

## ------------------------------------------------------------------------
cbind(
  "Avg Income"     = .data$X1,
  "Risk Aversion"  = .data$X2,
  "Life Insurance" = .data$Y)

## ----avplot, fig.height=4, fig.width=7-----------------------------------
A <- resid(lm(Y ~ X2, .data))
B <- resid(lm(X1 ~ X2, .data))
par(mfrow = c(1, 2), pch = 19)

plot(resid(fit) ~ X1, .data, xlab = "X1", ylab = "Residual")
title("(a) Residual Plot against X1")
abline(0, 0, lty = 2)

plot(A ~ B, xlab = "e(X1|X2)", ylab = "e(Y|X2)")
title("(b) Added-Variable Plot for X1")
abline(lm(A ~ B))
abline(0, 0, lty = 2)

## ------------------------------------------------------------------------
.data <- CH07TA01
names(.data) <- c("X1", "X2", "X3", "Y")
fit <- lm(Y ~ X1 + X2, .data)

## ----resid-plot, message=FALSE, fig.height=4, fig.width=7----------------
avPlots(fit)
residualPlots(fit, fitted = FALSE, quadratic = FALSE, test = FALSE)

## ------------------------------------------------------------------------
fit <- lm(Y ~ X1 + X2, data =
  data.frame(
    'Y' = c(301, 327, 246, 187),
    'X1' = c(14, 19, 12, 11),
    'X2' = c(25, 32, 22, 15)))

cbind(
  'X1'   = model.frame(fit)$X1,
  'X2'   = model.frame(fit)$X2,
  'Y'    = model.frame(fit)$Y,
  'Yhat' = fitted(fit), 
  'e'    = resid(fit),
  'h'    = hatvalues(fit),
  's'    = deviance(fit) * (1 - hatvalues(fit)))

## ------------------------------------------------------------------------
plot(X2 ~ X1, model.frame(fit), pch = 19, xlim = c(10, 20), ylim = c(10, 37),
     xlab = "X1", ylab = "X2", main = "FIGURE 10.6")
text(X2 ~ X1, model.frame(fit), labels = round(hatvalues(fit), 4), pos = 3)

points(mean(X2) ~ mean(X1), model.frame(fit), pch = 22)
text(mean(X2) ~ mean(X1), model.frame(fit), 
     cex = 0.75, pos = 4, labels = "(X1bar, X2bar)")

## ------------------------------------------------------------------------
fit <- lm(Y ~ X1 + X2, .data)
cbind(
  'e' = round(resid(fit), 3),
  'h' = round(hatvalues(fit), 3),
  't' = round(rstudent(fit), 3))

## ------------------------------------------------------------------------
plot(X2 ~ X1, .data, type = 'n', xlim = c(13, 33), ylim = c(40, 60),
     xlab = "Triceps Skinfold Thickness", ylab = "Thigh Circumference")
text(X2 ~ X1, .data, labels = seq(nrow(.data)), cex=.75)

## ------------------------------------------------------------------------
cbind(
  "DFFITS"  = round(dffits(fit), 4),
  "D"       = round(cooks.distance(fit), 4),
  "DFBETA0" = round(dfbetas(fit)[,1], 4),
  "DFBETA1" = round(dfbetas(fit)[,2], 4),
  "DFBETA2" = round(dfbetas(fit)[,3], 4)) 

## ----influence, fig.height=4, fig.width=7--------------------------------
par(mfrow = c(1, 2), pch = 19)
plot(resid(fit) ~ fitted(fit), cex = cooks.distance(fit)*10,
     xlim = c(10, 30), ylim = c(-4.5, 4.5),
     xlab = "YHAT", ylab = "Residual")
title("(a) Proportional Influence Plot")
plot(seq(nrow(.data)), cooks.distance(fit), type = "o", lwd = 2,
     xlab = "Case Index Number", ylab = "Cook's Distance D")
title("(b) Index Influence Plot")

## ------------------------------------------------------------------------
fit <- lm(scale(Y) ~ scale(X1) + scale(X2) + scale(X3), .data)
cbind(
  "Beta*" = round(as.vector(coef(fit)[-1]), 4),
  "VIF"   = round(vif(lm(Y ~ ., .data)), 2))

## ------------------------------------------------------------------------
.data <- CH09TA01
names(.data) <- c(paste("X", seq(8), sep = ""), "Y" , "lnY")
fit <- lm(lnY ~ X1 + X2 + X3 + X8, .data)
vif(fit)  # For page 412

## ------------------------------------------------------------------------
par(mfrow = c(2, 2), pch = 19)

plot(resid(fit) ~ fitted(fit), ylim = c(-0.7, 0.7),
     xlab = "Predicted Value", ylab = "Residual")
title("(a) Residual Plot against Predicted")

plot(resid(fit) ~ .data$X5, ylim = c(-0.7, 0.7),
     xlab = "X5", ylab = "Residual")
title("(b) Residual Plot against X5")

plot(resid(fit) ~ resid(update(fit, X5 ~ .)), ylim = c(-0.7, 0.7),
     xlab = "e(X5|X1238)", ylab = "e(Y'|X1238)")
title("(c) Added-Variable Plot for X5")
abline(lm(resid(fit) ~ resid(update(fit, X5 ~ .))) )

qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
qqline(resid(fit))
title("(d) Normal Probability Plot")

## ------------------------------------------------------------------------
par(mfrow = c(2, 2), pch = 19)

plot(seq(nrow(.data)), rstudent(fit), type = "o",
     main = "(a) Studentized Deleted Residuals",
     xlab = "Case Index", ylab = "t")
plot(seq(nrow(.data)), hatvalues(fit), type = "o",
     main = "(b) Leverage Values",
     xlab = "Case Index", ylab = "h")
plot(seq(nrow(.data)), cooks.distance(fit), type = "o",
     main = "(c) Cook's Distance",
     xlab = "Case Index", ylab = "D")
plot(seq(nrow(.data)), dffits(fit), type = "o",
     main = "(d) DFFITS values",
     xlab = "Case Index", ylab = "DFFITS")

## ------------------------------------------------------------------------
cases <- c(17, 23, 28, 32, 38, 42, 52)
cbind(
  'e'      = round(resid(fit), 4),
  't'      = round(rstudent(fit), 4),
  'h'      = round(hatvalues(fit), 4),
  'D'      = round(cooks.distance(fit), 4),
  'DFFITS' = round(dffits(fit), 4))[cases, ]

