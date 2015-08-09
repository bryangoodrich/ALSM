## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(MPV)  # For PRESS wrapper function
library(leaps)  # bestsubsets

data("CH09TA01", package = "ALSM")
data("CH09TA05", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH09TA01
names(.data) <- c(paste("X", 1:8, sep = ""), "Y" , "lnY")
data.names <- c("Bloodclot", "Progindex", "Enzyme", "Liver",
                "Age", "Gender", "Alc.Mod", "Alc.Heavy",
                "Survival", "LnSurvival")

## ------------------------------------------------------------------------
{
  print(head(.data))
  cat("... \n")
  print(tail(.data))
}

## ------------------------------------------------------------------------
with(.data,
{
    stem(X1, 2)
    stem(X2, 4)
    stem(X3, 4)
    stem(X4)
    cor(.data[1:4])
})

pairs(.data[1:4], labels = data.names[1:4])

with(.data, 
{
    par(mfrow = c(2, 2))
    boxplot(X1, main = "Bloodclot")
    boxplot(X2, main = "Progindex")
    boxplot(X3, main = "Enzyme"   )
    boxplot(X4, main = "Liver"    )
})

## ------------------------------------------------------------------------
par(mfrow = c(2, 2))

fit <- lm(Y ~ X1 + X2 + X3 + X4, .data)
plot(fitted(fit), resid(fit), xlab = "Predicted value", ylab = "Residual")
title("(a) Residual Plot for Y")
abline(0, 0)

qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
qqline(resid(fit))
title("(b) Normal Plot for Y")

fit <- lm(lnY ~ X1 + X2 + X3 + X4, .data)
plot(fitted(fit), resid(fit), xlab = "Predicted value", ylab = "Residual")
title("(c) Residual Plot for lnY")
abline(0, 0)

qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
qqline(resid(fit))
title("(d) Normal Prot for lnY")

## ------------------------------------------------------------------------
cor(model.frame(fit))
pairs(model.frame(fit), labels = data.names[c(10, 1:4)], pch = 18, col = "grey50")

## ----best, message=FALSE-------------------------------------------------
best <- function(model, ...) 
{
  subsets <- regsubsets(formula(model), model.frame(model), ...)
  subsets <- with(summary(subsets),
                  cbind(p = as.numeric(rownames(which)), which, rss, rsq, adjr2, cp, bic))
  
  return(subsets)
}  

round(best(fit, nbest = 4), 4)

subsets <- regsubsets(formula(fit), model.frame(fit), nbest = 4)

plot(subsets)
plot(subsets, scale = "Cp")
plot(subsets, scale = "adjr2")
plot(subsets, scale = "r2")

## ------------------------------------------------------------------------
x <- best(fit, nbest = 6)

par(mfrow = c(2, 2), pch = 19)

plot(rsq ~ p, x,   xlab = "(a)", ylab = "Rsq", col = "gray50")
lines(1:4, tapply(x[, "rsq"], x[, "p"], max), lwd = 2)

plot(adjr2 ~ p, x, xlab = "(b)", ylab = "Adj Rsq", col = "gray50")
lines(1:4, tapply(x[, "adjr2"], x[, "p"], max), lwd = 2)

plot(cp ~ p, x,    xlab = "(c)", ylab = "Cp", col = "gray50")
lines(1:4, tapply(x[, "cp"], x[, "p"], min), lwd = 2)

plot(bic ~ p, x,   xlab = "(d)", ylab = "BIC", col = "gray50")
lines(1:4, tapply(x[, "bic"], x[, "p"], min), lwd = 2)

## ------------------------------------------------------------------------
fit <- lm(lnY ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8, .data)
x   <- best(fit, nbest = 20)

par(mfrow = c(2, 2), pch = 19)
plot(rsq ~ p, x, xlab = "(a)", ylab = "Rsq", col = "gray50")
lines(1:8, tapply(x[, "rsq"], x[, "p"], max), lwd = 2)

plot(adjr2 ~ p, x, xlab = "(b)", ylab = "Adj Rsq", col = "gray50")
lines(1:8, tapply(x[, "adjr2"], x[, "p"], max), lwd = 2)

plot(cp ~ p, x, xlab = "(c)", ylab = "Cp", col = "gray50")
lines(1:8, tapply(x[, "cp"], x[, "p"], min), lwd = 2)

plot(bic ~ p, x, xlab = "(d)", ylab = "BIC", col = "gray50")
lines(1:8, tapply(x[, "bic"], x[, "p"], min), lwd = 2)

## ------------------------------------------------------------------------
round(best(fit, nbest = 2), 4)  # Fundamentally same info as in TABLE 9.3

## ------------------------------------------------------------------------
BIC <- log(nrow(model.frame(fit)))
fit <- lm(lnY ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8, .data)  # Full Model
step(fit, direction = "forward")      # AIC = -160.77 Drop None
step(fit, dir = "forward",  k = BIC)  # BIC = -142.87 Drop None
step(fit, dir = "backward")           # AIC = -163.83 Drop X4 and X7 
step(fit, dir = "backward", k = BIC)  # BIC = -153.41 Drop X4 through X7 -- This model matches the book's search

## ------------------------------------------------------------------------
# Import the training and validation data sets, respectively
.data <- do.call("rbind.data.frame", list(CH09TA01, CH09TA05))
names(.data) <- c(paste("X", seq(8), sep = ""), "Y", "lnY")

# Add a factor to which data set the observation belongs
.data <- transform(.data, class = gl(2, 54, labels = c("training", "validation")))

## ------------------------------------------------------------------------
{
  print(head(.data))
  cat("... \n")
  print(tail(.data))
}

## ----reg-result, message=FALSE-------------------------------------------
newsummary <- function(model)
{
    list('coefs'    = round(t(summary(model)$coef[, 1:2]), 4),
         'criteria' = cbind('SSE'   = anova(model)["Residuals", "Sum Sq"],
                            'PRESS' = PRESS(model),
                            'MSE'   = anova(model)["Residuals", "Mean Sq"],
                            'Rsq'   = summary(model)$adj.r.squared))
}

newsummary(lm(lnY ~ X1 + X2 + X3 + X8,           .data, class == "training"))
newsummary(lm(lnY ~ X1 + X2 + X3 + X8,           .data, class == "validation"))

newsummary(lm(lnY ~ X1 + X2 + X3 + X6 + X8,      .data, class == "training"))
newsummary(lm(lnY ~ X1 + X2 + X3 + X6 + X8,      .data, class == "validation"))

newsummary(lm(lnY ~ X1 + X2 + X3 + X5 + X6 + X8, .data, class == "training"))
newsummary(lm(lnY ~ X1 + X2 + X3 + X5 + X6 + X8, .data, class == "validation"))

