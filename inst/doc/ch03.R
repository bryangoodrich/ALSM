## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(MASS)
library(car)
library(lmtest)

data("CH01TA01", package = "ALSM")
data("CH03TA01", package = "ALSM")
data("CH03TA04", package = "ALSM")
data("CH03TA07", package = "ALSM")
data("CH03TA08", package = "ALSM")
data("CH03TA10", package = "ALSM")

## ------------------------------------------------------------------------
.data  <- CH01TA01
names(.data) <- c('x', 'y')
fit <- lm(y ~ x, .data)

## ------------------------------------------------------------------------
par(mfrow = c(2, 2), pch = 19)
stem(.data$x, scale=3)  # Printed to screen, not graphics window

dotchart(.data$x, xlab = "Lot Size", main = "(a) Dot Plot")

plot(.data$x, type = "b", lty = 2, xlab = "Run", ylab = "Lot Size")
title("(b) Sequence Plot")

boxplot(.data$x, horizontal = TRUE, xlab = "Lot Size", main = "(d) Box Plot")

## ------------------------------------------------------------------------
par(mfrow = c(2, 2), pch = 19)

plot(.data$x, resid(fit), xlab = "Lot Size", ylab = "Residual")
title("(a) Residual Plot against x")

plot(resid(fit), type = "b", lty = 2, xlab = "Run", ylab = "Residual")
title("(b) Sequence Plot")

boxplot(resid(fit), horizontal = TRUE, xlab = "Residual")
title("(c) Box Plot")

qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
title("(d) Normal Probability Plot")
qqline(resid(fit))

## ------------------------------------------------------------------------
.data <- CH03TA01
names(.data) <- c('y', 'x')
fit <- lm(y ~ x, .data)

## ------------------------------------------------------------------------
plot(y ~ x, .data, ylim = c(0, 8), main = "(a) Scatter Plot")
abline(fit)

plot(.data$x, resid(fit), main = "(b) Residual Plot")
abline(0, 0)

## ------------------------------------------------------------------------
cbind(
  "Increase in Ridership" = .data$y,
  "Maps Distributed"      = .data$x,
  "Fitted Values"         = round(fitted(fit), 2),
  "Residuals"             = round(resid(fit), 2))

## ------------------------------------------------------------------------
.data  <- CH01TA01
names(.data) <- c('x', 'y')
fit <- lm(y ~ x, .data)

## ------------------------------------------------------------------------
# Expected value under normality comes from equation (3.6)
cbind(
  "Residual"                   = round(resid(fit), 2),
  "Rank"                       = rank(resid(fit)),
  "Exp. Value under Normality" = round(sqrt(deviance(fit) / df.residual(fit)) * 
                                       qnorm((rank(resid(fit)) - 0.375) / (nrow(data) + 0.25)), 2))

## ----message=FALSE-------------------------------------------------------
levene <- function(model, alpha = 0.05)
{
  f <- function(x) 
  {
    within(x, {
      dsqdev = scale(abs(e - median(e)), T, F)^2
      d      = abs(e - median(e))
    })
  }  # end f
  
  data <- model.frame(model)  # Grab the data used in the model
  data <- transform(data,     # Append the residuals and splitting factor
    e = resid(model),
    group = cut(x, 2, labels = LETTERS[1:2]))

  # Split by group factor and add last TABLE 3.3 columns. Notice the syntax for within()
  data.split <- with(data, split(subset(data, select = -group), group))
  data.split <- lapply(data.split, f)

  
  # Define the relevant variables for hypothesis test and return object
  SSd     <- lapply(data.split, function(x) sum(x$dsqdev))
  dbar    <- lapply(data.split, function(x) mean(x$d))
  n       <- c(Total = nrow(data), lapply(data.split, nrow))   # List of n's
  s       <- sqrt( (SSd$A + SSd$B) / (n$Total - 2) )           # sqrt of (3.9a)
  tstar   <- (dbar$A - dbar$B) / (s * sqrt(1/n$A + 1/n$B))     # (3.9)
  tval    <- qt(1 - alpha/2, n$Total - 2)
  p.value <- 2 * pt(-abs(tval), n$Total - 1)                   # = 0.0495
  
  # Print conclusion of the decision rule
  if (abs(tstar) <= tval) {
    print("error variance is constant")
  } else 
    print("error variance is not constant")

  # return split data and defined variables
  data <- list(
    'data'    = data.split,
    'p.value' = p.value,
    'tstar'   = tstar,
    'tval'    = tval,
    'dbar'    = unlist(dbar),
    'SSd'     = unlist(SSd), 
    'n'       = unlist(n),
    's'       = s)

  return(data)
}  # end levene function

levene(fit)
with(.data, leveneTest(y, cut(x, 2)))  # F value is within acceptance region

## ----message=FALSE-------------------------------------------------------
bptest(fit, student = FALSE)
qchisq(1 - 0.05, 1) # BP falls within the chisq acceptance region

## ------------------------------------------------------------------------
.data <- CH03TA04
names(.data) <- c('x', 'y')
fit <- lm(y ~ x, .data)

## ------------------------------------------------------------------------
as.table(cbind(
  "Size of Min. Deposit"   = .data$x,
  "Number of New Accounts" = .data$y
))

anova(fit)

## ------------------------------------------------------------------------
plot(y ~ x, .data)
abline(fit)

## ------------------------------------------------------------------------
tab           <- do.call("cbind", unstack(.data, y ~ x))
tab[2, 4]     <- NA  # The singular "150" level value is duplicated; remove it
tab           <- rbind(tab, colMeans(tab, na.rm = TRUE))
dimnames(tab) <- list(Replicate = c(1:2, "Mean"), Deposit = dimnames(tab)[[2]])
(tab          <- as.table(tab))

## ------------------------------------------------------------------------
fit.aov <- anova(update(fit, . ~ . + factor(x)))
as.table(cbind(
  'SS' = c('SSR'   =     fit.aov[1,   2], 
           'SSE'   = sum(fit.aov[2:3, 2]),
           'SSLF'  =     fit.aov[2,   2],
           'SSPE'  =     fit.aov[3,   2],
           'Total' = sum(fit.aov[1:3, 2])), 

  'Df' = c(              fit.aov[1,   1],
                     sum(fit.aov[2:3, 1]),
                         fit.aov[2,   1],
                         fit.aov[3,   1],
                     sum(fit.aov[1:3, 1])),  

  'MS' = c(              fit.aov[1,   3],
                     sum(fit.aov[2:3, 2]) / sum(fit.aov[2:3, 1]),
                         fit.aov[2,   3],
                         fit.aov[3,   3],
                         NA)  
))

## ------------------------------------------------------------------------
.data <- CH03TA07
names(.data) <- c('x', 'y')
fit <- lm(y ~ x, .data)

## ------------------------------------------------------------------------
(.data <- transform(.data, sqrtx = sqrt(x)))

## ------------------------------------------------------------------------
par(mfrow = c(2, 2), pch = 19)
plot(y ~ x, .data, xlab = "Days", ylab = "Performance")
title("(a) Scatter Plot")

plot(y ~ sqrt(x), .data, xlab = expression(sqrt(x)), ylab = "Performance")
title(expression(paste("(b) Scatter Plot against ", sqrt(x))))

plot(resid(fit) ~ sqrt(x), .data, xlab = expression(sqrt(x)), ylab = "Residual")
title(expression(paste("(c) Residual Plot against ", sqrt(x))))

qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
title("(d) Normal Probability Plot")

## ------------------------------------------------------------------------
.data <- CH03TA08
names(.data) <- c("x", "y", "z")
fit <- lm(z ~ x, .data)

## ------------------------------------------------------------------------
with(.data, data.frame("Age" = x, "Plasma" = y, "Transform" = z))

## ------------------------------------------------------------------------
par(mfrow = c(2, 2), pch = 19)
plot(y ~ x, .data, xlab = "Age", ylab = "Plasma Level")
title("(a) Scatter Plot")

plot(z ~ x, .data, xlab = "Age")
title("(b) Scatter Plot with Y' = log(Y)")

plot(.data$x, resid(fit)*100, xlab = "Age", ylab = "Residual x 100")
title("(c) Residual Plot against X")
abline(0, 0)

qqnorm(resid(fit), xlab = "Expected Valued", ylab = "Residual", main = "")
title("(d) Normal Probability Plot")

## ----box-cox, message=FALSE, fig.width=10, fig.height=4------------------
par(mfrow = c(1, 3))
boxcox(lm(y ~ x, .data))
boxCox(lm(y ~ x, .data))   
boxcox.sse <- function(lambda, model)
{
  x  <- model.frame(model)$x
  y  <- model.frame(model)$y
  K2 <- prod(y)^( 1/length(y) )            # (3.36a)
  K1 <- 1 / (lambda * K2^(lambda - 1))     # (3.36b)
  ifelse(lambda != 0,                      # (3.36)
    assign("W", K1 * (y^lambda - 1)),
    assign("W", K2 * log(y)))

  # Deviance = Residual Sum of Squares
  return(deviance(lm(W ~ x)))  
}

lambda <- seq(-2, 2, by = 0.1)
SSE = sapply(lambda, boxcox.sse, lm(y ~ x, .data))
plot(lambda, SSE, type = "l", xlab = expression(lambda))
abline(v = -0.5, lty = 3)
cbind('lambda' = lambda, 'SSE' = SSE)

## ------------------------------------------------------------------------
.data <- CH01TA01
names(.data) <- c('x', 'y')
fit  <- lm(y ~ x, .data)

## ------------------------------------------------------------------------
with(.data, scatter.smooth(x, y))
title("Lowess Curve and Linear Regression Confidence Bands")

plot(y ~ x, .data, xlab = "Lot Size", ylab = "Hours")
title("Lowess Curve and Linear Regression Confidence Bands")
with(.data, lines(loess.smooth(x, y), col = "red"))

# Gather confidence bands, ordered by x, and add lines to plot
ci <- cbind(model.frame(fit), predict(fit, int = "c"))[order(.data$x), ]
lines(lwr ~ x, ci, col = "blue", lty = "dashed" )
lines(upr ~ x, ci, col = "blue", lty = "dashed" )

## ------------------------------------------------------------------------
.data <- CH03TA10
names(.data) <- c("y", "x")

## ------------------------------------------------------------------------
cbind("Plutonium Activity" = .data$x, "Alpha Count Rate" = .data$y)

## ----warning=FALSE-------------------------------------------------------
plot(y ~ x, .data, pch = 19, xlab = "pCi/g", ylab = "#/sec")
with(.data, lines(loess.smooth(x, y)))  # Warnings may occur

c(Stat  = bptest(lm(y ~ x, .data), student = FALSE)$statistic,
  Chisq = pchisq(.95, 1))  # BP falls outside acceptance region: Reject H0

## ------------------------------------------------------------------------
.data <- .data[-24, ]  # Remove outlier: Record 24
.data <- transform(.data, sqrty = sqrt(y), sqrtx = sqrt(x))

# Linear Models Summay and Anova
summary(lm(y ~ x, .data))
anova(lm(y ~ x, .data))
summary(lm(sqrty ~ x, .data))
anova(lm(sqrty ~ x, .data))
summary(lm(sqrty ~ sqrtx, .data))
anova(lm(sqrty ~ sqrtx, .data)) 

pplot <- function(formula, data) 
{
  fit <- lm(formula, data)
  plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual")
  title("(b) Residual Plot")
  abline(0, 0)
  
  qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
  qqline(resid(fit))
  title("(c) Normal Probability Plot")
}

par(mfrow = c(2, 2), pch = 19)
pplot(y ~ x, .data)          # Untransformed Plots (FIGURE 3.21)
pplot(sqrty ~ x, .data)      # Transformed Response Plots (FIGURE 3.22)
pplot(sqrty ~ sqrtx, .data)  # Transformed Response and Predictor Plots (FIGURE 3.23)


par(mfrow = c(1,1))

fit <- lm(sqrty ~ sqrtx, .data)
scatter.smooth(.data$sqrtx, .data$sqrty, main = "(d) Confidence Band", pch = 19,
               xlab = expression(sqrt(x)), ylab = expression(sqrt(y)))

newdata = data.frame(sqrtx = sort(.data$sqrtx))
pp <- predict(fit, newdata, int = "c")
lines(newdata$sqrtx, pp[, 2], col = "blue")
lines(newdata$sqrtx, pp[, 3], col = "blue")

