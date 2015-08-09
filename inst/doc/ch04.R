## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
data("CH01TA01", package = "ALSM")
data("CH04TA02", package = "ALSM")

## ------------------------------------------------------------------------
.data  <- CH01TA01
names(.data) <- c('x', 'y')
fit  <- lm(y ~ x, .data)

## ------------------------------------------------------------------------
confint(fit, level = (1 - 0.1/2))

## ------------------------------------------------------------------------
ci.wh <- function(model, newdata, alpha = 0.1)
{
  df    <- nrow(model.frame(model)) - length(coef(model))  # 23
  W     <- sqrt( 2 * qf(1 - alpha, 2, df) )                # 2.2580
  ci    <- predict(model, newdata, se.fit = TRUE)   
  x <- cbind(
    'x'   = newdata,
    's'   = ci$se.fit,
    'fit' = ci$fit,
    'lwr' = ci$fit - W * ci$se.fit,
    'upr' = ci$fit + W * ci$se.fit)
  
  return(x)
}

new <- data.frame(x = c(30, 65, 100))
ci.wh(fit, new)

## ------------------------------------------------------------------------
predict(fit, new, int = "c", level = (1 - 0.1/nrow(new)), se.fit = TRUE)

## ------------------------------------------------------------------------
ci.sim <- function(model, newdata, type = c("B", "S"), alpha = 0.05)
{
  g  <- nrow(newdata)
  CI <- predict(model, newdata, se.fit = TRUE)
  M  <- ifelse(match.arg(type) == "B",
          qt(1 - alpha / (2*g), model$df),              # B = (4.9a)
          sqrt( g * qf( 1 - alpha, g, model$df)))       # S = (4.8a)

  spred <- sqrt( CI$residual.scale^2 + (CI$se.fit)^2 )  # (2.38) 
  x <- data.frame(
    "x"     = newdata,
    "spred" = spred,
    "fit"   = CI$fit,
    "lower" = CI$fit - M * spred,
    "upper" = CI$fit + M * spred)
  
  return(x)
}

new <- data.frame(x = c(80, 100))
ci.sim(fit, new, type = "S")
ci.sim(fit, new, type = "B")

## ------------------------------------------------------------------------
.data <- CH04TA02
names(.data) <- c('x', 'y')

## ------------------------------------------------------------------------
fit <- lm(y ~ 0 + x, .data)   # alternate way to exclude intercept: lm(y ~ x - 1)
tab <- transform(.data, 
  x   = round(x, 0),
  y   = round(y, 0),
  xy  = round(x * y, 0),
  xsq = round(x * x, 0),
  fit = round(fitted(fit), 2),
  e   = round(resid(fit), 2))

rbind(tab, Total = colSums(tab))

plot(.data$x, fitted(fit), pch = 19, xlab = "Work United Performed", ylab = "Variable Labor Costs")
abline(fit)

# Some of the results from pp 163-4 are included or derivable from these
summary(fit)
anova(fit)
confint(fit)

