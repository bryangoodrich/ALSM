## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## -----------------------------------------------------------------------------
library(lattice)
library(scatterplot3d)
library(rgl)
library(Rcmdr)

data("CH06FI05", package = "ALSM")

## -----------------------------------------------------------------------------
.data <- CH06FI05
names(.data) <- c("x1", "x2", "y")
fit <- lm(y ~ x1 + x2, .data)

## -----------------------------------------------------------------------------
with(.data, pairs(data.frame("SALES" = y, "TARGETPOP" = x1, "DISPOINC" = x2)))
with(.data, cor  (data.frame("SALES" = y, "TARGETPOP" = x1, "DISPOINC" = x2)))

## -----------------------------------------------------------------------------
cbind(
  "TARGETPOP (X1)" = .data$x1,
  "DISPOINC (X2)"  = .data$x2,
  "SALES (Y)"      = .data$y,
  "FITTED"         = round(fitted(fit),2),
  "RESIDUAL"       = round(resid(fit), 2))

summary(fit)
anova(fit)

## -----------------------------------------------------------------------------
library(lattice)
cloud(y ~ x1 + x2, .data, scales = list(arrows = FALSE), 
      xlab = "TARGETPOP", ylab = "DISPOINC", zlab = "SALES",
      main = "(a) Before Spinning")

cloud(y ~ x1 + x2, .data, scales = list(arrows = FALSE),
      xlab = "TARGETPOP", ylab = "DISPOINC", zlab = "SALES",
      main = "(b) After Spinning", screen = list(z = -40, x = -60, y = 0))

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

# The resolution ('res') parameter controls how 'fine' the plane is
lm3d(fit, res = 30, 
     ticktype = "detailed", shade = 0.5, expand = 2/3,
     xlab = "TARGETPOP", ylab = "DISPOINC", zlab = "SALES",
     theta = 310, phi = 30)

## ----note, eval=FALSE----------------------------------------------------
#  # Change "len = 24" lower or higher for larger or smaller rotations
#  for (t in seq(0, 360, len = 24))
#  {
#    lm3d(fit, 25, ..., theta = t)  # change phi for z-rotation
#    Sys.sleep(0.1)                 # How many seconds to delay before rotating
#  }

## ------------------------------------------------------------------------
par(mfrow = c(2, 2), pch = 19)
plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual")
title("(a) Residual Plot against Y")

plot(resid(fit) ~ x1, .data, xlab = "Targtpop", ylab = "Residual")
title("(b) Residual Plot against X1")

plot(resid(fit) ~ x2, .data, xlab = "Dispoinc", ylab = "Residual")
title("(c) Residual Plot against X2")

plot(resid(fit) ~ I(x1 * x2), .data, xlab = "X1X2", ylab = "Residual")
title("(d) Residual Plot against X1X2")

## ----add-diag, fig.height=4, fig.width=7---------------------------------
par(mfrow = c(1, 2), pch = 19)
plot(abs(resid(fit)) ~ fitted(fit), xlab = "Fitted", ylab = "Absresid")
title("(a) Plot of Absolute Residuals against Y")

qqnorm(resid(fit), main = "", xlab = "Expected", ylab = "Residual")
qqline(resid(fit))
title("(b) Normal Probability Plot")

## ------------------------------------------------------------------------
confint(fit)

## ------------------------------------------------------------------------
predict(fit, data.frame(x1 = 65.4, x2 = 17.6), interval = "confidence")

## ------------------------------------------------------------------------
ci.sim <- function(model, newdata, type = c("B", "S"), alpha = 0.05) 
{
  g  <- nrow(newdata)
  CI <- predict(model, newdata, se.fit = TRUE)
  M  <- ifelse(match.arg(type) == "B",
          qt(1 - alpha / (2*g), model$df),              # B = (4.9a)
          sqrt(g * qf(1 - alpha, g, model$df)))         # S = (4.8a)
  
  spred <- sqrt( CI$residual.scale^2 + (CI$se.fit)^2 )  #     (2.38) 
  x <- data.frame(
    "x"     = newdata,
    "spred" = spred,
    "fit"   = CI$fit,
    "lower" = CI$fit - M * spred,
    "upper" = CI$fit + M * spred)
  
  return(x)
}

newdata <- data.frame( x1 = c(65.4, 53.1), x2 = c(17.6, 17.7) )
ci.sim(fit, newdata, "B", 0.1)  # Bonferroni Prediction
ci.sim(fit, newdata, "S", 0.1)  # Scheffe Prediction

