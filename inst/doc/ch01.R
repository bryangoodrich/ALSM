## ----set-global-opts, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)
options(width = 100)

## -------------------------------------------------------------------------------------------------
data("CH01TA01", package = "ALSM")

## -------------------------------------------------------------------------------------------------
.data        <- CH01TA01     # Assign a data set to an R object
names(.data) <- c("x", "y")  # Give useful, albeit arbitrary, names to our variables

## -------------------------------------------------------------------------------------------------
xx <- scale(.data[, 'x'], TRUE, FALSE)  # Center the x values
yy <- scale(.data[, 'y'], TRUE, FALSE)  # Center the y values

# Create a table for our output
tab <- transform(.data,                 
  x      = x,
  y      = y,
  xdif   = xx,
  ydif   = yy,
  crp    = round(xx * yy),
  sqdevx = xx^2,
  sqdevy = round(yy^2))

# Append summary rows to our table
tab <- rbind(                       
  tab,
  Total = round(apply(tab, 2, sum), 0),
  Mean  = c(colMeans(.data), rep("", 5)))

print(tab, quote = FALSE)

## -------------------------------------------------------------------------------------------------
fit <- lm(y ~ x, data = .data)  # Fit the linear regression model
summary(fit)                 # Output summary information for this model

plot(y ~ x, .data,  xlab = "Lot Size", ylab = "Hours", pch = 20)
title("(a) Scatter Plot")  # Add a title to the current device plot

# notice the different ways by which to issue plot commands
with(.data, plot(x, y, xlab = "Lot Size", ylab = "Hours", pch = 20))
title("(b) Fitted Regression Line")  # Add a title to the current device plot
abline(fit)                          # Add a trend line for this fitted model

## -------------------------------------------------------------------------------------------------
tab <- cbind(
  "Lot Size (X)"       = .data$x,
  "Work Hours (Y)"     = .data$y,
  "Est. Mean Response" = round(fitted(fit),   2),
  "Residuals"          = round( resid(fit),   2),
  "Sq.Residuals"       = round( resid(fit)^2, 1))

rbind(tab, Totals = colSums(tab))

