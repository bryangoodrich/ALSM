# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input Diastolic Blood Pressure (DBP) data
df <- get("CH11TA01", envir = env)
names(df) <- c("X", "Y")



# TABLE 11.1     (p 427)
# Weighted Least Squares--Blood Pressure Example
fit <- lm(Y ~ X, df)
w <- fitted(lm(abs(resid(fit)) ~ X, df))^(-2)
fit <- update(fit, weights = w)
 
cbind(
  "Subject" = seq(nrow(df)),
  "Age"     = df$X,
  "DBP"     = df$Y,
  "e"       = round(resid(lm(Y ~ X, df)), 2),
  "|e|"     = round(abs(resid(lm(Y ~ X, df))), 2),
  "s"       = round(sqrt(1 / w), 4),
  "w"       = round(fit$weights, 5))



# FIGURE 11.1     (p 428)
# Diagnostic Plots Detecting Unequal Variance--Blood Pressure Example
fit <- lm(Y ~ X, df)
par(mfrow = c(1, 3), pch = 19)

plot(Y ~ X, df, xlab = "Age", ylab = "Blood Pressure")
title("(a) Scatter Plot")
abline(fit)

plot(resid(fit) ~ X, df, xlab = "Age", ylab = "Residual")
title("(b) Residual Plot Against X")
abline(0, 0)

plot(abs(resid(fit)) ~ X, df, xlab = "Age", ylab = "Absolute Residual")
title("(c) Absolute Residual Plot Against X")
abline(lm(abs(resid(fit)) ~ X, df))



# Input Body Fat Data (See Ch. 7)
df <- get("CH07TA01", envir = env)
names(df) <- c("x1", "x2", "x3", "y")
fit <- lm(y ~ x1 + x2 + x3, df)



# TABLE 11.2 and TABLE 11.3     (p 434)
# Ridge Estimated Standardized Regression Coefficients for Different
# Biasing Constants c--Body Fat Eample with Three Predictor Variables.
#
# VIF Values for Regression Coefficients and R-squared for Different
# Biasing Constants c--Body Fat Example with Three Predictor Variables.
#
# R has the lm.ridge (MASS) function, but it produces different results.
# I've tried uses methods in other packages, like genridge or ElemStateLearn,
# but the results were pretty much the same. 
c <- c(0, .002, .004, .006, .008, .01, .02, .03, .04, .05, .1, .5, 1)
fun1 <- function(c, rxx, rxy) solve(rxx + c * diag(ncol(rxx)), rxy)  # (11.34)
fun2 <- function(c, rxx) 
{
  bias = solve(rxx + c*diag(ncol(rxx)))
  diag(bias %*% rxx %*% bias)                                        # (11.36)
}  # end fun2

data.frame(lambda = c, 
           t(sapply(c, fun1, rxx = cor(df)[-4, -4], rxy = cor(df)[-4, 4])),
           VIF = t(sapply(c, fun2, rxx = cor(df)[-4, -4])))



# FIGURE 11.3     (p 435)
# Ridge Trace of Estimated Standardized Regression
# Coefficients--Body Fat Example with Three Predictor Variables.
#
# One could also plot the object returned by lm.ridge (in library MASS).
# It does lack the logarithmic scale on x, however. Simply use a code like:
#   plot(lm.ridge(y ~ ., df, lambda = seq(0, 1, by = 0.001)))
lambda <- seq(0.001, 1, by = 0.001)
f <- function(c) {solve(cor(df)[-4, -4] + c * diag(3), cor(df)[-4, 4])}
tab <- data.frame(c = seq(0, 1, by = 0.001), 
                  t(sapply(lambda, f)))

par(lwd = 2)
plot(x1 ~ c, tab, type = "l", log = "x", xaxt = "n", lty = 5,
     xlim = c(0.001, 1), ylim = c(-1, 3), xlab = "", ylab = "")
axis(1, c(0.001, 0.01, .10, 1.00))
lines(x2 ~ c, tab, lty = 3)
lines(x3 ~ c, tab, lty = 1)
abline(0,0)
abline(v = 0.02, lwd = 1)  # Authors choice of c
text(0.02, -1, pos=4, labels = "0.02")



# Input the Mathematics Proficiency Data
df <- get("CH11TA04", envir = env)
names(df) <- c("State", "Y", "X1", "X2", "X3", "X4", "X5")



# FIGURE 11.5     (p 442)
# Comparison of Lowess, Ordinary Least Squares Fits, and
# Robust Quadratic Fits--Mathematics Proficiency Example.
#
# The author confuses the variables of interest. They begin by regressing on
# X2 (HOMELIB), but he talks about X3 (READING). I will maintain using X2 as
# the variable of interest.
#
# Note that there is an R robust regression function rlm (MASS). It uses both
# the huber (psi = psi.huber) and bisquare (psi.bisquare) weighing functions.
# The latter method requires library lqs. For more see
#   http://cran.r-project.org/doc/contrib/Fox-Companion/appendix-robust-regression.pdf
library(MASS)
par(mfrow = c(3, 2), pch = 19)
xlab <- expression(X[2])

# Linear Regression Fit
fit  <- lm(Y ~ X2, df)
outliers <- subset(df, State %in% c("Guam", "D.C.", "Virgin_Islands"))

plot(Y ~ X2, df, xlab = xlab, ylab = "Y")                          # (a)
lines(loess.smooth(df$X2, df$Y), lty = 2)
title("(a) Lowess and Linear Regression Fits")
abline(fit)
with(outliers, text(Y ~ X2, labels = State, pos = 4))

plot(resid(fit) ~ X2, df, xlab = xlab, ylab = "Residual")          # (b)
title("(b) Residuals from Linear Regression")
abline(0, 0)

# OLS Quadratic Fit
fit <- lm(Y ~ x + xsq, transform(df,
  x   = scale(X2, T, F),
  xsq = scale(X2, T, F)^2))
x <- data.frame(X2 = df$X2, Y = fitted(fit))

plot(Y ~ X2, df, xlab = expression(X[2]), ylab = "Y")              # (c)
title("(c) OLS Quadratic Fit")
lines(Y ~ X2, x[order(x$X2), ])

plot(cooks.distance(fit), xlab = "Index", ylab = "D", type = "o")  # (d)
title("(d) Cook's Distance--OLS Quadratic Fit")

# Robust Quadratic Fit
fit <- rlm(formula(fit), model.frame(fit))
x <- data.frame(X2 = df$X2, Y = fitted(fit))

plot(Y ~ X2, df, xlab = expression(X[2]), ylab = "Y")              # (e)
title("(e) Robust Quadratic Fit")
lines(Y ~ X2, x[order(x$X2), ])

plot(fit$w, xlab = "Index", ylab = "W", type = "o")                # (f)
title("(f) Robust Weights")



# TABLE 11.4     (p 443)
# Data Set--Mathematics Proficiency Example.
with(df, cbind("State"    = State,
               "MATHPROF" = Y,
               "PARENTS"  = X1,
               "HOMELIB"  = X2,
               "READING"  = X3,
               "TVWATCH"  = X4,
               "ABSENCES" = X5)) 



# TABLE 11.5     (p 444)
# Iteratively Huber-Reweighted Least Squares
# Calculations--Mathematics Proficiency Example.
#
# This is displayed purely for pedagogical reasons. The final results are
# similar to, per the above example, rlm(formula(fit)). Just compare the
# last column (w7) below with rlm(formula(fit) and model.frame(fit))$w.
fit <- lm(Y ~ x + xsq, transform(df,
  x   = scale(X2, T, F),
  xsq = scale(X2, T, F)^2))

IRLS <- function(model, i = 1) 
{
  e   <- resid(model)
  MAD <- median(abs(e - median(e))) / 0.6745               # (11.46)
  u   <- e / MAD                                           # (11.47)
  w   <- apply(data.frame(u), 1, 
    function(x) if (abs(x) <= 1.345) 1 else 1.345/abs(x))  # (11.44)

  model <- update(model, weights = w)
  if (i > 1) return(IRLS(model, i-1)) else return(model)
}

round(
  cbind("e0" = resid(fit),
        "u0" = resid(fit) / (median(abs(resid(fit) - median(resid(fit)))) / 0.6745),
        "w1" = IRLS(fit, 1)$weights,
        "e1" = resid( IRLS(fit, 1) ),
        "w2" = IRLS(fit, 2)$weights,
        "e2" = resid( IRLS(fit, 2) ),
        "w7" = IRLS(fit, 7)$weights,
        "e7" = resid( IRLS(fit, 7) )),
  4)



# FIGURE 11.6     (p 446)
# Scatter Plot Matrix with Lowess Smooths, and
# Correlation Matrix--Mathematics Proficiency Example.
#
# The use of "pairs" could produce these results, but I want to demonstrate
# the facility of the scatterplotMatrix function (car). It has an "spm" use
# for convenience. Also, the author's plot isn't right for TVWATCH. If you
# plot the scatterplot of it against any of the other variables, you do not
# get the results in this scatterplot matrix. The rest of the scatterplots
# look fine, and so does the correlation table.
library(car)
spm(df[-1], span = 0.75, diagonal = "none", reg.line = FALSE, spread = FALSE,
  var.labels = c("MATHPROF", "PARENTS", "HOMELIB", 
                 "READING", "TVWATCH", "ABSENCES"))
cor(df[-1])



# TABLE 11.6     (p 447)
# Diagnostics for First-Order Model with All Five
# Explanatory Variables--Mathematics Proficiency Example.
fit <- lm(Y ~ X1 + X2 + X3 + X4 + X5, df)
data.frame("State" = levels(df$State),
           "h"     = round(hatvalues(fit), 2),
           "t"     = round(rstudent(fit), 2),
           "D"     = round(cooks.distance(fit), 2))



# UNPRESENTED FIGURE     (p 447)
# Residual plots against each independent variable
# and dependent variable--Mathematics Proficiency Example.
#
# The author indicates that the residual plots present no strong indication
# of nonconstancy of the error variance for the states aside from the
# outliers.
par(mfrow = c(3, 2))
apply(model.frame(fit)[-1], 2, function(x)
  {plot(resid(fit) ~ x, ylab = "Residual"); abline(0, 0)})



# FIGURE 11.7     (p 448)
# Best Subsets Regression--Mathematics Proficiency Example.
#
# In the earlier chapter I used the regsubsets function (leaps) to produce
# the possible submodels and their selection criteria, displaying the best
# based on the number of the "best" parameter. Another approach is to use the
# bestglm function (bestglm). The regsubsets will produce an output similar
# to that found in FIGURE 11.7, but bestglm offers some of its own
# flexibility. See more at the document:
#   http://cran.r-project.org/web/packages/bestglm/vignettes/bestglm.pdf
library(bestglm)
bestglm(model.frame(fit)[c(2:6, 1)], IC = "AIC")$Subsets  # Author's choice is
bestglm(model.frame(fit)[c(2:6, 1)], IC = "BIC")$Subsets  # 3 on each search

summary(rlm(Y ~ X2 + X3 + X4, df))  # (11.53) 
summary(lm (Y ~ X2 + X3 + X4, df))  # (11.54)

with(summary(regsubsets(formula(fit), df, nbest=2)),  # Output similar to
  cbind(which[, -1], rss, rsq, adjr2, cp, bic))        # FIGURE 11.7




# Input Life Insurance Data
df <- get("CH11TA07", envir = env)
names(df) <- c("X1", "X2", "Y")




# TABLE 11.7     (p 451)
# Lowess Calculations for Non-parametric Regression Fit
# at Xh1 = 30, Xh2 = 3--Life Insurance Example.
fit <- lm(Y ~ X1 + X2, df)  # (6.1)
lo <- function(model, newdata = data.frame(X1 = 30, X2 = 3), q = 0.5) 
{
  x<- transform(model.frame(model),
    d = sqrt(scale(X1, newdata$X1, sd(X1))^2 + scale(X2, newdata$X2, sd(X2))^2))
  dq <- with(x, sort(d)[nrow(x) * q])
  x <- transform(x, w = ifelse(d > dq, 0, (1-(d/dq)^3)^3))
  with(x, update(model, weights = w))
}

transform(df, 
  d = sqrt(scale(X1, 30, sd(X1))^2 + scale(X2, 3, sd(X2))^2),
  w = lo(fit)$weights)

predict(lo(fit), data.frame(X1 = 30, X2 = 3))  # 4.664607

# Compare
fit <- loess(Y ~ X1 + X2, df, span = 0.5, degree = 1)
predict(fit, data.frame(X1 = 30, X2 = 3))  # returns 4.2894 with normalize = FALSE



# FIGURE 11.8     (p 452)
# Contour and Conditioning Plots for Lowess Nonparametric
# Regression--Life Insurance Example.
#
# Loess regressions can have a variety of results due to small calibrations
# in the possible parameters. While I cannot fit a model to exactly produce
# what the book expects, I demonstrate how to do contour plots, extend that
# to a 3D surface, and do the associated coplots manually. There does
# exist a `coplot` function, but you will have to figure it out.
fit <- loess(Y ~ X1 + X2, df, span = 0.5, degree = 1)

x1 <- seq(min(df$X1), max(df$X1), len=25)    # Equally spaced sequence of values from X1
x2 <- seq(min(df$X2), max(df$X2), len=25)    # Equal length sequence like above for X2
newdata <- expand.grid(X1=x1, X2=x2)         # Create a matrix grid of all combinations
y  <- matrix(predict(fit, newdata), 25, 25)  # Prediction matrix for all pairs in grid
contour(x1, x2, y)

# Better!
persp(x1, x2, y, theta=270, phi=0,  ticktype="detailed", expand=2/3, shade=0.5)
persp(x1, x2, y, theta=300, phi=30, ticktype='detailed', expand=2/3, shade=0.5)

par(mfrow=c(1,3), pch=20)
plot(x1, predict(fit, data.frame(X1=x1, X2=3)), type='l', ylim=c(0, 360))
plot(x1, predict(fit, data.frame(X1=x1, X2=6)), type='l', ylim=c(0, 360))
plot(x1, predict(fit, data.frame(X1=x1, X2=9)), type='l', ylim=c(-100, 360))



# Input Steroid Level Data
df <- get("CH11TA08", envir = env)
names(df) <-c("Y", "X")



# TABLE 11.8     (p 453)
# Data Set and 5-Region Regression Tree Fit--Steriod Level Example.
#
# The cut function will easily factor Age into the regions (R) of interest.
# However, since these are left-closed, right-opened, the appropriate
# parameters are specified. The book says to not include Age = 25 (it's right
# open). This doesn't seem right, and cut doesn't allow it. The mean score is
# not too different.
df <- transform(df, 
  R = cut(X, c(8, 9, 10, 13, 14, 25), include.lowest = TRUE, right = FALSE))
with(df, tapply(Y, R, mean))



# FIGURE 11.9     (p 454)
# Fitted Regression Tree, Residual Plot, and
# Regression Tree Diagram--Steroid Level Example.
#
# While the above method demonstrates an approach to categorizing our bins in
# a regression tree, the authors have not demonstrated much in the way of how
# to construct the trees. Therefore, I will simply use R facilities to
# generate the authors output as close as possible. For more on the topic of
# data mining and regression tress I suggest the author's recommendation
# (11.15), Hastie et al. "The Elements of Statistical Learning."
#
# R has a regression tree function rpart (rpart). It takes in a regression
# formula or object and generates a tree. To manage the parameters of the tree,
# you manipulate the control parameter using the rpart function
# rpart.control. The selection of parameters in this case were trial-and-
# error. Since there were at least two Y values in each bin, I set the
# minsplit = 2. The cp value defaults to 0.01, which results in breaking up
# the [14, 25] bin. It turns out that cp = 0.013 maintained the author's
# results. The return object of rpart can be used just like an lm object as
# demonstrated below.
library(rpart)

fit <- rpart(Y ~ X, df, control = rpart.control(minsplit = 2, cp = 0.013))

# FIGURE 11.9a
plot(Y ~ X, df, xlab = "Age", ylab = "Steroid Level")
lines(predict(fit, data.frame(X = sort(X))) ~ sort(X), df, type = 's')  # Step plot

# FIGURE 11.9b 
plot(resid(fit) ~ predict(fit), xlab = "Predicted", ylab = "Residual")

# FIGURE 11.9c
plot(fit); text(fit, use.n = TRUE, xpd = TRUE)



# FIGURE 10.10     (p 455)
# Growing the Regression Tree--Steroid Level Example
#
# To grow the tree is basically to manipulate the parameters on the rpart
# function. Therefore, I will construct a few different models and plot them
# like the above step function. To faciliate that, I will write a convenience
# wrapper that plots the results. All parameters go straight to rpart.
rplot <- function(formula, data, ...)
{
  df <- model.frame(formula, data)
  model <- rpart(Y ~ X, df, ...)
  plot(formula, df, xlab = "AGE", ylab = "Steroid Level")
  lines(predict(model, data.frame(X = sort(X))) ~ sort(X), df, type = 's')
}

par(mfrow = c(2, 2))
rplot(Y ~ X, df, minsplit = 2, cp = 0.095)
rplot(Y ~ X, df, minsplit = 2, cp = 0.04)
rplot(Y ~ X, df, minsplit = 2, cp = 0.03)
rplot(Y ~ X, df, minsplit = 2, cp = 0.015)



# Import University Admissions data (Appendix C.4)
df <- get("APPENC04", envir = env)
names(df) <- c("id", "gpa", "rank", "score", "year")



# FIGURE 11.12     (p 457)
# Regression Tree Results--University Admissions Example
# TO BE COMPLETED



# Input the Toluca Company Data
df <- get("CH11TA09", envir = env)
names(df) <- c("X", "Y")
fit <- lm(Y ~ X, df)



# TABLE 11.9     (p 461)
# Bootstrapping with Fixed X Sampling--Toluca Company Example
#
# Columns (5) and (6) will be different since (5) is selected by a random
# sample from the residuals of the original model. The given trial will not
# be part of the following 1000 bootstraps. This is here for completeness.
estar <- sample(resid(fit), size = nrow(df), replace = TRUE)
cbind(data,
  Yh    = fitted(fit),
  ei    = resid(fit),
  estar = estar,
  ystar = fitted(fit) + estar)



# FIGURE 11.13     (p 461)
# Histogram of Bootstrap Estimates b1*--Toluca Company Example
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
    d1  = b - p[[1]],
    d2  = p[[2]] - b,
    lwr = 2*b - p[[2]],
    upr = 2*b - p[[1]])

  return (list(coefs = coefs, statistics = statistics, confint = confint))
}

z <- bootstrap(fit)
z$statistics; z$confint
hist(z$coefs, breaks = 40, freq = F, main = "", xlab = "Bootstrap b1*")

rm(bootstrap, z)



# Input Blood Pressure data
df <- get("CH11TA10", envir = env)
names(df) <- c("X", "Y")
fit  <- lm(Y ~ X, df)



# TABLE 11.10     (p 462)
# Bootstrapping with Random X Sampling--Blood Pressure Example
xstar <- sample(df$X, size = nrow(df), replace = TRUE)
ystar <- sample(df$Y, size = nrow(df), replace = TRUE)
estar <- abs(resid( lm(ystar ~ xstar) ))
cbind(data,
  'X*' = xstar,
  'Y*' = ystar,
  'e*' = round( resid( lm(ystar ~ xstar) ), 2),
  's*' = round( fitted( lm(estar ~ xstar) ), 2),
  'w*' = round( 1 / (fitted( lm(estar ~ xstar) ))^2, 4)) 



# FIGURE 11.14     (p 463)
# Histogram of Boostrap Estimates b1*--Blood Pressure Example
#
# I first code an algorithm similar to the one used in the previous example.
# I then finish by showing the R approach using the boot function (boot). To
# use boot requires you to define a statistic. In this case, we define a
# function boot.coef that extracts the coefficient of interest as in the
# other defined functions. The boot function then calculates the relevant
# information, and we can use the boot.ci function to generate confidence
# intervals using a number of methods.
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
    d1  = b - p[[1]],
    d2  = p[[2]] - b,
    lwr = 2*b - p[[2]],
    upr = 2*b - p[[1]])

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
z$statistics; z$confint
hist(z$coefs, breaks = 40, freq = F, main = "", xlab = "Bootstrap b1*")

# Using the boot package
z <- boot(data = model.frame(fit), statistic = boot.coef, R = 1000)
print(z); boot.ci(z); 
hist(z$t, breaks = 40, freq = F, main = "", xlab = "Bootstrap b1*")



# Input Traffic Data
#
# R defaults to making the first level (in this case 1) the reference level
# in a linear model. This means that instead of, say, having dummy variables
# for CONTROL = 1 and making CONTROL = 2 absent, it will force CONTROL1 out
# and include CONTROL2. Therefore, we will control for this by looping
# through the qualitative variables and releveling their reference factor to
# the last one (i.e., the size of the levels).
df <- get("CH11TA11", envir = env)
names(df) <- c("AADT", "CTYPOP", "LANES", "WIDTH", 
               "CONTROL", "CLASS", "TRUCK", "LOCALE")
df <- transform(df,
                CONTROL = factor(CONTROL),
                CLASS   = factor(CLASS),
                TRUCK   = factor(TRUCK),
                LOCALE  = factor(LOCALE))
fit <- lm(AADT ~ ., df)



# TABLE 11.11     (p 465)
# Data--MNDOT Traffic Estimation Example
{
  print(head(df))
  cat("...\n")
  print(tail(df))
}



# FIGURE 11.15     (p 466)
# Scatter Plot Matrix--MNDOT Traffic Estimation Example
library(car)
spm(df, diagonal = "none", span = 0.75, reg.line = F, by.groups = TRUE)



# FIGURE 11.16     (p 467)
# All-Possible-Regressions Output--MNDOT Traffic Estimation Example
#
# Since the data was input with the qualitative variables as factor, the
# regsubsets function will automatically split them into dummy variables.
# However, vif function cannot account for the factors giving a different
# result than the book indicates, but it is still indicative nonetheless.
library(leaps)
library(car)

# Some Mentioned Diagnostics
plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual", pch = 19)
vif(fit)
max(cooks.distance(fit), na.rm = TRUE)

# Best Subsets Analysis
best <- regsubsets(formula(fit), df, force.in = c(1, 2), nbest=5, nvmax = 7)
best$xnames <- paste("X", 0:13, sep = "")
with(summary(best), cbind(
    which[,-c(1:3)],
    Cp     = round(cp, 4),
    Rsq    = round(rsq, 4),
    AdjRsq = round(adjr2, 4))) 



# BEST SUBSETS     (p 468)
#
# To collapse CLASS2 and CLASS4 into one group we need only change the levels
# attributed to those factors. I will choose CLASS = 2 as that group. Since
# CLASS was originally releveled so CLASS = 4 came first, we can simply
# rename the first element of its levels as below.
#
# The names for the pool of variables will go as follows:
#   intercept         X0      CTYPOP            X1      LANES           X2
#   CONTROL1          X3      CLASS1            X4      CLASS3          X5
#   CTYPOP2           X6      LANES2            X7      CTYPOP:LANES    X8
#   CTYPOP:CONTROL1   X9      CTYPOP:CLASS1     X10     CTYPOP:CLASS3   X11
#   LANES:CONTROL1    X12     LANES:CLASS1      X13     LANES:CLASS3    X14
#   CONTROL1:CLASS1   X15     CONTROL1:CLASS3   X16
levels(df$CLASS)[2] <- "4"
df <- transform(df, 
  CTYPOP2 = scale(CTYPOP, T, F)^2,
  LANES2  = scale(LANES, T, F)^2)

pool <- AADT ~ (CTYPOP + LANES + CONTROL + CLASS)^2 + CTYPOP2 + LANES2
best <- regsubsets(pool, df, force.in = c(1, 2), nbest = 5)

best$xnames <- paste("X", 0:16, sep = "")
with(summary(best), cbind(
    which[,-c(1:3)],
    Cp     = round(cp, 4),
    Rsq    = round(rsq, 4),
    AdjRsq = round(adjr2, 4)))



# FIGURE 11.17     (p 468)
# Plots of Studentized Residuals versus Fitted
# Values--MNDOT Traffic Estimation Example
par(mfrow = c(1, 2), pch = 19)

fit <- lm(AADT ~ CTYPOP + LANES + CONTROL + CLASS, df)
plot(rstudent(fit) ~ fitted(fit), xlab = "Fitted Values", ylab = "Studentized Residual")
title("(a) First-Order OLS Model")

fit <- lm(AADT ~ CTYPOP + LANES + LANES2 + CONTROL + CTYPOP:CONTROL, df)
plot(rstudent(fit) ~ fitted(fit), xlab = "Fitted Values", ylab = "Studentized Residual")
title("(b) Second-order OLS Model")



# FIGURE 11.18     (p 469)
# Weighted Least Squares Regression Results--MNDOT Traffic Estimation Example
#
# The weights can be extrated using the method used in TABLE 11.1 above. The
# standard deviation function will just be modifed as the authors detail. The
# R values come out as the authors describe. They repeated this four times in
# all. The summaries show that the results (particularly for CONTROL1) do not
# come out the same at all (CONTROL1 is negative).
fit <- lm(AADT ~ CTYPOP + LANES + LANES2 + CONTROL + CTYPOP:CONTROL, df)
summary(fit)  # Compare before iterations
w   <- fitted(lm(abs(resid(fit)) ~ CTYPOP + LANES, df))^(-2)
fit <- update(fit, weights = w)                           # 1st iteration
w   <- fitted(lm(abs(resid(fit)) ~ CTYPOP + LANES, df))^(-2)
fit <- update(fit, weights = w)                           # 2nd iteration
w   <- fitted(lm(abs(resid(fit)) ~ CTYPOP + LANES, df))^(-2)
fit <- update(fit, weights = w)                           # 3nd iteration
w   <- fitted(lm(abs(resid(fit)) ~ CTYPOP + LANES, df))^(-2)
fit <- update(fit, weights = w)                           # 4nd iteration

summary(fit)  # Compare after iterations



# FIGURE 11.19     (p 470)
# Residual Plots for Final Weighted Least Squares
# Regression Fit--MNDOT Traffic Estimation Example
par(mfrow = c(1, 2), pch = 19)
plot(rstudent(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Student. Residual")
title("(a) Residual Plot against Fitted")

qqnorm(rstudent(fit), xlab = "Expected", ylab = "Student. Residual")
qqline(rstudent(fit))
title("(b) Normal Probability Plot")

