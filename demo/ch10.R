# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input Life Insurance Data
df <- get("CH10TA01", envir = env)
names(df) <- c("X1", "X2", "Y")
fit <- lm(Y ~ ., df)



# TABLE 10.1     (p 387)
# Basic Data--Life Insurance Example
cbind(
  "Avg Income"     = df$X1,
  "Risk Aversion"  = df$X2,
  "Life Insurance" = df$Y
);



# FIGURE 10.3     (p 387)
# Residual Plot and Added-Variable Plot--Life Insurance Example
# Note: The 'car' package contains avPlots(fit) for this, among
# other useful diagnostic functions. See below.
# At this time, Quick-R contains some examples at
# http://www.statmethods.net/stats/rdiagnostics.html
A <- resid(lm(Y ~ X2, df))
B <- resid(lm(X1 ~ X2, df))
par(mfrow = c(1, 2), pch = 19)

plot(resid(fit) ~ X1, df, xlab = "X1", ylab = "Residual")
title("(a) Residual Plot against X1")
abline(0, 0, lty = 2)

plot(A ~ B, xlab = "e(X1|X2)", ylab = "e(Y|X2)")
title("(b) Added-Variable Plot for X1")
abline(lm(A ~ B))
abline(0, 0, lty = 2)

rm(A, B)



# Input the Body Fat Data
df <- get("CH07TA01", envir = env)
names(df) <- c("X1", "X2", "X3", "Y")
fit <- lm(Y ~ X1 + X2, df)



# FIGURE 10.4     (p 389)
# Residual Plots and Added-Variable Plots--Body
# Fat Example with Two Predictor Variables
# Note: The above manual approach can be extended
# to this example, but no more pedagogical benefit
# derives. The reader should be able to match the 
# output to that in the book.
library(car)
avPlots(fit)
residualPlots(fit, fitted = FALSE, quadratic = FALSE, test = FALSE)



# TABLE 10.2     (p 393)
# Illustration of Hat Matrix
#
# See Chapter 5 tutorial on how to do the matrix algebra. This will all be
# done using the R utilities. Note, you can use model.matrix(fit) to obtain
# the X matrix for calculating the hat matrix and error variance accordingly:
#   e.g., X %*% solve(t(X) %*% X) %*% t(X)
fit <- lm(Y ~ X1 + X2, data =
  data.frame(
    Y = c(301, 327, 246, 187),
    X1 = c(14, 19, 12, 11),
    X2 = c(25, 32, 22, 15))
); ## end lm fit

cbind(
  X1   = model.frame(fit)$X1,
  X2   = model.frame(fit)$X2,
  Y    = model.frame(fit)$Y,
  Yhat = fitted(fit), 
  e    = resid(fit),
  h    = hatvalues(fit),
  s    = deviance(fit) * (1 - hatvalues(fit))
); ## end cbind



# FIGURE 10.6     (p 398)
# Illustration of Leverage Values as Distance Measures--Table 10.2 Example
plot(X2 ~ X1, model.frame(fit), pch = 19, xlim = c(10, 20), ylim = c(10, 37),
     xlab = "X1", ylab = "X2", main = "FIGURE 10.6")
text(X2 ~ X1, model.frame(fit), labels = round(hatvalues(fit), 4), pos = 3)
points(mean(X2) ~ mean(X1), model.frame(fit), pch = 22)
text(mean(X2) ~ mean(X1), model.frame(fit), 
     cex = 0.75, pos = 4, labels = "(X1bar, X2bar)")



# TABLE 10.3     (p 397)
# Residuals, Diagonal Elements of the Hat Matrix, and Studentized
# Deleted Residuals--Body Fat Example with Two Predictor Variables
fit <- lm(Y ~ X1 + X2, df)
cbind(
  e = round(resid(fit), 3),
  h = round(hatvalues(fit), 3),
  t = round(rstudent(fit), 3)
);



# FIGURE 10.7     (p 399)
# Scatter Plot of Thigh Circumference against Triceps Skinfold
# Thickness--Body Fat Example with Two Predictor Variables
plot(X2 ~ X1, df, type = 'n', xlim = c(13, 33), ylim = c(40, 60),
     xlab = "Triceps Skinfold Thickness", ylab = "Thigh Circumference")
text(X2 ~ X1, df, labels = seq(nrow(df)), cex=.75)



# TABLE 10.4     (p 402)
# DFFITS, Cook's Distances, and DFBETAS--Body
# Fat Example with Two Predictor Variables
# Note: Compare with the R function influence.measures(fit)
cbind(
  DFFITS  = round(dffits(fit),         4),
  D       = round(cooks.distance(fit), 4),
  DFBETA0 = round(dfbetas(fit)[,1],    4), 
  DFBETA1 = round(dfbetas(fit)[,2],    4), 
  DFBETA2 = round(dfbetas(fit)[,3],    4)
); 



# FIGURE 10.8     (p 404)
# Proportional Influence Plot (Points Proportional
# in Size to Cook's Distance Measure) and Index Influence
# Plot--Body Fat Example with Two Predictor Variables
par(mfrow = c(1, 2), pch = 19)
plot(resid(fit) ~ fitted(fit), cex = cooks.distance(fit)*10,
     xlim = c(10, 30), ylim = c(-4.5, 4.5),
     xlab = "YHAT", ylab = "Residual")
title("(a) Proportional Influence Plot")
plot(seq(nrow(df)), cooks.distance(fit), type = "o", lwd = 2,
     xlab = "Case Index Number", ylab = "Cook's Distance D")
title("(b) Index Influence Plot")



# TABLE 10.5     (p 409)
# Variance Inflation Factors--Body Fat Example with Three Predictor Variables
#
# R has a VIF function in the "car" package that takes in an lm object.
# To get the standardized coefficients we update the model by generically
# scaling the values. The maximum VIF is clear, and an average can easily be
# computed. I leave it to the interested student to perform that simple task.
fit <- lm(scale(Y) ~ scale(X1) + scale(X2) + scale(X3), df)
cbind(
  "Beta*" = round(as.vector(coef(fit)[-1]), 4),
  "VIF"   = round(vif(lm(Y ~ ., df)), 2)
);  # end cbind



# Input the Surgical Unit Data
# Note: R > 2.15 has paste0. No need for sep = "" anymore.
df <- get("CH09TA01", envir = env)
names(df) <- c(paste("X", seq(8), sep = ""), "Y" , "lnY")
fit <- lm(lnY ~ X1 + X2 + X3 + X8, df)
vif(fit)  # For page 412



# FIGURE 10.9     (p 411)
# Residual and Added-Variable Plots for Surgical
# Unit Example--Regression Model (10.45)
# Note: The 'car' package functions used earlier can
# be used to explore the unshown information discussed
# on page 410. This example, however, requires relating
# added-values from variables not within the model. Thus,
# it cannot be used with residualPlots or avPlots. It may
# be possible to fit a full model and use these functions
# with appropriately selected 'terms' parameters, though.
par(mfrow = c(2, 2), pch = 19)

plot(resid(fit) ~ fitted(fit), ylim = c(-0.7, 0.7),
     xlab = "Predicted Value", ylab = "Residual")
title("(a) Residual Plot against Predicted")

plot(resid(fit) ~ df$X5, ylim = c(-0.7, 0.7),
     xlab = "X5", ylab = "Residual")
title("(b) Residual Plot against X5")

plot(resid(fit) ~ resid(update(fit, X5 ~ .)), ylim = c(-0.7, 0.7),
     xlab = "e(X5|X1238)", ylab = "e(Y'|X1238)")
title("(c) Added-Variable Plot for X5")
abline(lm(resid(fit) ~ resid(update(fit, X5 ~ .))) )

qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
title("(d) Normal Probability Plot")
qqline(resid(fit))



# FIGURE 10.10     (p 413)
# Diagnostic Plots for Surgical Unit Example
par(mfrow = c(2, 2), pch = 19)

plot(seq(nrow(df)), rstudent(fit), type = "o",
     main = "(a) Studentized Deleted Residuals",
     xlab = "Case Index", ylab = "t")
plot(seq(nrow(df)), hatvalues(fit), type = "o",
     main = "(b) Leverage Values",
     xlab = "Case Index", ylab = "h")
plot(seq(nrow(df)), cooks.distance(fit), type = "o",
     main = "(c) Cook's Distance",
     xlab = "Case Index", ylab = "D")
plot(seq(nrow(df)), dffits(fit), type = "o",
     main = "(d) DFFITS values",
     xlab = "Case Index", ylab = "DFFITS")



# TABLE 10.6     (p 413)
# Various Diagnostics for Outlying Cases--Surgical Unit Example
cases <- c(17, 23, 28, 32, 38, 42, 52)
cbind(
  e      = round(resid(fit), 4),
  t      = round(rstudent(fit), 4),
  h      = round(hatvalues(fit), 4),
  D      = round(cooks.distance(fit), 4),
  DFFITS = round(dffits(fit), 4)
)[cases, ]

rm(cases)



# Clean up the R environment from this session
rm(df, fit)

rm(env)  # Ignore this if you're moving on to the next chapter.
