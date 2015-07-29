# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input the Surgical Unit Data
df <- get("CH09TA01", envir = env)
names(df) <- c(paste("X", 1:8, sep = ""), "Y" , "lnY")
data.names <- c("Bloodclot", "Progindex", "Enzyme", "Liver",
                "Age", "Gender", "Alc.Mod", "Alc.Heavy",
                "Survival", "LnSurvival")



# TABLE 9.1     (p 350)
# Potential Predictor Variables and Response Variable--Surgical Unit Example
head(df); cat("... \n"); tail(df)



# Diagnostic Results Not Shown in Book     (p 351)
#
# Included here are the stem-and-leaf plots for each of the predictor
# variables chosen to be used at this time (i.e., the first four). Also
# included is the scatterplot matrix and correlation matrix.
#
# Since I think the boxplot is more visually telling than the stem-and-leaf
# plots, I have also included a 2x2 graphics output of the boxplots for each
# of the predictor variables.
with(df, {
  stem(X1, 2)
  stem(X2, 4)
  stem(X3, 4)
  stem(X4)
  cor(df[1:4])
  }   
);  # end with
pairs(df[1:4], labels = data.names[1:4])

with(df, {
  par(mfrow = c(2, 2))
  boxplot(X1, main = "Bloodclot")
  boxplot(X2, main = "Progindex")
  boxplot(X3, main = "Enzyme"   )
  boxplot(X4, main = "Liver"    )
  }
);  # end with



# FIGURE 9.2     (p 351)
# Some Preliminary Residual Plots--Surgical Unit Example
par(mfrow = c(2, 2))

fit <- lm(Y ~ X1 + X2 + X3 + X4, df)
plot(fitted(fit), resid(fit), xlab = "Predicted value", ylab = "Residual")
title("(a) Residual Plot for Y")
abline(0, 0)

qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
title("(b) Normal Plot for Y")
qqline(resid(fit))

fit <- lm(lnY ~ X1 + X2 + X3 + X4, df)
plot(fitted(fit), resid(fit), xlab = "Predicted value", ylab = "Residual")
title("(c) Residual Plot for lnY")
abline(0, 0)

qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", main = "")
title("(d) Normal Prot for lnY")
qqline(resid(fit))



# FIGURE 9.3     (p 352)
# Scatter Plot Matrix and Correlation Matrix when
# Response Variable is lnY--Surgical Unit Example
cor(model.frame(fit))
pairs(model.frame(fit), labels = data.names[c(10, 1:4)], pch = 18, col = "grey50")



# TABLE 9.2     (p 353)
# Selected Criteria Values for Model Selection for All
# Possible Regression Models--Surgical Unit Example
#
# R contains many ways by which to access the relevant selection criteria
# values. There is also the regsubsets function in the "leaps" library that
# calculates most of these values for all the subsets you specify. However,
# if you desired a way to obtain the values for each combination of models,
# then you could use an algorithm similar to the follow:
#
#   models <- list()
#   for(i in 1:p) models[[i]] <- t(combn(1:p, i))
#   models <- lapply(models, function(set)
#     apply(set, 1, function(row) lm(y ~ ., data[row])))
#   unlist(lapply(models, function(set) lapply(set, sumfun)))
#
# The requirements are that p is the number of parameters to combine, the
# data needs to be arranged so the first columns are the predictors and "y"
# names the response variable (use transform on the full model's frame), and
# "sumfun" names either your own or some defined function to operate on lm
# objects. In this case, there are functions for AIC (stats), PRESS (MPV),
# and others. Since the "models" object contains a list of lists of the
# possible combinations of p predictors, the use of "unlist" is to return the
# raw values as a vector. Then you can cbind them into a new dataframe to
# summarize results. Note, for BIC set the AIC k value to log(n).
#
# Since regsubsets provides most of the criteria in TABLE 9.2, I will simply
# make use of that output. I will also provide some plotting with regsubsets.
# Since this procedure will be repeated a few times, I will make a "wrapper"
# function to generate the table of interest.
library(leaps)

best <- function(model, ...) {
  subsets <- regsubsets(formula(model), model.frame(model), ...)
  subsets <- with(summary(subsets),
    cbind(p = as.numeric(rownames(which)), which, rss, rsq, adjr2, cp, bic)
  );  # end with summary subsets. 
  
  return(subsets)
}  # end best subsets function

round(best(fit, nbest = 4), 4)

subsets <- regsubsets(formula(fit), model.frame(fit), nbest = 4)

par(mfrow = c(2, 2))
plot(subsets)                   # Shaded means included in model.
plot(subsets, scale = "Cp")     # See help page for 
plot(subsets, scale = "adjr2")  # plot.regsubsets for details.
plot(subsets, scale = "r2")     # Useful when table is large.

rm(subsets)



# FIGURE 9.4     (p 356)
# Plot of Variables Selection Criteria--Surgical Unit Example
# The 1:4 is hard-coded into the 'lines' calls. This could be
# automated with something like this: seq(unique(x[, 'p']))
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



# FIGURE 9.5     (p 362)
# Plot of Variable Selection Criteria with All
# Eight Predictors--Surgical Unit Example
fit <- lm(lnY ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8, df)
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

rm(x)



# FIGURE 9.6     (p 363)
# Output for "Best" Two Subsets for Each
# Subset Size--Surgical Unit Example
round(best(fit, nbest = 2), 4) # Fundamentally same info as in TABLE 9.3
rm(best)



# FIGURE 9.7     (p 366)
# Forward Stepwise Regression Output--Surgical Unit Example
#
# R comes with a number of stepwise functions. The most basic of these 
# is "step" contained in the base packages of R. The MASS package 
# has "stepAIC" which provides for a wider range of object classes. 
# The Rcmdr package has the "stepwise" function that appears to pose an
# easier interface to using, at least, the BIC instead of AIC. As below,
# you need to define the 'k' parameter appropriately to use the BIC in
# your stepwise regression methods. Even more, libraries like the wle
# package has the mle.stepwise function that will be ignored here. 
#
# If one wanted to run the process themselves, there exists the intuitive
# functions add1 and drop1 that provide for an easy way to accomplish that.
# They provide RSS and AIC to make comparisons at each step. For repeated
# use, the 'x' parameter can be passed a model matrix, but no checks are
# done on its validity. 
BIC <- log(nrow(model.frame(fit)))
fit <- lm(lnY ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8, df)  # Full Model
step(fit, direction = "forward")      # AIC = -160.77 Drop None
step(fit, dir = "forward",  k = BIC)  # BIC = -142.87 Drop None
step(fit, dir = "backward")           # AIC = -163.83 Drop X4 and X7 
step(fit, dir = "backward", k = BIC)  # BIC = -153.41 Drop X4 through X7
                                      # This model matches the book's search
rm(BIC)



# Input the Surgical Unit Validation Data
df <- do.call("rbind.data.frame", list(
  get("CH09TA01", envir = env),  # Training Data   (n = 54)
  get("CH09TA05", envir = env))  # Validation Data (n = 54)
);  # end do.call
names(df) <- c(paste("X", seq(8), sep = ""), "Y", "lnY")
df <- transform(df, class = gl(2, 54, labels = c("training", "validation")))



# TABLE 9.5     (p 374)
# Potential Predictor Variables and
# Response Variable--Surgical Unit Example
head(df); cat("... \n"); tail(df)



# TABLE 9.4     (p 373)
# Regression Results for Candidate Models (9.21), (9.22), and (9.23)
# Based on Model-Building and Validation Data Sets--Surgical Unit Example
#
# Instead of a table this process will create a list object containing the
# results provided for each of the model fits for each of the data sets.
library(MPV)  # For PRESS wrapper function
newsummary <- function(model) {
  list(
    coefs = round(t(summary(model)$coef[, 1:2]), 4),
    criteria = cbind(
      SSE   = anova(model)["Residuals", "Sum Sq"],
      PRESS = PRESS(model),
      MSE   = anova(model)["Residuals", "Mean Sq"],
      Rsq   = summary(model)$adj.r.squared
    )  # end cbind
  )  # end return list
}  # end model summary function

newsummary(lm(lnY ~ X1 + X2 + X3 + X8,           df, class == "training"))
newsummary(lm(lnY ~ X1 + X2 + X3 + X8,           df, class == "validation"))

newsummary(lm(lnY ~ X1 + X2 + X3 + X6 + X8,      df, class == "training"))
newsummary(lm(lnY ~ X1 + X2 + X3 + X6 + X8,      df, class == "validation"))

newsummary(lm(lnY ~ X1 + X2 + X3 + X5 + X6 + X8, df, class == "training"))
newsummary(lm(lnY ~ X1 + X2 + X3 + X5 + X6 + X8, df, class == "validation"))

rm(newsummary)



# Clean up the R environment from this session
rm(df, fit, data.names)

rm(env)  # Ignore this if you're moving on to the next chapter.
