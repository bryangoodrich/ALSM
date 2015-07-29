# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input the Power Cells Data
df <- get("CH08TA01", envir = env)
names(df) <- c("Y", "X1", "X2")
df <- transform(df, 
  x1   = scale(X1, scale = 0.4),
  x2   = scale(X2, scale = 10),
  x1sq = scale(X1, scale = 0.4)^2,
  x2sq = scale(X2, scale = 10)^2,
  x1x2 = scale(X1,, 0.4) * scale(X2,, 10)
);  # end data transformation
fit <- lm(Y ~ ., df)



# Table 8.1     (p 300)
# Data--Power Cells Example
model.frame(df)



# FIGURE 8.5     (p 303)
# Diagnostic Residual Plots--Power Cells Example
par(mfrow = c(2, 2), pch = 19)
plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual")
title("(a) Residual Plot Against Fitted Values")

plot(resid(fit) ~ x1, df, xlab = "x1", ylab = "Residual")
title("(b) Residual Plot against x1")

plot(resid(fit) ~ x2, df, xlab = "x2", ylab = "Residual")
title("(c) Residual Plot against x2")

qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
title("(d) Normal Probability Plot")
qqline(resid(fit))



# FIGURE 8.4     (p 301-5)
# Regression Output for Second-Order Polynomial Model--Power Cells Example
#
# A lot of analysis was performed over these pages, both shown and not shown.
# For that reason I will follow a similar process. However, results that can
# be easily obtained by summary functions will be ignored. The F-test that is
# performed is not that in the summary results (either here or in the SAS
# output the authors provide). The summary results show the F-test as to
# whether or not all the predictor variables are relevant. This statistic is
# 10.57, which is still less than the 19.2. Recall from Chapter 3 that the
# alr3 library contains a function for an expanded anova table.
library(alr3)
fit <- lm(Y ~ ., df)                                 # Restate the model
summary(fit)
pureErrorAnova(fit)                                  # alr3 function
anova(lm(Y ~ 1, df), fit)                            # Explicit F-Test

# Calculate F*     (p 302)
pureErrorAnova(fit)[" Lack of fit", "Mean Sq"] /     # (6.68b)
pureErrorAnova(fit)[" Pure Error",  "Mean Sq"]       # 1.820

with(df, rbind(                                      # Correlation (p 301)
  "X and Xsq" = round(c("1" = cor(X1, X1^2), "2" = cor(X2, X2^2)), 3),
  "x and xsq" = round(c("1" = cor(x1, x1^2), "2" = cor(x2, x2^2)), 3))
);

anova(lm(Y ~ x1 + x2, df), fit)                      # partial F-test
fit <- lm(Y ~ x1 + x2, df)                           # refit model (p 304)

# Diagnostic Plots
par(mfrow = c(2, 2), pch = 19)                       # Same as FIGURE 8.5
plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual")
title("(a) Residual Plot Against Fitted Values")

plot(resid(fit) ~ x1, df, xlab = "x1", ylab = "Residual")
title("(b) Residual Plot against x1")

plot(resid(fit) ~ x2, df, xlab = "x2", ylab = "Residual")
title("(c) Residual Plot against x2")

qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
title("(d) Normal Probability Plot")
qqline(resid(fit))
                                                     # Bonferroni procedure on
confint(lm(Y ~ X1 + X2, df))                         # original coefs. (p 305)



# FIGURE 8.6     (p 304)
# Plot of Fitted Response Plane--Power Cells Example
library(Rcmdr)
with(df,
  scatter3d(X2, Y, X1, xlab = "Chage Rate", ylab = "Cycles", zlab = "Temp"))



# Input the Body Fat Data
df <- get("CH07TA01", envir = env)
names(df) <- c("X1", "X2", "X3", "Y")
df <- transform(df,
  x1 = scale(X1, T, F),
  x2 = scale(X2, T, F), 
  x3 = scale(X3, T, F)
);  # end transformation
fit <- lm(Y ~ (X1 + X2 + X3)^2, df)



# Implementation of Interaction Regression Models                   (p 312-3)
# Body Fat Example
cor(model.matrix(fit)[, -1])         # Interactions, high correlation (p 312)

fit <- lm(Y ~ (x1 + x2 + x3)^2, df)  # fit centered values instead
cor(model.matrix(fit)[, -1])         # less correlation, but still present
summary(fit); anova(fit)
anova(lm(Y ~ x1 + x2 + x3, df), fit) # partial F-test: F = 0.53 < qf(.95, 3,13)



# Input the Insurance Innovation Data
df <- get("CH08TA02", envir = env)
names(df) <- c("Y", "X1", "X2")
df <- transform(df, X2 = factor(X2, labels = c("Mutual", "Stock")))
fit <- lm(Y ~ X1 + X2, df)



# FIGURE 8.11     (p 316)
# Plot of Y against X1 at both levels of indicator
# variable--Insurance Innovation Example
plot(Y ~ X1, df, pch = 19, ylim = c(0, 40), xlim = c(0, 310),
     xlab = "Size of Firm", ylab = "Months Elapsed")
abline(coef(fit)[[1]],                   # Intercept: Mutual Firm
       coef(fit)[[2]],)                  # Slope:     Mutual Firm
abline(coef(fit)[[1]] + coef(fit)[[3]],  # Intercept: Stock  Firm
       coef(fit)[[2]], col = "red")      # Slope:     Stock  Firm



# TABLE 8.2     (p 317)
# Data and Indicator Coding--Insurance Innovation Example
data.frame(                                  # Note: factor levels while cbind
    "Months"      = df$Y,                    # will show R storage values of
    "Firm Size"   = df$X1,                   # 1 and 2. However, matrix.model
    "Firm Type"   = df$X2,                   # shows R fits factors correctly.
    "Indicator"   = model.matrix(fit)[, 3],  
    "Interaction" = df$X1 * (unclass(df$X2) - 1)  
);



# TABLE 8.3     (p 317)
# Regression Results for Fit of Regression
# Model (8.33)--Insurance Innovation Example
#
# Also included are the results for the confidence intervals of the
# parameters. The formal test is inherent in the summary as the p-values and
# printed stars demonstrate these coeffficients are significant at an alpha
# of 0.01. This is confirmed by the respective confidence interval not
# including zero.
summary(fit); anova(fit); confint(fit)



# FIGURE 8.12     (p 318)
# Fitted Regression Function for Regression
# Model (8.33)--Insurance Innovation Example
fit <- lm(Y ~ X1 + X2, df)
plot(Y ~ X1, df, col = X2, pch = 19, ylim = c(0, 40), xlim = c(0, 310),
     xlab = "Size of Firm", ylab = "Months Elapsed")
abline(coef(fit)[[1]],                   # Intercept: Mutual Firm
       coef(fit)[[2]])                   # Slope:     Mutual Firm
abline(coef(fit)[[1]] + coef(fit)[[3]],  # Intercept: Stock  Firm
       coef(fit)[[2]], col = "red")      # Slope:     Stock  Firm



# FIGURE 8.14 and FIGURE 8.15     (p 325)
# Plot of Y against X1 at both levels of indicator
# variable--Insurance Innovation Example
#
# Note that the author does an illustration that exaggerates the change in
# in slope. As is evidenced from this example, the slope is relatively
# unchanged so that the two lines appear parallel. The difference between the
# two is the coefficient of X1X2 = -0.000417
#
# Here I set a new model that nests the X2 levels. By excluding the intercept
# term, the coefficients 1 and 3 apply to the first X2 level and the 
# coefficients 2 and 4 apply to the second X2 level.
fit <- lm(Y ~ X2/X1 - 1, df)  # Using nested model for ease of presentation
plot(Y ~ X1, df, col = X2, pch = 19, ylim = c(0, 40), xlim = c(0, 310),
     xlab = "Size of Firm", ylab = "Months Elapsed")
abline(coef(fit)[[1]],               # Intercept: Mutual Firm
       coef(fit)[[3]])               # Slope:     Mutual Firm
abline(coef(fit)[[2]],               # Intercept: Stock  Firm
       coef(fit)[[4]], col = "red")  # Slope:     Stock  Firm



# TABLE 8.4     (p 327)
# Regression Results for Fit of Regression Model (8.49)
# with Interaction Term--Insurance Innovation Example
fit <- lm(Y ~ X1 * X2, df)  # Going back to intended model
summary(fit); anova(fit); anova(lm(Y ~ X1 + X2, df), fit)



# Input the Soap Production Lines Data
df <- get("CH08TA05", envir = env)
names(df) <- c("Y", "X1", "X2")
df <- transform(df, X2 = factor(X2, labels = c("Line1", "Line2")))



# TABLE 8.5     (p 330)
# Data--Soap Production Lines Example
data.frame(
  "Scrap"      = df$Y,
  "Line_Speed" = df$X1,
  "Indicator"  = df$X2
);



# FIGURE 8.16     (p 331)
# Symbolic Scatter Plot--Soap Production Lines Example
fit <- lm(Y ~ X2/X1 - 1, df)  # Using nested model for ease of presentation
plot(Y ~ X1, df, col = X2, ylim = c(100, 500), xlim = c(100, 350), pch = 19,
     xlab = "Line Speed", ylab = "Amount of Scrap")
abline(coef(fit)[[1]],               # Intercept: Production Line 1
       coef(fit)[[3]])               # Slope:     Production Line 1
abline(coef(fit)[[2]],               # Intercept: Production Line 2
       coef(fit)[[4]], col = "red")  # Slope:     Production Line 2



# FIGURE 8.17     (p 332)
# Residual Plots against Yhat--Soap Production Lines Example
#
# The author mentions a residual plot against X2 and a normal probability
# plot (page 331), but did not include them. I have included them below.
# Note that 'fit' is still referring to the fitted model, but the results
# are the same as lm(Y ~ X1*X2, df). 
plot(fitted(fit), resid(fit), col = df$X2, pch = 19,
     xlim = c(100, 500), ylim = c(-40,40),
     ylab = "Residual", xlab = "Fitted")
title("(a) Production Line 1")
abline(0, 0)

plot(jitter(as.numeric(df$X2), 0.1), resid(fit), pch = 19, 
     xlab = "Production Line", ylab = "Residual", sub = "jittered added to X2")

qqnorm(resid(fit), xlab = "Expected", ylab = "Residual", main = "")
title("Normal Probability Plot")
qqline(resid(fit))



# TABLE 8.6     (p 332)
# Regression Results for Fit of Regression
# Model (8.55)--Soap Production Lines Example
fit <- lm(Y ~ X1 * X2, df)  # No longer using the nested model
summary(fit); anova(fit)
   


# Clean up the R environment from this session
rm(df, fit)

rm(env)  # Ignore this if you're moving on to the next chapter.
