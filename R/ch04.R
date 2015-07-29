# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input the Toluca Company Data
df  <- get("CH01TA01", envir = env)
names(df) <- c('x', 'y')
fit  <- lm(y ~ x, df)



# Bonferroni Joint Confidence Intervals     (p 156)
#
# As the authors explain, we're looking for B = t(1 - 0.10/4; 23) for the 90%
# joint confidence interval. But notice that the confint funtion uses a 95%
# confidence interval t(1 - 0.05/2; 23) = B! Therefore, all we require in
# this instance is to use confint(fit) to get the correct result.
#
# The reason this works is because the Bonferroni correction is just a
# modification to the t-statistic. The confint function uses the level
# parameter to produce the alpha = 0.05 in the above case. Therefore, what we
# require is a level to turn 1 - 0.05/2 into 1 - 0.1/4, or a level such that:
#   1 - a/2 = 1 - (1 - level)/2 = 1 - alpha/4,   alpha = 0.1, a = (1 - level)
#
# The 4, as we shall see later, is really just 2g, where g is the number of
# confidence coefficients in the family, in this case 2. It turns out that a
# level = 1 - alpha/g = 1 - 0.1/2 will give us the desired result:
#   1 - a/2 = 1 - (1 - level)/2 = 1 - (1 - (1 - alpha/g))/2 = 1 - alpha/2g
#
# In this case, 1 - alpha/g = 0.95, the default level of confint.
confint(fit, level = (1 - 0.1/2))  # = confint(fit) since g = 2



# Simultaneous Estimation of Mean Responses     (p 158)
# The Working-Hotelling Procedure (90% confidence coefficient)
ci.wh <- function(model, newdata, alpha = 0.1) {
  df    <- nrow(model.frame(model)) - length(coef(model))  # 23
  W     <- sqrt( 2 * qf(1 - alpha, 2, df) )                # 2.2580
  ci    <- predict(model, newdata, se.fit = TRUE)   
  x <- cbind(
    x   = newdata,
    s   = ci$se.fit,
    fit = ci$fit,
    lwr = ci$fit - W * ci$se.fit,
    upr = ci$fit + W * ci$se.fit
  )  # end cbind
  
  return(x)
} # end wh function

new <- data.frame(x = c(30, 65, 100))
ci.wh(fit, new)



# Simultaneous Estimation of Mean Responses     (p 159)
# The Bonferroni Procedure (90% confidence coefficient)
#
# As in the above Bonferroni correction, we require level = (1 - alpha/g).
# For clarity, part of the return object includes se.fit (set to true), which
# are the standard errors, and residual.scale, which is equal to the positive
# squart root of the MSE. This will be used in calculations below.
predict(fit, new, int = "c", level = (1 - 0.1/nrow(new)), se.fit = TRUE)



# Simultaneous Prediction Intervals for New Observations     (p 160)
# The Scheffe and Bonferroni Procedures (95% confidence coefficient)
ci.sim <- function(model, newdata, type = c("B", "S"), alpha = 0.05) {
  g  <- nrow(newdata)
  CI <- predict(model, newdata, se.fit = TRUE)
  M  <- ifelse(match.arg(type) == "B",
          qt(1 - alpha / (2*g), model$df),              # B = (4.9a)
          sqrt( g * qf( 1 - alpha, g, model$df) )       # S = (4.8a)
        ); ## End ifelse type
  spred <- sqrt( CI$residual.scale^2 + (CI$se.fit)^2 )  # (2.38) 
  x <- data.frame(
    "x"     = newdata,
    "spred" = spred,
    "fit"   = CI$fit,
    "lower" = CI$fit - M * spred,
    "upper" = CI$fit + M * spred
  )  # end data.frame
  
  return(x)
} # end sim function

new <- data.frame(x = c(80, 100))
ci.sim(fit, new, type = "S")
ci.sim(fit, new, type = "B")

rm(ci.sim, ci.wh, new)



# Input the Warehouse data
df <- get("CH04TA02", envir = env)
names(df) <- c('x', 'y')



# TABLE 4.2 and FIGURE 4.1     (pp 162-4)
# Scatter Plot and Fitted Regression through Origin--Warehouse Example
#
# Included in the summary and anova output are the values calculated on
# pages 163 and 164.
fit <- lm(y ~ 0 + x, df)   # alternate way to exclude intercept: lm(y ~ x - 1)
tab <- transform(df, 
  "x"   = round(x, 0),
  "y"   = round(y, 0),
  "xy"  = round(x * y, 0),
  "xsq" = round(x * x, 0),
  "fit" = round(fitted(fit), 2),
  "e"   = round(resid(fit), 2)
); 
rbind(tab, Total = colSums(tab))

plot(df$x, fitted(fit), pch = 19, 
     xlab = "Work United Performed", ylab = "Variable Labor Costs")
abline(fit)

# Some of the results from pp 163-4 are included or derivable from these
summary(fit); anova(fit); confint(fit)



# Clean up the R environment from this session
rm(df, fit, tab)

rm(env)  # Ignore this if you're moving on to the next chapter.

