# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input the Toluca Company Data
df        <- get("CH01TA01", envir = env)
names(df) <- c("x", "y")
fit       <- lm(y ~ x, df)



# TABLE 2.1 and FIGURE 2.2     (p 46)
# Results for Toluca Company Example Obtained in Chapter 1
# Regression Output--Toluca Company Example
summary(fit); anova(fit)   

round(
  c(n    = nrow(fit$model),              # 25
    xbar = mean(fit$model[[2]]),         # 70
    b0   = coef(fit)[[1]],               # 62.366
    b1   = coef(fit)[[2]],               # 3.570
    SSE  = anova(fit)["Residuals", 2],   # 54,825.460
    MSE  = anova(fit)["Residuals", 3],   # 2,383.716
    SSx  = sum(scale(x, T, F)^2),        # 19,800
    SSTO = sum(anova(fit)[,2])           # 307,203
  ), 2
); 



# FUNCTION DEFINITION ci
#
# Description
# Predicts values based on a linear model object
#
# Usage
# ci(model, Xh = 0, alpha = 0.1, m = 1,
#    interval = c("confidence", "prediction", "band"))
#
# Arguments
# model    Object of class "lm".
# Xh       An optional vector of prediction values. Default results in
#          in evaluating the interval about the intercept coefficient.
# alpha    Significance level for prediction interval.
# m        Number of new observations on response for a given Xh.
# interval Type of interval calculation.
#
# Notes
# This function mimics the predict.lm function in the stats library. The aim 
# is to follow the calculations provided for the examples in Chapter 2. Since
# the various equations (2.23, 2.25, 2.30, 2.33, 2.36, 2.38a, and 2.39) share
# the same general form, this function could easily generalize to each of the
# provided examples. In particular, the estimated standard deviation equation
# follows the general form in this function where the constant c defaults to
# 0 under confidence intervals, 1 under predictions. The constant m defaults
# to 1 unless otherwise specified to be looking for numerous mean responses.
# Additionally, an interval of "band" is provided. This option uses the
# Working-Hotelling 1-alpha confidence measure from equations (2.40) and
# (2.40a). Everything else is the same as a confidence interval.
#
# The following examples will make use of this pedagogical function alongside
# the use of the R function that provides the same results. It is left to the
# reader to understand the R function arguments provided.
ci <- function(model, Xh = 0, alpha = 0.1, m = 1,
               interval = c("confidence", "prediction", "band")) {

  int <- match.arg(interval)
  c        <- ifelse(int == "prediction", 1, 0)
  x        <- model$model[[2]]
  n        <- length(x)                                   
  df       <- summary(model)$df[2]                        
  SSx      <- sum( scale(x, T, F)^2 )                     
  MSE      <- anova(model)["Residuals", 3]                
  dev      <- Xh - mean(x)
  s        <- sqrt( MSE * (c/m + 1/n + dev^2 / SSx) )     
  
  ifelse(int == "band", 
    val    <- sqrt( 2 * qf(1 - alpha, 2, df) ),
    val    <- qt(1 - alpha/2, df)
  ); ## end ifelse

  pred    <- coef(fit)[[1]] + coef(fit)[[2]] * Xh
  lwr     <- pred - val * s
  upr     <- pred + val * s

  return( cbind(pred, lwr, upr) )
}  # end confidence interval function



# EXAMPLE     (p 49)
# Confidence Interval for beta0
ci(fit)
confint(fit, parm = 1, level = 0.9)



# EXAMPLE     (pp 54-5)
# Confidence Interval for E{Yh}
newdata = c(65, 100)
ci(fit, Xh = newdata)
predict(fit, data.frame(x = newdata), interval = "confidence", level = 0.9)



# EXAMPLE     (p 59)
# Prediction Interval for Yh(new) when Parameters Unknown
newdata = 100
ci(fit, Xh = newdata, interval = "prediction")
predict(fit, data.frame(x = newdata), interval = "prediction", level = 0.9)



# EXAMPLE     (p 61)
# Prediction of Mean of m New Observations Given Xh
ci(fit, newdata, interval = "prediction", m = 3)
predict(fit, data.frame(x = newdata), interval = "prediction", 
        level = 0.9, weights = 3)




# FIGURE 2.6     (p 62)
# Confidence Bands for Regression Line--Toluca Company Example
#
# Since this user-defined function ci behaves similar to predict, it should
# come as no surprise that this approach to regression bands can be performed
# using predict.
#
# The book here desires a Working-Hotelling confidence interval. It results
# in a slightly less confident (i.e., wider) band about the regression line.
# The reason for this has to do with the family confidence interval that is
# discussed in Chapter 4 of the text. The Chapter 4 walk-through will include
# the manual calculations using predict that are required to calculate the
# W statistic. Until that time, the basic approach will be used here.
ci(fit, 100, interval = "band")  # Example case when Xh = 100

plot(y ~ x, df, xlab = "Lot Size X", ylab = "Hours Y")
  ## same as below: curve(predict(fit, data.frame(x = x)), add = TRUE) 
curve(ci(fit, x, interval = "band")[, "pred"], add = TRUE)
curve(ci(fit, x, interval = "band")[,  "lwr"], add = TRUE, col = "red")
curve(ci(fit, x, interval = "band")[,  "upr"], add = TRUE, col = "red")
curve(predict(fit, data.frame(x = x), int = "c")[, 2], add = TRUE, col = "blue")
curve(predict(fit, data.frame(x = x), int = "c")[, 3], add = TRUE, col = "blue")

rm(ci, bands, newdata)



# EXAMPLE     (p 71)
# F-test of beta1
{
  MSR  <- anova(fit)["x",         "Mean Sq"]   # 252,377.600
  MSE  <- anova(fit)["Residuals", "Mean Sq"]   #   2,383.716
  F    <- MSR / MSE                            #     105.876
  fval <- qf(0.95, 1, 23)                      #       4.279
  ifelse(abs(F) <= fval,
    print("Conclude H0"),                      # Fail to reject H0
    print("Conclude Ha")                       # Reject H0
  ); ## end ifelse
  rm(MSR, MSE, F, fval)
}



# EXAMPLE     (pp 75-6)
# Coefficient of Determination and Coefficient of Correlation
round(c(
  SSR  = anova(fit)["x", "Mean Sq"],         # 252,377
  SSTO = sum( anova(fit)["Sum Sq"] ),        # 307,203
  R2   = summary(fit)$r.squared,             # 0.8215 (= SSR / SSTO)
  r    = with(df, cor(x, y))                 # 0.9064 (= sqrt(R2) )
), 4)



# Input the Sales Marketing Data
df <- get("CH02TA04", envir = env)
y1 <- df[, 1]
y2 <- df[, 2]


# TABLE 2.4     (pp 87-9)
# Data on Population, Expenditures and Their Ranks--Sales Marketing Example
cbind(
  Population  = y1,
  Expenditure = y2,
  R1          = rank(y1),
  R2          = rank(y2)
); 

# Calculate Spearman Correlation Information
c(                                                # Spearman Rho Calculation
  rho =  sum( scale(rank(y1), T, F) * scale(rank(y2), T, F) ) /
         sqrt(sum(scale(rank(y1), T, F)^2) * sum(scale(rank(y2), T, F)^2))
); 
cor.test(y1, y2, method = "spearman")$estimate    # Spearman Rho = 0.895
c(rho = cor(rank(y1), rank(y2)))                  # Alternative R method
c(                                                # t* calculation
  tstar = cor(rank(y1), rank(y2)) * sqrt(length(y1) - 2) / 
          sqrt(1 - cor(rank(y1), rank(y2))^2)
);

rm(y1, y2)



# Clean up the R environment from this session
rm(df, fit)

rm(env)  # Ignore this if you're moving on to the next chapter.
