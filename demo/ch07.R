# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input the Body Fat Data
df <- get("CH07TA01", envir = env)
names(df) <- c("x1", "x2", "x3", "y")
fit <- lm(y ~ x1 + x2 + x3, df)



# TABLE 7.1     (p 257)
# Basic Data--Body Fat Example
with(df, cbind(
  "Triceps"  = x1,
  "Thigh"    = x2,
  "Midarm"   = x3,
  "Body Fat" = y
));



# TABLE 7.2     (p 257)
# Regression Results for Several Fitted Models--Body Fat Example
#
# One thing to realize with R anova reports is that it doesn't give you an
# SSR for the entire regression. Instead, it gives you the SS for the first
# term in the model, then the extra from the addition of the second term. So,
# with their example of y ~ x1 + x2 the SSR = 385.44; however, the anova
# output provided would need to be summed to obtain that answer. Therefore,
# what we require is:
#
#   SSR(x1) + SSR(x2|x1) = 352.27 + 33.17 = 385.44 = SSR(x1, x2)
#
# Thus, the R anova summary provides the extra SS directly, depending on how
# the terms were entered into the model.
summary(lm(y ~ x1,      df)); anova(lm(y ~ x1,      df))
summary(lm(y ~ x2,      df)); anova(lm(y ~ x2,      df))
summary(lm(y ~ x1 + x2, df)); anova(lm(y ~ x1 + x2, df))
summary(fit);                 anova(fit)



# TABLE 7.4     (p 262)
# ANOVA Table with Decomposition of SSR--Body Fat Example
#
# For numeric subscripting, the following will be utilized
#   Number    Rows                      Columns
#   1         x1        (x1)            Df
#   2         x2        (x2 | x1)       Sum Sq
#   3         x3        (x3 | x1, x2)   Mean Sq
#   4         Residuals (pure error)
fit.aov <- anova(fit)
round(as.table(
  cbind(
  SS = c("SSR(x1, x2, x3)" = sum(fit.aov[1:3, 2]),
         "SSR(x1)"         = fit.aov[1, 2],
         "SSR(x2|x1)"      = fit.aov[2, 2],
         "SSR(x3|x1, x2)"  = fit.aov[3, 2],
         "SSE"             = fit.aov[4, 2],
         "Total"           = sum(fit.aov[, 2])
         ),  # end SS column

  Df = c(                    sum(fit.aov[1:3, 1]),
                             fit.aov[1, 1],
                             fit.aov[2, 1],
                             fit.aov[3, 1],
                             fit.aov[4, 1],
                             sum(fit.aov$Df)
         ),  # end df column

  MS = c(                    sum(fit.aov[1:3, 2]) / sum(fit.aov[1:3, 1]),
                             fit.aov[1, 3],
                             fit.aov[2, 3],
                             fit.aov[3, 3],
                             fit.aov[4, 3],
                             NA
         )  # end MS column
  )  # end cbind
), 2);  # end table and round

rm(fit.aov)



# Uses of Extra Sums of Squares in Tests for Regression Coefficients    (p 264)
# Test Whether A Single Coefficient Can Be Dropped (i.e., bk = 0)
#
# Since the above example demonstrates how to extract each of the components
# in an extra SS analysis, I will assume the reader can figure out how to
# compute the arithmetic with those components. Instead, here will be two
# ways to perform the same analysis. The anova function can compare models in
# precisely this way. R also has the drop1 function that does the same test
# for each of possible "drop 1 term" models.
anova(update(fit, . ~ . - x3), fit)  # the "." indicates "old model terms"
drop1(fit, test = "F")               # Alternative for each term



# Uses of Extra Sums of Squares in Tests for Regression Coefficients    (p 265)
# Test Whether Several Coefficients Can Be Dropped
anova(lm(y ~ x1, df), fit)



# Coefficients of Partial Determination and Correlation     (p 270-1)
sign.1 <- sign(coef(lm(y ~ x1 + x2, df)))[["x1"]]
sign.2 <- sign(coef(lm(y ~ x1 + x2, df)))[["x2"]]
sign.3 <- sign(coef(fit))[["x3"]]

rbind(
  "2|1"  =  c("Rsq" = anova(fit)["x2", 2] / 
                      anova(lm(y ~ x1, df))["Residuals", 2],
              "r"   = sign.2 * sqrt(anova(fit)["x2", 2] /
                      anova(lm(y ~ x1, df))["Residuals", 2])
  ),  # end 2|1
  "3|12" = c(         anova(fit)["x3", 2] / 
                      anova(lm(y ~ x1 + x2, df))["Residuals", 2],
                      sign.3 * sqrt(anova(fit)["x3", 2] /
                      anova(lm(y ~ x1 + x2, df))["Residuals", 2])
  ),  # end 3|12
  "1|2"  = c(         anova(lm(y ~ x2 + x1, df))["x1", 2] /
                      anova(lm(y ~ x2, df))["Residuals", 2],
                      sign.1 * sqrt(anova(lm(y ~ x2 + x1, df))["x1", 2] / 
                      anova(lm(y ~ x2, df))["Residuals", 2])
  )  # end 1|2
);  # end rbind

rm(sign.1, sign.2, sign.3)



# Input the Dwaine Studios Data
df <- get("CH06FI05", envir = env)
names(df) <- c("x1", "x2", "y")



# TABLE 7.5     (pp 276-7)
# Correlation Transformation and Fitted Standardized
# Regression Model--Dwaine Studios Example
with(df, cbind(Sales = y, Population = x1, Income = x2))  # Table 7.5(a)

fit <- lm(y ~ 0 + x1 + x2, data =
  transform(df,   # Standardized Data Transform
    y  = c(scale(y )) / sqrt(nrow(df) - 1),
    x1 = c(scale(x1)) / sqrt(nrow(df) - 1),
    x2 = c(scale(x2)) / sqrt(nrow(df) - 1))
);

model.frame(fit)  # Table 7.5(b) Transformed Data

coef(fit)
with(df,          # obtain regular coefficients
  coef(fit) * c(sd(y) / sd(x1), sd(y) / sd(x2)))  # (7.53)



# Input the Work Crew Productivity Data
df <- get("CH07TA06", envir = env)
names(df) <- c("x1", "x2", "y")



# TABLE 7.6     (p 279)
# Uncorrelated Predictor Variables--Work Crew Productivity Example
cbind("Crew Size"         = df$x1,
      "Bonus Pay"         = df$x2,
      "Crew Productivity" = df$y
);



# TABLE 7.7     (p 280)
# Regression Results when Predictor Variables Are
# Uncorrelated--Work Crew Productivity Example
anova(lm(y ~ .,  df))
anova(lm(y ~ x1, df))
anova(lm(y ~ x2, df))



# Input the Body Fat Data
df <- get("CH07TA01", envir = env)
names(df) <- c("x1", "x2", "x3", "y")



# FIGURE 7.3     (p 284)
# Scatter Plot Matrix and Correlation Matrix
# of the Predictor Variables--Body Fat Example
with(df, pairs(cbind(x1, x2, x3), pch=19))
with(df, cor(cbind(x1, x2, x3)))



# Clean up the R environment from this session
rm(df, fit)

rm(env)  # Ignore this if you're moving on to the next chapter.
