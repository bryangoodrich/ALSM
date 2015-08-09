## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
data("CH07TA01", package = "ALSM")
data("CH06FI05", package = "ALSM")
data("CH07TA06", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH07TA01
names(.data) <- c("x1", "x2", "x3", "y")
fit <- lm(y ~ x1 + x2 + x3, .data)

## ------------------------------------------------------------------------
with(.data, cbind("Triceps"  = x1,
               "Thigh"    = x2,
               "Midarm"   = x3,
               "Body Fat" = y))

## ------------------------------------------------------------------------
summary(lm(y ~ x1, .data))
anova(lm(y ~ x1, .data))
summary(lm(y ~ x2, .data))
anova(lm(y ~ x2, .data))
summary(lm(y ~ x1 + x2, .data))
anova(lm(y ~ x1 + x2, .data))
summary(fit)
anova(fit)

## ------------------------------------------------------------------------
fit.aov <- anova(fit)
tab <- as.table(cbind(
  'SS' = c("SSR(x1, x2, x3)" = sum(fit.aov[1:3, 2]),
         "SSR(x1)"           = fit.aov[1, 2],
         "SSR(x2|x1)"        = fit.aov[2, 2],
         "SSR(x3|x1, x2)"    = fit.aov[3, 2],
         "SSE"               = fit.aov[4, 2],
         "Total"             = sum(fit.aov[, 2])),
  
  'Df' = c(                    sum(fit.aov[1:3, 1]),
                               fit.aov[1, 1],
                               fit.aov[2, 1],
                               fit.aov[3, 1],
                               fit.aov[4, 1],
                               sum(fit.aov$Df)),
  
  'MS' = c(                    sum(fit.aov[1:3, 2]) / sum(fit.aov[1:3, 1]),
                               fit.aov[1, 3],
                               fit.aov[2, 3],
                               fit.aov[3, 3],
                               fit.aov[4, 3],
                               NA)
))

round(tab, 2)

## ------------------------------------------------------------------------
anova(update(fit, . ~ . - x3), fit)  # the "." indicates "old model terms"
drop1(fit, test = "F")               # Alternative for each term

## ------------------------------------------------------------------------
anova(lm(y ~ x1, .data), fit)

## ------------------------------------------------------------------------
sign.1 <- sign(coef(lm(y ~ x1 + x2, .data)))[["x1"]]
sign.2 <- sign(coef(lm(y ~ x1 + x2, .data)))[["x2"]]
sign.3 <- sign(coef(fit))[["x3"]]

rbind(
  "2|1"  =  c("Rsq" = anova(fit)["x2", 2] / anova(lm(y ~ x1, .data))["Residuals", 2],
              "r"   = sign.2 * sqrt(anova(fit)["x2", 2] / anova(lm(y ~ x1, .data))["Residuals", 2])),
  
  "3|12" = c(         anova(fit)["x3", 2] / anova(lm(y ~ x1 + x2, .data))["Residuals", 2],
                      sign.3 * sqrt(anova(fit)["x3", 2] / anova(lm(y ~ x1 + x2, .data))["Residuals", 2])),
  
  "1|2"  = c(         anova(lm(y ~ x2 + x1, .data))["x1", 2] / anova(lm(y ~ x2, .data))["Residuals", 2],
                      sign.1 * sqrt(anova(lm(y ~ x2 + x1, .data))["x1", 2] / anova(lm(y ~ x2, .data))["Residuals", 2]))
)

## ------------------------------------------------------------------------
.data <- CH06FI05
names(.data) <- c("x1", "x2", "y")

## ------------------------------------------------------------------------
with(.data, cbind(Sales = y, Population = x1, Income = x2))  # Table 7.5(a)

# Data used in this model is a standardized transform, performed inline
fit <- lm(y ~ 0 + x1 + x2, data =
  transform(.data,   
    y  = c(scale(y )) / sqrt(nrow(.data) - 1),
    x1 = c(scale(x1)) / sqrt(nrow(.data) - 1),
    x2 = c(scale(x2)) / sqrt(nrow(.data) - 1)))

model.frame(fit)  # Table 7.5(b) Transformed Data

coef(fit)

# Return to regular coefficients
with(.data, coef(fit) * c(sd(y) / sd(x1), sd(y) / sd(x2)))  # (7.53)

## ------------------------------------------------------------------------
.data <- CH07TA06
names(.data) <- c("x1", "x2", "y")

## ------------------------------------------------------------------------
cbind("Crew Size"         = .data$x1,
      "Bonus Pay"         = .data$x2,
      "Crew Productivity" = .data$y)

## ------------------------------------------------------------------------
anova(lm(y ~ .,  .data))
anova(lm(y ~ x1, .data))
anova(lm(y ~ x2, .data))

## ------------------------------------------------------------------------
.data <- CH07TA01
names(.data) <- c("x1", "x2", "x3", "y")

## ------------------------------------------------------------------------
with(.data, pairs(cbind(x1, x2, x3), pch=19))
with(.data, cor(cbind(x1, x2, x3)))

