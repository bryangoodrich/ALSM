## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
data("CH16TA01", package = "ALSM")
data("CH17TA02", package = "ALSM")
data("CH17TA04", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH17TA02
names(.data) <- c("y", "x1", "x2")
fit <- lm(y ~ factor(x1), .data)

## ------------------------------------------------------------------------
addmargins(xtabs(y ~ x2 + x1, .data), 1, mean)
cat("Ybar.. = ", mean(.data$y))
anova(fit)
anova(fit)[1, 3] / anova(fit)[2,3]  # F* = MSTR / MSE = 866.12
qf(0.95, 3, 36)                     # F(0.95, 3, 36) = 2.87 < F*

## ------------------------------------------------------------------------
.data <- CH16TA01
names(.data) <- c("y", "x1", "x2")

## ------------------------------------------------------------------------
with(.data, stripchart(tapply(y, x1, mean), pch = 20, ylim = c(1, 2), xlim = c(12, 30), xlab = "Cases Sold"))
with(.data, text(tapply(y, x1, mean), rep(1, 4), labels = seq(4), pos = 3))

## ----plots, fig.width=12, fig.height=6-----------------------------------
# The factor means to be plotted
means <- with(.data, tapply(y, x1, mean))

# The radius portion of the 95% confidence interval for each factor
diff  <- with(.data, qt(0.975, 15) * sqrt(10.55 / tapply(y, x1, length)))

par(mfrow = c(1, 2))
barplot(means, xlab = "Design", ylab = "Cases Sold")
title("(a) Bar Graph")

plot(means, type = "o", pch = 19, ylim = c(0, 30), xlab = "Design", ylab = "Cases Sold")
abline(h = mean(.data$y), lty = 2, col = "gray40", lwd = 2)
title(main = "(b) Main Effects Plot")

bar <- barplot(means, xlab = "Design", ylab = "Cases Sold", ylim = c(0, 32))
arrows(bar, means+diff, bar, means-diff, angle = 90, code = 3)
title("(a) Bar-Interval Graph")

plot(means, pch = 19, ylim = c(0, 30), xlim = c(.5, 4.5), xlab = "Design", ylab = "Cases Sold")
arrows(seq(4), means+diff, seq(4), means-diff, angle = 90, code = 3)
abline(h = mean(.data$y), lty = 2, col = "gray40", lwd = 2)
title("(b) Interval Plot")

## ------------------------------------------------------------------------
means <- with(.data, tapply(y, x1, mean))         # Defined from earlier
D = means['3'] - means['4']                    # (17.11)
se = sqrt(10.55 * (1/4 + 1/5))                 # (17.14)
c("lwr" = D, "upr" = D) +                      # (17.16) +/- "t(0.975, 15) x se"
  c(-qt(0.975,15) * se, qt(0.975, 15) * se)
cat('t* =', D / se)                            # (17.18)

## ------------------------------------------------------------------------
lengths <- with(.data, tapply(y, x1, length))
means   <- with(.data, tapply(y, x1, mean))
cont    <- c(1/2, 1/2, -1/2, -1/2)              # contrasts
L  <- sum(cont * means)                         # (17.20)
se <- sqrt(10.55 * sum(cont^2 / lengths))       # (17.22)
c('lwr' = L, 'upr' = L) +                       # (17.24)
  c(-qt(.975, 15)*se, qt(.975, 15)*se)
cat('t* =', L / se)                             # (17.26)

## ------------------------------------------------------------------------
.data <- CH17TA02
names(.data) <- c("y", "x1", "x2")

## ------------------------------------------------------------------------
fit <- aov(y ~ factor(x1), .data)
TukeyHSD(fit)
plot(TukeyHSD(fit))  # Diff 4-1 contains 0

## ------------------------------------------------------------------------
.data <- CH16TA01
names(.data) <- c("y", "x1", "x2")

## ------------------------------------------------------------------------
Scheffe <- function(formula, data, cont, conf.level = 0.90, MSE) 
{
    f  <- function(contrast, mean)   {contrast %*% mean}                # (17.41)
    h  <- function(contrast, n, MSE) {sqrt(MSE * sum(contrast^2 / n))}  # (17.42)
    
    .data <- model.frame(formula, data)
    y  <- .data[, 1]
    x  <- .data[, 2]
    r  <- nrow(cont)

    means <- tapply(y, x, mean)
    L     <- apply(cont, 2, f, mean = means)
    se    <- apply(cont, 2, h, n = tapply(y, x, length), MSE = MSE)
    S     <- sqrt((r-1) * qf(conf.level, r-1, nrow(.data)-r))               # (17.43a)
    L + S*cbind('lwr' = -se, 'upr' = se)                                 # (17.43)
}
cont <- matrix(c(0.5,  0.5, -0.5, -0.5,
                 0.5, -0.5,  0.5, -0.5,
                   1,   -1,    0,    0,
                   0,    0,    1,   -1), nrow = 4)
Scheffe(y ~ x1, .data, cont = cont, conf.level = 0.90, MSE = 10.55)

TukeyHSD(aov(y ~ factor(x1), .data), conf.level = 0.90)                 # For Example 2--Unequal Sample Sizes (p 750)

## ------------------------------------------------------------------------
Bonferroni <- function(formula, data, cont, conf.level = 0.975, MSE) {
    f  <- function(contrast, mean)   {contrast %*% mean}                # (17.41)
    h  <- function(contrast, n, MSE) {sqrt(MSE * sum(contrast^2 / n))}  # (17.42)
    
    .data <- model.frame(formula, data)
    y  <- .data[, 1]
    x  <- .data[, 2]
    r  <- nrow(cont)
    g  <- ncol(cont)

    means <- tapply(y, x, mean)
    L     <- apply(cont, 2, f, mean = means)
    se    <- apply(cont, 2, h, n = tapply(y, x, length), MSE = MSE)
    B     <- qt(1 - (1 - conf.level)/(2*g), nrow(.data)-r)                  # (17.46a)
    L + B*cbind('lwr' = -se, 'upr' = se)                                 # (17.46)
}

# Define contrast matrix
Bonferroni(y ~ x1, .data, cont[, 1:2], conf.level = 0.975, MSE = 10.55)
Scheffe(y ~ x1, .data, cont = cont[, 1:2], conf.level = 0.90, MSE = 10.55)  # Compare with Scheffe
Bonferroni(y ~ x1, .data, cont, conf.level = 0.975, MSE = 10.55)  # Compare to Scheffe example above
Scheffe(y ~ x1, .data, cont = cont, conf.level = 0.90, MSE = 10.55)  # From previous example


## ------------------------------------------------------------------------
r = 4
u = sum(means) / r                          # (17.48a)
n = with(.data, tapply(y, x1, length))         # lengths ni
s = vector("numeric", r)
for(i in seq(r))
  s[i] = (10.55/n[i]) * ((r-1)/r)^2 + (10.55/r^2) * sum(1/n[-i])  # (17.49)

plot(x = seq(r), means, pch = 20, ylim = c(10, 30), xlab = "Levels of Design", ylab = "Mean", xaxt = 'n')
axis(1, seq(4))
segments(seq(4), u, seq(4), means)
lines(seq(1, 4.5, 0.5), rep(u+s, each = 2), type = "S")
lines(seq(1, 4.5, 0.5), rep(u-s, each = 2), type = "S")
abline(h = u)

## ------------------------------------------------------------------------
.data <- CH17TA04
names(.data) <- c("Y", "X1", "X2")
.data <- transform(.data, X1 = factor(X1, labels = c(' 6 hours', ' 8 hours', '10 hours', '12 hours')))

## ------------------------------------------------------------------------
fit <- aov(Y ~ X1, .data)
addmargins(xtabs(Y ~ X1 + X2, .data, sparse = TRUE), 2, list(list('Mean' = mean, 'SD' = sd)))
model.tables(fit, "means")                      # Alternative way to see table of means
summary(fit)                                    # ANOVA summary
summary.lm(fit)                                 # Coefficients are mean difference from "6 hours" mean
TukeyHSD(fit)                                   # 1st 3 'diff' match lm model coefficients
plot(TukeyHSD(fit))
cat("q* =", qtukey(0.95, 4, 24))
diff(with(.data, tapply(Y, X1, mean)))             # Diminishing returns to training (p 764)

## ------------------------------------------------------------------------
x1 <- rep(seq(6, 12, 2), each = 7)                  # numeric encoding of X1
xx <- seq(min(x1), max(x1), len = 200)              # For curve plotting below
.data <- transform(.data, x = c(scale(x1, TRUE, FALSE)))  # Center X1
.data <- transform(.data, X = x1)                         # Store the numeric encoding
(.data <- transform(.data, xsq = x*x))                    # TABLE 17.5

plot(.data$Y ~ x1, pch = 20, ylim = c(30, 65), xlab = "Hours of Training", ylab = "Number of Acceptable Units")
lines(xx, coef(lm(Y ~ X + I(X^2), .data)) %*% rbind(1, xx, xx^2))  # Uses (17.51)

## ------------------------------------------------------------------------
fit <- lm(Y ~ x + xsq, .data)  # (17.52)

anova(fit)                  # (a) Regression Model (17.51)
anova(lm(Y ~ X1, .data))       # (b) ANOVA (17.50)
anova(fit, lm(Y ~ X1, .data))  # (c) ANOVA for Lack of Fit

