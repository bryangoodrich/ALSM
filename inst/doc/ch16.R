## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
library(lattice)       # xyplot
library(RColorBrewer)  # brewer.pal

data("CH16TA01", package = "ALSM")

## ------------------------------------------------------------------------
curve(dnorm(x, mean = 58, sd = 4), 40, 110, ylim = c(0, .2), xlab = "", ylab = "", col = "green")
text(58, .11, "Type 2", cex = .75)

curve(dnorm(x, mean = 70, sd = 4), add = TRUE, col = "blue")
text(70, .11, "Type 1", cex = .75)

curve(dnorm(x, mean = 84, sd = 4), add = TRUE, col = "red")
text(84, .11, "Type 4", cex = .75)

curve(dnorm(x, mean = 90, sd = 4), add = TRUE)
text(90, .11, "Type 3", cex = .75)
points(c(51, 78), c(0, 0), pch = 19)

## ------------------------------------------------------------------------
.data  <- CH16TA01
names(.data) <- c("y", "x1", "x2")

## ------------------------------------------------------------------------
cbind('Table' = addmargins(xtabs(y ~ x1 + x2, .data), 2),
      'Mean'  = tapply(.data$y, .data$x1, mean),
      'n'     = tapply(.data$y, .data$x1, function(r) sum(r > 0)))
                                    
with(.data, c("Y.."   = sum(tapply(y, x1, sum)), 
           "Ybar." = mean(tapply(y, x1, mean)),
           "n.."   = sum(tapply(y, x1, function(r) sum(r > 0)))))

## ------------------------------------------------------------------------
pal <- brewer.pal(5, "Set1")
xyplot(y ~ factor(x1), .data, groups = x2, auto.key = list(columns = 5), 
       par.settings = simpleTheme(col = pal, pch = 19), 
       xlab = "Package Design", ylab = "Cases Sold", main = "Summary Plot")

## ------------------------------------------------------------------------
.data <- transform(.data, u = unlist(tapply(y, x1, scale, scale = FALSE)))
addmargins(xtabs(u ~ x1 + x2, .data, sparse = TRUE), 2)  # Their sums are as 0 as it gets in R.

## ------------------------------------------------------------------------
.data <- transform(.data, x1 = factor(x1))
fit <- lm(y ~ x1 - 1, .data)  # This is the cell means model

anova(fit)
summary(aov(y ~ factor(x1) - 1, .data))  # Same as anova(fit)

summary(fit)  # Notice that the coefficients are just the group means
summary.lm(aov(fit))  # Same as summary(fit)

confint(fit)
summary(fit)$f[1] <= qf(1-.05, 4-1, 19-4)  # F-test; Conclude H0?

## ------------------------------------------------------------------------
# ANOVA as Regression Model (16.79)
contrasts(.data$x1) <- matrix(c(1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 1, -1), 4, 3)
fit <- lm(y ~ x1, .data)
model.matrix(fit)
summary(fit)
anova(fit)


# ANOVA as Factor Effects Model with Weighted Means (16.82)
contrasts(.data$x1) <- matrix(c(1, 0, 0, -1, 0, 1, 0, -1, 0, 0, 1, -0.8), 4, 3)
fit <- lm(y ~ x1, .data)
model.matrix(fit)
summary(fit)
anova(fit)


# ANOVA as Cell Means Model (16.85)
contrasts(.data$x1) <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 4, 4)
fit <- lm(y ~ x1 - 1, .data)  # This is the original model fitted
model.matrix(fit)
summary(fit)
anova(fit)

## ------------------------------------------------------------------------
remainder <- function(x, set) set[!set %in% x]
f <- function(Y, X) {
  Y <- matrix(Y)                                # Turn row-vector into column
  p <- ncol(X)
  n <- nrow(X)
  J <- matrix(1, n, n)                          # (5.18)
  H <- X %*% solve(t(X) %*% X) %*% t(X)         # (5.73a)
  SSE <- t(Y) %*% (diag(n) - H) %*% Y           # (5.89b)
  SSR <-  t(Y) %*% (H - (1/n)*J) %*% Y          # (5.89c)
  fstar <- (SSR / (p - 1)) / (SSE / (n - p))    # (6.39b)
}

base <- c(1.1, 0.5, -2.1, 4.2, 3.7, 0.8, 3.2, 2.8, 6.3)
t2   <- t12 <- t123 <- list()
y    <- NULL
X    <- cbind(
  X1 = c(1, 1, 1, 0, 0, 0, 0, 0, 0),
  X2 = c(0, 0, 0, 1, 1, 1, 0, 0, 0),
  X3 = c(0, 0, 0, 0, 0, 0, 1, 1, 1))

t1   <- t(combn(base, 3))
seq6 <- t(combn(base, 3, remainder, set = base))

for (i in 1:84)  t2[[i]] <- t(combn(seq6[i, ], 3))
for (i in 1:84) t12[[i]] <- cbind(t1[i, 1], t1[i, 2], t1[i, 3], t2[[i]])
for (i in 1:84)
  t123[[i]] <- cbind(t12[[i]], t(apply(t12[[i]], 1, remainder, set = base)))
for (i in 1:84) y <- rbind(y, t123[[i]])

fstar <- apply(y, 1, function(Y) f(Y, X))

hist(fstar, freq = FALSE, ylim = c(0, 1), col = "gray90", main = "")
curve(df(x, 2, 6), add = TRUE, lwd = 2)

# LAST EXAMPLE FOR CHAPTER, fyi
# BIG OUTPUT TO FOLLOW! 
cbind(y, data.frame(f = fstar))

