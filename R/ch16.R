################################################################################
## FIGURE 16.2                                                         (p 680) #
## Analysis of Variance Model Representation--Incentive Pay Example            #
##                                                                             #
## This example is not computationally interesting, but I present a brute      #
## force way of depicting the representation in R given the relevant           #
## information provided on page 683.                                           #
################################################################################
curve(dnorm(x, mean = 58, sd = 4), 40, 110, ylim = c(0, .2),
  xlab = "", ylab = "", col = "green")
text(58, .11, "Type 2", cex = .75) 
curve(dnorm(x, mean = 70, sd = 4), add = TRUE, col = "blue")
text(70, .11, "Type 1", cex = .75)
curve(dnorm(x, mean = 84, sd = 4), add = TRUE, col = "red")
text(84, .11, "Type 4", cex = .75)
curve(dnorm(x, mean = 90, sd = 4), add = TRUE)
text(90, .11, "Type 3", cex = .75)
points(c(51, 78), c(0, 0), pch = 19)


################################################################################
## Input Kenton Food Company Data                                              #
################################################################################
data <- read.table("data/CH16TA01.txt", 
  colClasses = c("numeric", "factor", "factor"),
  col.names = c("sold", "package", "store"))
attach(data)
fit  <- lm(sold ~ package - 1)  ## This is the cell means model discussed below


################################################################################
## TABLE 16.1                                                          (p 686) #
## Number of Cases Sold by Stores for Each of Four                             #
## Package Designs--Kenton Food Company Example                                #
##                                                                             #
## Since xtabs automatically fills in the empty cells to 0, using length to    #
## obtain counts will not work in the addmargins function list. Instead, I do  #
## a separate call leaving out the response term of the formula. Using the     #
## single-sided formula produces a table of counts. In this case, we can then  #
## sum the rows to get 'n'.                                                    #
################################################################################
addmargins(xtabs(sold ~ package + store), 2, list(list(Sum = sum, Mean = mean)))
addmargins(xtabs(     ~ package + store), 2, list(n = sum))
coef(fit)     ## Coefficients in cell means model are the factor means


################################################################################
## FIGURE 16.3                                                         (p 686) #
## Plot of Number of Cases Sold by Package Design--Kenton Food Company Example #
################################################################################
library(lattice)
xyplot(sold ~ package, data = data, groups = store,
  pch = 19, auto.key = list(columns = 5))


################################################################################
## TABLE 16.2                                                          (p 689) #
## Residuals--Kenton Food Company Example                                      #
################################################################################
round(addmargins(xtabs(resid(fit) ~ package + store), 2), 2)


################################################################################
## FIGURE 16.5                                                         (p 695) #
## Output for Single-Factor Analysis of Variance--Kenton Food Company Example  #
##                                                                             #
## Note that "Root Mean Square Error" is just the "Residual Standard Error" in #
## the aov output. Also, "C. Total" is just the aggregate.                     #
################################################################################
anova(fit); aov(fit); summary(fit); confint(fit)
summary(fit)$f[1] <= qf(1-.05, 4-1, 19-4)  ## F-test; Conclude H0?


################################################################################
## TABLE 16.4                                                       (p 707-12) #
## Regression Approach to the Analysis of                                      #
## Variance--Kenton Food Company Example                                       #
## EXAMPLE                                                                     #
## Factor Effects Model with Weighted Means--Kenton Food Company Example       #
## EXAMPLE                                                                     #
## Cell Means Model--Kenton Food Company Example                               #
################################################################################
X1 = X2 = X3 = ifelse(package == 4, -1, 0)
X1[package == 1] <- 1
X2[package == 2] <- 1
X3[package == 3] <- 1

cbind(i = package, j = store, Y = sold, X1, X2, X3)[order(package), ]
(fit <- lm(sold ~ X1 + X2 + X3)); summary(fit); anova(fit)


  ## Factor effects model with weighted means
X1 = ifelse(package == 1, 1, 0); X1[package == 4] <- -5/5
X2 = ifelse(package == 2, 1, 0); X2[package == 4] <- -5/5
X3 = ifelse(package == 3, 1, 0); X3[package == 4] <- -4/5

cbind(i = package, j = store, Y = sold, X1, X2, X3)[order(package), ]
(fit <- lm(sold ~ X1 + X2 + X3)); summary(fit); anova(fit)


  ## Cell means model
fit <- lm(sold ~ package - 1)  ## This is the original model fitted
summary(fit); anova(fit)

rm(X1, X2, X3)


################################################################################
## TABLE 16.5                                                          (p 715) #
## Randomization Samples and Test Statistics--Quality Control Example          #
## FIGURE 16.8                                                                 #
## Randomization Distribution of F* and Corresponding                          #
## F Distribution--Quality Control Example                                     #
##                                                                             #
## Since there is no algorithm to compute this example we had to devise one.   #
## It should come as rather straight-forward. The Xi's are as in the above     #
## examples. The 'y' will hold the 1,680 cases of 9-sequences consisting of    #
## the response variables. The 'ti' implies the treatment group. In this case  #
## t1 is the first group (3-sequence) and t12 is the composite of t1 and t2.   #
## The 'remainder' function is a wrapper for grabing a subset of 'set' based   #
## on those values not in 'x'. The 'seq6' is the 6-sequence remainder after t1 #
## is defined. The whole process took less than 10 seconds on a 2.4 GHz        #
## processor. As for the output, the columns are arbitrarily labeled 1-9.      #
## Clearly they represent the three treatment groups based on groups of three. #
## The function 'f' uses the matrix algebra discussed in Ch. 5. It is possible #
## to get away with merely fitting an 'lm' object, and then extract the        #
## f-statistic in a single call. However, this requires a lot of additional    #
## work for each of the 1680 rows. It took somewhere between 30-60 seconds to  #
## produce the same result.                                                    #
################################################################################
remainder <- function(x, set) set[!set %in% x]
f <- function(Y, X) {
  Y <- matrix(Y)                                ## Turn row-vector into column
  p <- ncol(X);     n <- nrow(X)
  J <- matrix(1, n, n)                          ## (5.18)
  H <- X %*% solve(t(X) %*% X) %*% t(X)         ## (5.73a)
  SSE <- t(Y) %*% (diag(n) - H) %*% Y           ## (5.89b)
  SSR <-  t(Y) %*% (H - (1/n)*J) %*% Y          ## (5.89c)
  fstar <- (SSR / (p - 1)) / (SSE / (n - p))    ## (6.39b)
}

base <- c(1.1, 0.5, -2.1, 4.2, 3.7, 0.8, 3.2, 2.8, 6.3)
t2   <- t12 <- t123 <- list()
y    <- NULL
X    <- cbind(
  X1 = c(1, 1, 1, 0, 0, 0, 0, 0, 0), 
  X2 = c(0, 0, 0, 1, 1, 1, 0, 0, 0),
  X3 = c(0, 0, 0, 0, 0, 0, 1, 1, 1)
);

t1   <- t(combn(base, 3)) 
seq6 <- t(combn(base, 3, remainder, set = base))

for (i in 1:84)  t2[[i]] <- t(combn(seq6[i, ], 3))
for (i in 1:84) t12[[i]] <- cbind(t1[i, 1], t1[i, 2], t1[i, 3], t2[[i]])
for (i in 1:84) 
  t123[[i]] <- cbind(t12[[i]], t(apply(t12[[i]], 1, remainder, set = base)))
for (i in 1:84) y <- rbind(y, t123[[i]])

fstar <- apply(y, 1, function(Y) f(Y, X))

cbind(y, data.frame(f = fstar)) 
hist(fstar, freq = FALSE, ylim = c(0, 1), col = "gray90", main = "")
curve(df(x, 2, 6), add = TRUE, lwd = 2)

rm(base, fstar, i, remainder, seq6, t1, t2, t12, t123, f, X, y)


################################################################################
## Clean up the R environment from this session                                #
################################################################################
detach(data)
rm(data, fit)


