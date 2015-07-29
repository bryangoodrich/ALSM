################################################################################
## FIGURE 13.1                                                         (p 512) #
## Plots of Exponential and Logistic Response Functions                        #
################################################################################
par(mfrow = c(1, 2))
curve(100 - 50*exp(-2*x),
  xlim = c(0, 4), ylim = c(50, 100), xlab = "", ylab = "")
curve(10/(1 + 20*exp(-2*x)),
  xlim = c(0, 4), ylim = c(0, 10), xlab = "", ylab = "")


################################################################################
## Input Severely Injured Patients Data                                        #
##                                                                             
## Choice of start values discussed on (p 521).                                #
################################################################################
data <- read.table("data/CH13TA01.txt", col.names = c("index", "days"))
attach(data)
fit <- nls(index ~ a * exp(b * days), start = c(a = 56.6646, b = -0.03797))


################################################################################
## TABLE 13.1                                                          (p 515) #
## Data--Severely Injured Patients Example                                     #
################################################################################
cbind(Patient = 1:15, "Days Hospitalized" = days, "Prognostic Index" = index)


################################################################################
## FIGURE 13.2                                                         (p 515) #
## Scatter Plot and Fitted Nonlinear Regression                                #
## Function--Severely Injured Patients Example                                 #
################################################################################
plot(index ~ days, pch = 16)
curve(coef(fit)[1] * exp(coef(fit)[2] * x), add = TRUE)


################################################################################
## TABLE 13.2                                                        (p 522-4) #
## Y(0) and D(0) Matrices--Severely Injured Patients Example                   #
## TABLE 13.3                                                                  #
## Gauss-Newton Method Iterations and Final Nonlinear                          #
## Least Squares Estimates--Severely Injured Patients Example                  #
##                                                                             #
## Below is a function suited for this example. It outputs part (a) under the  #
## list items "Estimate". Under "LS" we have a column of the final chosen 'g'  #
## and a column for the MSE and final SSE. In "se" is a reproduction of the    #
## matrix for part (c). Included are also the fitted values, which can be      #
## matched against fitted(fit), and the Y and D matrices for TABLE 13.2.       #
## However, these are the final matrices, not the first. The approach here can #
## manually be applied to reproduce TABLE 13.2 exactly. Note, the standard     #
## errors in (b) can be obtained from (c)                                      #
################################################################################
nlm <- function(x, y, g, FUN = function(x, a, b) a * exp(b * x), TOL = 1e-5) {
  n <- length(y)
  SSE <- NULL
  i <- 1
  G <- matrix(nrow = 50, ncol = 2)

  repeat {
    G[i, ] = g
    Y   = as.matrix(y - FUN(x, g[1], g[2]))               ## (13.25a)
    D   = cbind(exp(g[2] * x), g[1] * x * exp(g[2] * x))  ## (13.25b)
    b   = coef(lm(Y ~ D - 1))                             ## (13.25c)
    SSE = c(SSE, sum(Y^2))
    g <- g + b
    is.done = ((length(SSE) > 1) && (SSE[i-1] - SSE[i] < TOL))
    i = i + 1
    if (is.done) break
  } ## End repeat

  G <- G[!is.na(G[,1]), ]
  i = i-1
  tab <- list()
  tab$Estimates = cbind(G, SSE)
  tab$LS = cbind(G[i, ], rbind(SSE[i-1] / (n-2), SSE[i-1]))
  tab$se = (SSE[i] / (n-2)) * solve(t(D) %*% D)
  tab$fitted = G[i, 1] * exp(G[i, 2] * x)
  tab$Y = Y
  tab$D = D
  return(tab)
} ## end nonlinear model function

(coefs <- coef(lm(log(index) ~ days)))  ## Linear Transform Fit
(coefs <- c(exp(coefs[1]), coefs[2]))   ## Exp. Transformation
nlm(days, index, g = coefs)

rm(coefs, nlm)


################################################################################
## FIGURE 13.3                                                         (p 527) #
## Diagnostic Residual Plots--Severely Injured Patients Example                #
################################################################################
library(HH)     ## To test nonconstant error variance with Brown-Forsythe
par(ask = TRUE, pch = 16)
plot(resid(fit) ~ fitted(fit), xlab = "Fitted Values", ylab = "Residual",
  main = "(a) Residual Plot against Fitted"); abline(0, 0)
qqnorm(resid(fit), xlab = "Expected Value", ylab = "Residual", 
  main = "(b) Normal Probability Plot"); abline(0, 0)

hov(resid(fit) ~ cut(days, 2))       ## BF test. P-value = 0.64
deviance(fit) / df.residual(fit)     ## (13.31)
vcov(fit)                            ## (13.32b)


################################################################################
## FIGURE 13.4                                                         (p 531) #
## Bootstrap Sampling Distribution--Severely Injured Patients Example          #
##                                                                             #
## The plotting limits could be reduced so the output looks more like those in #
## the text.                                                                   #
################################################################################
library(boot)
boot.coef <- function(data, indices, g) {
  x   <- data[indices, 2]
  y   <- data[indices, 1]
  mod <- nls(y ~ a * exp(b*x), start = c(a = g[[1]], b = g[[2]]))
  return(coef(mod))
}  ## end bootstrap statistic function

(z <- boot(data, boot.coef, R = 1000, g = coef(lm(log(index) ~ days))))

  ## NLS coefficients and Bootstrap CIs
coef(fit); boot.ci(z, index = 1); boot.ci(z, index = 2)

par(ask = TRUE)
g0 = expression(g[0]^symbol("\052")); g1 = expression(g[1]^symbol("\052"))
hist(z$t[, 1], 40, freq = FALSE, xlab = g0, main = expression(
  "(a) Histogram of Bootstrap Estimates " * g[0]^symbol("\052")))
hist(z$t[, 2], 40, freq = FALSE, xlab = g1, main = expression(
  "(b) Histogram of Bootstrap Estimates " * g[1]^symbol("\052")))

s = sqrt(diag(vcov(fit)))         ## Simultaneous Interval Estimation
(coef(fit) + 2.16 * c(lwr = -s, upr = s))[c(1,3, 2,4)]

rm(s, z, boot.coef, g0, g1)


################################################################################
## Input Learning Curve Data                                                   #
################################################################################
detach(data)
data <- read.table("data/CH13TA04.txt", 
  col.names = c("location", "week", "efficiency"))
attach(data)


################################################################################
## TABLE 13.4                                                          (p 534) #
## Data--Learning Curve Example                                                #
################################################################################
xtabs(efficiency ~ week + location)  ## Easier to view than "print(data)"


################################################################################
## TABLE 13.5                                                          (p 535) #
## Nonlinear Least Squares Estimates and Standard Deviations                   #
## and Bootstrap Results--Learning Curve Example                               #
################################################################################
library(boot)
boot.coef <- function(data, indices, g) {
  y   <- data[indices, 1]
  x1  <- data[indices, 2]
  x2  <- data[indices, 3]
  mod <- nls(y ~ a  + b*x1 + c*exp(d*x2), start = g)
  return(coef(mod))
}  ## end bootstrap statistic function

g <- c(a = 1.025, b = -0.0459, c = -0.5, d = -0.122)
fit <- nls(efficiency ~ a + b*location + c*exp(d*week), start = g)
z <- boot(data[c(3, 1:2)], boot.coef, R = 1000, g = g)
cbind(
    g   = g,
  "(1)" = summary(fit)$coef[, 1],
  "(2)" = summary(fit)$coef[, 2],
  "(3)" = apply(z$t, 2, mean),
  "(4)" = apply(z$t, 2, sd)
); 

  ## Residual Plots Not Shown
par(mfrow = c(2, 2))
plot(resid(fit) ~ fitted(fit), main = "Residuals aganist Fitted")
plot(resid(fit) ~ location, main = "Residuals against Location")
plot(resid(fit) ~ week, main = "Residuals against Week")
qqnorm(resid(fit)); qqline(resid(fit))

rm(boot.coef, g)


################################################################################
## FIGURE 13.5                                                         (p 534) #
## Scatter Plot and Fitted Nonlinear Regression                                #
## Functions--Learning Curve Example                                           #
################################################################################
p <- coef(fit)
plot(efficiency ~ week, data[location == 0, ])
points(efficiency ~ week, data[location == 1, ], pch = 16)
curve(p[1] + p[2]*(0) + p[3] * exp(p[4] * x), to = 90, add = T) ## Location = 0
curve(p[1] + p[2]*(1) + p[3] * exp(p[4] * x), to = 90, add = T) ## Location = 1

rm(p)


################################################################################
## FIGURE 13.6                                                         (p 536) #
## Histograms of Bootstrap Sampling Distributions--Learning Curve Example      #
##                                                                             #
## Fullscreen if the plot is too narrow. We included a bottom row of density   #
## plots that give the same information, but in someways, an easier format to  #
## view. 
################################################################################
par(mfrow = c(2,4))
for(i in seq(4)) {
  hist(z$t[, i], 100, freq = FALSE)
  abline(v = mean(z$t[, i]), lwd = 2)}
for(i in seq(4)) plot(density(z$t[, i]))

rm(z, i)


################################################################################
## FIGURE 13.7                                                         (p 539) #
## Various Logistic Activation Functions for Single Predictor                  #
##                                                                             #
## The scale cannot possibly be the same for all lines plotted. When, on this  #
## notation, b = 0.1 the range of x values needs to be around -50:50 to see    #
## f(x) values near 0 and 1 as shown in (a).                                   #
################################################################################
f <- function(x, a, b) (1 + exp(-a - b*x))^(-1)

par(ask = TRUE, lwd = 2)
curve(f(x, a =  0, b = 0.1), -10, 10, ylim = c(0, 1))
curve(f(x, a =  0, b =   1), -10, 10, add = TRUE, lty = 5)
curve(f(x, a =  0, b =  10), -10, 10, add = TRUE, lty = 3)

curve(f(x, a =  0, b = -0.1), -10, 10, ylim = c(0, 1))
curve(f(x, a =  0, b =   -1), -10, 10, add = TRUE, lty = 5)
curve(f(x, a =  0, b =  -10), -10, 10, add = TRUE, lty = 3)

curve(f(x, a =  5, b = 1), -10, 10, ylim = c(0, 1))
curve(f(x, a =  0, b = 1), -10, 10, add = TRUE, lty = 5)
curve(f(x, a = -5, b = 1), -10, 10, add = TRUE, lty = 3)

rm(f)


################################################################################
## Input Ischemic Heart Disease (IHD) Data                                     #
################################################################################
detach(data)
data <- read.table("data/APPENC09.txt", 
  col.names = c("id", "cost", "age", "gender", "intervention", "drugs", 
                "visits", "complications", "comobidities", "duration")
);
attach(data)


################################################################################
## The rest of this tutorial will remain absent for now. The 'nnet' library    #
## listed below includes methods for developing neural network models. We      #
## leave it to the interested reader to learn more about this topic on their   #
## own, at least for now.                                                      #
################################################################################
library(nnet)





################################################################################
## Clean up the R environment from this session                                #
################################################################################
detach(data)
rm(data, fit)

