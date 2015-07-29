################################################################################
## Input Positively Autocorrelated Error Terms                                 #
## TABLE  12.1                                                         (p 482) #
## FIGURE 12.1                                                         (p 483) #
## Example of Positively Autocorrelated Error Terms                            #
##                                                                             #
## Since this data is entirely constructed I will define a function that will  #
## build the true regression output and the error terms given a set of         #
## distrubances (u). It will print out the error terms against time (Xt) while #
## accepting additional plotting parameters.                                   #
##                                                                             #
## Note that the interested reader could try other disturbances from the       #
## standard normal distribution by using rnorm, which defaults to mean = 0 and #
## sd = 1. Depending on the size required, the function will build a dataframe #
## of that size.                                                               #
################################################################################
error <- function(u, estart = 3.0, plot = FALSE, ...) {
  e  <- vector(mode = "numeric", length = length(u))  
  y  <- vector(mode = "numeric", length = length(u))
  Xt <- 0:10

  for(i in seq(u)) ifelse(i != 1, e[i] <- e[i-1] + u[i], e[i] <- estart)
  y <- 2 + 0.5 * Xt + e

  if (isTRUE(plot)) plot(e, pch = 19, ...)

  return( data.frame(u = u, e = e, y = y) )
} ## end function
 
u <- c(0, 0.5, -0.7, 0.3, 0, -2.3, -1.9, 0.2, -0.3, 0.2, -0.1)
(data <- error(u, xlab = "Time", plot = TRUE, ylim = c(-3, 4))); abline(0, 0)

rm(u)


################################################################################
## FIGURE 12.2                                                         (p 483) #
## Regression with Positively Autocorrelated Error Terms                       #
##                                                                             #
## Since (c) uses "different disturbances" the results will differ each time.  #
################################################################################
par(ask = TRUE, pch = 19)

with(data,
  plot(y, ylim = c(0, 10), main = "(a) True Regression Line e0 = 3")
); abline(2, 0.5)

with(data,
  plot(y, ylim = c(0, 10), main = "(b) Fitted Regression Line e0 = 3")
); abline( lm(y ~ x, data = transform(data, x = 0:10)) )

data <- error(rnorm(11), estart = -0.2)
with(data, 
  plot(y, ylim = c(0, 10), main = "(c) Fitted Regression Line e0 = -0.2")
); abline( lm(y ~ x, data = transform(data, x = 0:10)) )

rm(error)


################################################################################
## Input the Blaisdell Company Data                                            #
##                                                                             #
## Since R has extensive facilities for dealing with time series data, I will  #
## convert the data to a time series (ts) object with the appropriate time     #
## interval. For more see the Vito Ricci's reference card at:                  #
##   http://cran.r-project.org/doc/contrib/Ricci-refcard-ts.pdf                #
##                                                                             #
## It will be noted now that throughout the rest of this chapter there are two #
## references to many of the terms, either by t or t-1. Since R is inherently  #
## vectorized and each of these vectors include n-1 terms from the n-sized     #
## vector, we can unambiguously define t = 2:20 to which t-1 = 1:19 by the     #
## vectorized arithmetic. This will make a lot of the notation clearer.        #
################################################################################
data <- read.table("data/Ch12TA02.txt", col.names = c("Y", "X"))
attach(data)
fit  <- lm(Y ~ X)
t    <- 2:20


################################################################################
## FIGURE 12.3                                                         (p 489) #
## Residuals Plotted against Time--Blaisdell Company Example                   #
################################################################################
plot(resid(fit), xlab = "Time", ylim = c(-0.5, .5), pch = 19); abline(0,0)


################################################################################
## TABLE 12.2                                                          (p 489) #
## Data, Regression Results, and Durbin-Watson Test Calculations--Blaisdell    #
## Company Example (Company and Industry Sales Data Are Seasonally Adjusted)   #
##                                                                             #
## R contains a Durban-Watson (dwtest) test function (lmtest). See below. It   #
## is possible to cbind the table and the data object, but it would be         #
## superfluous. I then calculate D which can be checked against Table B.7 for  #
## n = 20, alpha = 0.1, and p = 2. The test conclusion is the same as dwtest.  #
################################################################################
library(lmtest)

round(as.table(cbind(
  '(1)' = Y,
  '(2)' = X,
  '(3)' = resid(fit),
  '(4)' = c(NA,  resid(fit)[t] - resid(fit)[t-1]),
  '(5)' = c(NA, (resid(fit)[t] - resid(fit)[t-1])^2),
  '(6)' = resid(fit)^2
)), 4); 

sum((resid(fit)[2:20] - resid(fit)[1:19])^2) / sum(resid(fit)^2)  ## D = 0.735

dwtest(fit)        ## alternative hypotheses defaults to p > 0


################################################################################
## TABLE 12.3                                                          (p 493) #
## Calculations for Estimating p with the Cochrane-Orcutt                      #
## Procedure--Blaisdell Company Example                                        #
################################################################################
as.table(cbind(
  '(1)' = resid(fit),
  '(2)' = c(NA, resid(fit)[t-1]),
  '(3)' = c(NA, resid(fit)[t-1] * resid(fit)[t]),
  '(4)' = c(NA, resid(fit)[t-1]^2)
)); 

cat("rho =", 
  sum(resid(fit)[t-1] * resid(fit)[t]) / sum(resid(fit)[t-1]^2), "\n"
); 


################################################################################
## TABLE 12.4                                                          (p 493) #
## Transformed Variables and Regression Results for First Iteration with       #
## Cochrane-Orcutt Procedure--Blaisdell Company Example                        #
################################################################################
library(lmtest)
cochrane.orcutt <- function(model) {
  x      <- model.matrix(model)[, -1]
  y      <- model.response(model.frame(model))
  e      <- resid(model)
  n      <- length(e)
  t      <- 2:n
  r      <- sum(e[t-1] * e[t]) / sum(e[t-1]^2)     ## (12.22)    Step 1
  y      <- y[t] - r * y[t-1]                      ## (12.18a)
  x      <- x[t] - r * x[t-1]                      ## (12.18b)
  model  <- lm(y ~ x)                              ## (12.19)    Step 2
  return(model)
}

cbind(
  model.frame(fit),
  rbind(c(NA, NA), model.frame(cochrane.orcutt(fit)))
); summary(cochrane.orcutt(fit)); anova(cochrane.orcutt(fit))

dwtest(cochrane.orcutt(fit))$p.value   ## DW = 1.6502, fail to reject H0

  ## Convert transformed coefficients back to original
  ## Standard error will be ignored.
r <- coef(lm(resid(fit)[2:20] ~ resid(fit)[1:19] - 1))    ## Alternate r calc.
coefs <- coef(cochrane.orcutt(fit))
c(b0 = coefs[1] / (1 - r), b1 = coefs[2])

rm(cochrane.orcutt, r, coefs)


################################################################################
## TABLE 12.5                                                          (p 495) #
## Hildreth-Lu Results--Blaisdell Company Example                              #
################################################################################
rho <- c(0.1, 0.3, 0.5, 0.7, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98)

hildreth.lu <- function(rho, model) {
  x <- model.matrix(model)[, -1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[t] - rho * y[t-1]
  x <- x[t] - rho * x[t-1]
  return(lm(y ~ x))
} ## end function

cbind(
  rho,
  SSE = apply(data.frame(rho), 1,
    function(r) anova(hildreth.lu(r, fit))[2, 2]
  ) ## end apply function
);     ## minimum SSE = 0.07167118 for rho = 0.96

summary(hildreth.lu(0.96, fit)); dwtest(hildreth.lu(0.96, fit))

cat("Y = ", coef(hildreth.lu(0.96, fit))[1] / (1 - 0.96), "+",  
  coef(hildreth.lu(0.96, fit))[2], "\n"                         ## (12.28)
); 

rm(rho, hildreth.lu)


################################################################################
## TABLE 12.6                                                          (p 497) #
## First Differences and Regression Results with First                         #
## Differences Procedure--Blaisdell Company Example                            #
################################################################################
fit <- lm(y ~ x - 1, data.frame(y = Y[t] - Y[t-1], x = X[t] - X[t-1]))
cbind(
  Y = Y, 
  X = X, 
  y = c(NA, Y[t] - Y[t-1]),
  x = c(NA, X[t] - X[t-1])
); 

cat("Y = ", mean(Y) - coef(fit)[1]*mean(X), "+", coef(fit)[1], "\n") ## (12.33)

rm(t)


################################################################################
## TABLE 12.7                                                          (p 498) #
## Major Regression Results for Three Transformation                           #
## Procedures--Blaisdell Company Example                                       #
##                                                                             #
## Ignored for lack of pedagogical value. Just read it.                        #
################################################################################

################################################################################
## FORECAST EXAMPLE                                                  (p 500-1) #
##                                                                             #
## This forecast was manually calculated from the Cochrane-Orcutt modified     #
## (but not transformed) regression model. The steps lack pedagogical value    #
## for using R since we're just using it as a calculator in that case.         #
##                                                                             #
## This chapter was a very brief introduction to time series analysis with R.  #
## I suggest looking at "Time Series Analysis and Its Applicatioins" by Robert #
## H. Shumway and David S. Stoffer. The text website for it can be found at:   #
##   http://www.stat.pitt.edu/stoffer/tsa3/                                    #
################################################################################


################################################################################
## Clean up the R environment from this session                                #
################################################################################
detach(data)
rm(data, fit)



