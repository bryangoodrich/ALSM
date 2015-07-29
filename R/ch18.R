################################################################################
## Input Rust Inhibitor Data                                                   #
################################################################################
data <- read.table("data/CH17TA02.txt",
  col.names = c("value", "brand", "unit")
);
attach(data)
fit <- lm(value ~ factor(brand)-1)


################################################################################
## TABLE 18.1                                                          (p 777) #
## Residuals--Rust Inhibitor Example                                           #
################################################################################
xtabs(resid(fit) ~ unit + brand)


################################################################################
## FIGURE 18.1                                                         (p 777) #
## Diagnostic Residual Plot--Rust Inhibitor Example                            #
################################################################################
par(ask = TRUE, pch = 16)

plot(resid(fit) ~ fitted(fit), ylab = "Residual", xlab = expression(hat(Y)),
main = expression(paste("(a) Residual against ", hat(Y))))
abline(0,0)

stripchart(split(resid(fit), brand), main = "(b) Aligned Residual Dot Plot",
  method = "stack", pch = 16)

qqnorm(resid(fit), xlab = "Exp Val", ylab = "Residual", 
  main = "(c) Normal Probability Plot")
qqline(resid(fit))


################################################################################
## FIGURE 18.2 - 18.5 Omitted; No data                                         #
################################################################################


################################################################################
## Input ABT Electronics Data                                                  #
################################################################################
detach(data)
data <- read.table("data/CH18TA02.txt", 
  col.names = c("pull", "flux", "joint")
);
attach(data)


################################################################################
## TABLE 18.2                                                          (p 783) #
## Solder Joint Pull Strength--ABT Electronics Example                         #
##                                                                             #
## R has a function for the H distribution used in this example. It comes from #
## the maxFratio (SuppDists) function. The df and k are reversed in its        #
## arguments from those used in the book. Note that (18.10) just is the        #
## variance (var) function applied to each factor level, as used below.        #
################################################################################
library(SuppDists)
round(addmargins(xtabs(pull ~ joint + flux), 1, 
  FUN = list(list(mean, median, s = var))), 3)

(Hstar <- max(tapply(pull, flux, s)) / min(tapply(pull, flux, var)))  ## (18.8)
Hstar <= print(qmaxFratio(0.95, 7, 5))                                ## (18.9)

rm(s, Hstar)


################################################################################
## FIGURE 18.6                                                         (p 784) #
## Dot Plots of Pull Strengths--ABT Electronics Example                        #
################################################################################
stripchart(pull ~ flux, pch = 16, method = "stack", xlim = c(0, 30),
  xlab = "Pull Strength", ylab = "Type")


################################################################################
## TABLE 18.3                                                          (p 785) #
## Absolute Deviations of Responses from                                       #
## Treatment Medians--ABT Electronics Example                                  #
##                                                                             #
## R does contain a function for the Brown-Forsythe test called hov (HH).      #
## Since this test is not nearly as simple as the Hartley test, which is one   #
## of its prime advantages, making use of hov is a time saver. The BF          #
## computations serve no pedagogical value at this time, sot hey will be       #
## omitted.                                                                    #
################################################################################
library(HH)
hov(pull ~ factor(flux))
data.frame(tapply(pull, flux, function(x) list( abs(x - median(x)) )))


################################################################################
## FIGURE 18.7                                                         (p 788) #
## Weighted Regression Output for Full and                                     #
## reduced Models--ABT Electronics Example                                     #
##                                                                             #
## While the model given in (18.19) is a regression on a single predictor with #
## no intercept, this is actually nothing more than a model with just the      #
## intercept term, because when an intercept term is included, it takes a      #
## column vector of '1' equal to Xij as its predictor. This is expressed by    #
## the fact the reduced model fitted below contains nothing but the intercept  #
## term.                                                                       #
################################################################################
w = 1 / tapply(pull, flux, var)  ## Factor level variances (18.10)
w = rep(w, each = 8)             ## Repeat for each treatment in a factor level
fit <- lm(pull ~ 0 + flux, data = transform(data, flux = factor(flux)), 
  weights = w)
summary(fit); anova(fit)
summary(lm(pull ~ 1, weights = w)); anova(lm(pull ~ 1, weights = w))
anova(lm(pull ~ 1, weights = w), fit)


################################################################################
## TABLE 18.4                                                          (p 788) #
## Data for Weighted Least Squares Regression--ABT Electronics Example         #
################################################################################
data.frame(
  i = flux,
  j = joint,
  Y = pull,
  model.matrix(fit),
  Weights = w,
  Reduced.Model = 1
);

rm(w)


################################################################################
## Input the Servo Data                                                        #
################################################################################
detach(data)
data <- read.table("data/CH18TA05.txt", 
  col.names = c("fail", "location", "interval")
);
attach(data)


################################################################################
## TABLE 18.5                                                          (p 791) #
## Time between Computer Failures at Three                                     #
## Locations (in hours)--Servo-Data Example                                    #
################################################################################
addmargins(xtabs(fail       ~ interval + location), 1, list(list(mean, var)))
addmargins(xtabs(rank(fail) ~ interval + location), 1, list(list(mean, var)))

as.table(cbind(
  round(tapply(fail, location, function(x) var(x) / mean(x)), 1),
  round(tapply(fail, location, function(x) sd(x)  / mean(x)), 2),
  round(tapply(fail, location, function(x) sd(x)  / mean(x)^2), 3)
));


################################################################################
## TABLE 18.6 Omitted.                                                 (p 792) #
## Normal Probability Plots for Original and                                   #
## Transformed Data--Servo-data Example                                        #
## TABLE 18.6 Omitted.                                                         #
## See TABLE 3.9 for details. boxcox (MASS) will be used instead.              #
##                                                                             #
## The multiple comparison intervals returned by TukeyHSD are somewhat         #
## different from the book. Almost exacting results can be obtained if one     #
## uses a conf.level = 0.9208, though. Still, the conclusions are the same.    #
################################################################################
library(MASS)
fit <- lm(fail ~ factor(location) - 1)
boxcox(fit)

par(mfrow = c(2, 2), pch = 16)
qqnorm(resid(fit)); qqline(resid(fit))

fit <- update(fit, log(.) ~ .)
qqnorm(resid(fit)); qqline(resid(fit))

  ## Some statistics (p 793)
anova(lm(log(fail) ~ 1), fit)     ## anova(reduced, full)
qf(0.90, 2, 12)                   ## Compare with above F*
summary(fit)                      ## Authors typo mean(location2) = 2.792
TukeyHSD(aov(log(fail) ~ factor(location)), conf.level = 0.90)


################################################################################
## EXAMPLE                                                           (p 795-8) #
## Nonparametric Rank F Test--Servo-Data Example                               #
##                                                                             #
## R contains a nonparametric multiple comparison function npmc (npmc). It is  #
## provided here for a contrast, but whether it is even applicable is unclear  #
## at this time.                                                               #
################################################################################
fit <- lm(rank(fail) ~ factor(location))
addmargins(xtabs(rank(fail) ~ interval + location), 1, list(list(mean, var)))
anova(fit); qf(0.90, 2, 12)
kruskal.test(fail ~ location)     ## See author's comments: Kurskal-Wallis test
qchisq(0.90, 2)                   ## As above, results are borderline, and this
                                  ## test fails to reject the null hypothesis.

  ## Multiple Pairwise Comparisons
B = qnorm(1-0.1/6)*sqrt((15*(15+1))/12 * (1/5 + 1/5))
m = tapply(rank(fail), location, mean)
round(rbind(
  "1-2" = c(diff = m[[1]] - m[[2]], (m[1] - m[2]) + c(lwr = -B,upr =  B)),
  "3-2" = c(diff = m[[3]] - m[[2]], (m[3] - m[2]) + c(lwr = -B,upr =  B)),
  "3-1" = c(diff = m[[3]] - m[[1]], (m[3] - m[1]) + c(lwr = -B,upr =  B))
), 1)
summary(npmc(data.frame(var = fail, class = location), alpha = 0.10, df = 0))

rm(B, m)


################################################################################
## Input Heart Transplant Data                                                 #
################################################################################
detach(data)
data <- read.table("data/CH18TA07.txt", col.names = c("time", "mismatch"))
data <- transform(data,
  mismatch = factor(mismatch, labels = c("Low", "Medium", "High"))
);
attach(data)
fit <- lm(time ~ mismatch - 1)


################################################################################
## TABLE 18.7                                                          (p 798) #
## Survival Times of Patients Following Heart                                  #
## Transplant Surgery--Heart Transplant Example                                #
################################################################################
split(time, mismatch)


################################################################################
## FIGURE 18.9                                                         (p 799) #
## Diagnostic Plots--Heart Transplant Example                                  #
################################################################################
library(HH)
library(MASS)

par(mfrow = c(2, 2))
stripchart(time ~ mismatch, method = "stack", pch = 16, xlab = "Survival Time",
  main = "(a) Dot Plots of Survival Times")
stripchart(rstudent(fit) ~ mismatch, method = "stack", pch = 16,
  main = "(b) Dot Plots of Studentized Residuals",
  xlab = "Studentized Residual")
qqnorm(rstudent(fit), pch = 16,
  xlab = "Expected Value", ylab = "Studentized Residual",
  main = "(c) Normal Probability Plot\n of Studentized Residuals")

hov(time ~ mismatch)     ## (p 799) test for constancy of variance
boxcox(fit)


################################################################################
## FIGURE 18.10                                                        (p 800) #
## Diagnostic Plots and ANOVA Table for                                        #
## Transformed Data--Heart Transplant Example                                  #
################################################################################
fit <- update(fit, log(.) ~ .)
summary(aov(log(time) ~ mismatch))

par(mfrow = c(2, 2))
stripchart(rstudent(fit) ~ mismatch, method = "stack", pch = 16,
  main = "(b) Dot Plots of Studentized Residuals",
  xlab = "Studentized Residual")
qqnorm(rstudent(fit), pch = 16,
  xlab = "Expected Value", ylab = "Studentized Residual",
  main = "(c) Normal Probability Plot\n of Studentized Residuals")


################################################################################
## Clean up the R environment from this session                                #
################################################################################
detach(data)
rm(data, fit)



