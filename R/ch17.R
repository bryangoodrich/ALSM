################################################################################
## TABLE 17.1                                                                  #
## Summary of Results--Kenton Food Company Example                             #
## Ignored due to lack of pedagogical benefit. See Chapter 16 for details.     #
################################################################################


################################################################################
## Input Rust Inhibitor Data                                                   #
################################################################################
data <- read.table("data/CH17TA02.txt",
  col.names = c("value", "brand", "unit")
);
attach(data)
fit <- lm(value ~ factor(brand))


################################################################################
## TABLE 17.2                                                          (p 735) #
## Data and Analysis of Variance Results--Rust Inhibitor Example               #
################################################################################
addmargins(xtabs(value ~ unit + brand), 1, mean); mean(value)
anova(fit)
anova(fit)[1,3] / anova(fit)[2,3]     ## F* = MSTR / MSE = 866.12
qf(0.95, 3, 36)                       ## F(0.95, 3, 36) = 2.87 < F*


################################################################################
## Input Kenton Food Company Data                                              #
################################################################################
detach(data)
data <- read.table("data/CH16TA01.txt", 
  colClasses = c("numeric", "factor", "factor"),
  col.names = c("sold", "package", "store"))
attach(data)


################################################################################
## FIGURE 17.1                                                         (p 736) #
## Line Plot of Estimated Factor Level Means--Kenton Food Company Example      #
##                                                                             #
## I don't know of a line plot function, but a makeshift way to generate an    #
## identical graph is to use the stripchart function (graphics).               #
################################################################################
stripchart(by(sold, package, mean), pch = 19,
  ylim = c(1,2), xlim = c(12, 30), xlab = "Cases Sold")
text(by(sold, package, mean), rep(1,4), labels = 1:4, pos = 3)


################################################################################
## FIGURE 17.2                                                         (p 736) #
## Bar Graph and Main Effects Plot of Estimated                                #
## Factor Level Means--Kenton Food Company Example                             #
## FIGURE 17.3                                                         (p 739) #
## Bar-Interval Graph and Interval Plot--Kenton Food Company Example           #
################################################################################
library(lattice)
means <- by(sold, package, mean)
diff = qt(0.975, 15) * sqrt(10.55 / by(sold, package, length))

par(ask = TRUE)
barplot(means, xlab = "Design", ylab = "Cases Sold",
  main = "(a) Bar Graph")

plot(means, type = "o", pch = 19, ylim = c(0, 30),
  xlab = "Design", ylab = "Cases Sold", main = "(b) Main Effects Plot")

bar <- barplot(means, xlab = "Design", ylab = "Cases Sold", ylim = c(0, 32),
  main = "(a) Bar-Interval Graph")
arrows(bar, means+diff, bar, means-diff, angle = 90, code = 3)

plot(means, pch = 19, ylim = c(0, 30), xlim = c(.5, 4.5),
  xlab = "Design", ylab = "Cases Sold", main = "(b) Interval Plot")
arrows(1:4, means+diff, 1:4, means-diff, angle = 90, code = 3)

rm(bar, diff) 


################################################################################
## EXAMPLE                                                          (p 739-41) #
## Inferences For Difference Between Two Factor Level Means                    #
################################################################################
means <- by(sold, package, mean)               ## Defined from earlier
D = as.vector(means[3] - means[4])             ## (17.11)
se = sqrt(10.55 * (1/4 + 1/5))                 ## (17.14)
c("lwr" = D, "upr" = D) +                      ## (17.16)
  c(-qt(0.975,15) * se, qt(0.975,15) * se)     ## "+/- t(0.975, 15) x se"
D / se                                         ## (17.18)

rm(D, se, means)


################################################################################
## EXAMPLE                                                          (p 741-43) #
## Inferences For Contrast of Factor Level Means                               #
################################################################################
lengths <- by(data, package, nrow)
means   <- by(sold, package, mean)
ci <- c(1/2, 1/2, -1/2, -1/2)                  ## contrasts
L  <- sum(ci * means)                          ## (17.20)
se <- sqrt(10.55 * sum(ci^2 / lengths))        ## (17.22)
c(lwr = L, upr = L) +                          ## (17.24)
  c(-qt(.975, 15)*se, qt(.975, 15)*se)
L / se                                         ## (17.26)

rm(ci, L, lengths, means, se)


################################################################################
## Input Rust Inhibitor Data                                                   #
################################################################################
detach(data)
data <- read.table("data/CH17TA02.txt",
  col.names = c("value", "brand", "unit")
);
attach(data)


################################################################################
## FIGURE 17.4                                                         (p 748) #
## Paired Comparison Plot--Rust Inhibitor Example                              #
##                                                                             #
## The TukeyHSD method (stats) performs a test for pairwise differences and    #
## contains a method for plotting those differences. It is not the same as     #
## this figure, but it is just as informative.                                 #
################################################################################
plot( TukeyHSD(aov(value ~ factor(brand))) )     ## Diff 4-1 contains 0


################################################################################
## TABLE 17.3                                                          (p 749) #
## Simultaneous Confidence Intervals and Tests for Pairwise                    #
## Differences Using the Tukey Procedure--Rust Inhibitor Example               #
##                                                                             #
## TukeyHSD does a systematic pairing of the possible combinations, and they   #
## do not match those in the table. However, the content is just as meaningful #
## for our purposes. Also, it gives the adjusted p-value, not q*. For access   #
## to the q* distribution, see ?qtukey for details.                            #
##                                                                             #
## Example 2--Unequal Sample Sizes (p 751) can be used in a similar manner     #
## with the following code (assuming data was input)                           #
##   TukeyHSD(aov(sold ~ factor(package)), conf.level = .90)                   #
################################################################################
TukeyHSD(aov(value ~ factor(brand)))


################################################################################
## Input Kenton Food Data                                                      #
################################################################################
detach(data)
data <- read.table("data/CH16TA01.txt", 
  colClasses = c("numeric", "factor", "factor"),
  col.names = c("sold", "package", "store"))
attach(data)


################################################################################
## EXAMPLE                                                           (p 753-5) #
## Scheffé Multiple Comparison Procedure                                       #
##                                                                             #
## As far as we know, R has no methods for this multiple comparison procedure. #
## Furthermore, as the authors point out (p 755), if ony pairwise comparisons  #
## are being made, Tukey gives narrower confidence limits. Since there is      #
## neither pedagogical nor practical benefit to this example, it will be       #
## ignored. The calculations were performed earlier, and the F quantile can be #
## calculated even easier. I leave it to the interested reader to concern      #
## themselves with the details of this example.                                #
################################################################################


################################################################################
## EXAMPLE                                                           (p 756-7) #
## Bonferroni Multiple Comparison Procedure                                    #
################################################################################
cont <- matrix(scan(), byrow = TRUE, ncol = 4)     ## Define contrast matrix
  0.5  0.5 -0.5 -0.5     
  0.5 -0.5  0.5 -0.5     

(B  <- qt(1-0.025/(2*2), 15))                            ## (17.46a)
se <- sqrt(10.55 * apply(cont, 1,                        ## (17.22)
  function(x) sum(x^2 / tapply(sold, package, length)))
); ## end sqrt
data.frame(                                     
  L   = cont %*% tapply(sold, package, mean),
  lwr = cont %*% tapply(sold, package, mean) - B*se,     ## (17.46)
  upr = cont %*% tapply(sold, package, mean) + B*se      ## (17.46)
);
abs(cont %*% tapply(sold, package,mean) / 1.5) <= B      ## (17.47) TRUE -> H0

rm(B, cont, L, se)


################################################################################
## FIGURE 17.5                                                         (p 759) #
## Analysis of Means Plot--Kenton Food Company Example                         #
##                                                                             #
## We cannot find any ANOM methods in R, but since the computations are        #
## straightforward, both the values of interest and graphic could be designed  #
## from the following code.                                                    #
################################################################################
r = 4
u = sum(tapply(sold, package, mean)) / r                          ## (17.48a)
l = tapply(sold, package, length)                                 ## lengths ni
s = vector("numeric", r)
for(i in seq(r))
  s[i] = (10.55/l[i]) * ((r-1)/r)^2 + (10.55/r^2) * sum(1/l[-i])  ##(17.49)

plot(tapply(sold, package, mean), pch = 19, ylim = c(10, 30))
segments(seq(4), u, seq(4), tapply(sold, package, mean))
lines(seq(1, 4.5, 0.5), rep(u+s, each = 2), type = "S")
lines(seq(1, 4.5, 0.5), rep(u-s, each = 2), type = "S")
abline(h = u)

rm(i, l, r, s, u)


################################################################################
## Input Piecework Trainees Data                                               #
################################################################################
detach(data)
data <- read.table("data/CH17TA04.txt", 
  colClasses = c("numeric", "factor", "factor"),
  col.names  = c("units", "hours", "employee"))
levels(data$hours) <- c(6, 8, 10, 12)
attach(data)


################################################################################
## TABLE 17.4                                                          (p 762) #
## Data--Piecework Trainees Example                                            #
################################################################################
xtabs(units ~ hours + employee)


################################################################################
## FIGURE 17.6                                                         (p 763) #
## Computer Output--Piecework Trainees Example                                 #
##                                                                             #
## The sort of SPSS output in this figure is not reproducible in whole by R,   #
## but R does contain the necessary methods to present the given results.      #
##                                                                             #
## The homogenous subsets can be checked via TukeyHSD, but to replicate a      #
## similar summary table, you can use the multcompTs (multcompView) function.  #
## Its input is vector of values, either logical or numeric (p values), among  #
## other things. See the documentation for details. The output matrix          #
## specifies group pairs by "0" and non-groups pairs by "-1". The use of       #
## multcompLetters is to display those entities that are grouped.
##                                                                             #
## The TukeyHSD output can be read to see that all pairs are significantly     #
## different. The use of multcompTs matrix having all off-diag "-1" indicates  #
## the same thing. Finishing with multcompLetters on the same output, which    #
## automatically converts the vector of p-values into logicals based on the    #
## "compare" parameter (see help), which defaults to "< 0.05", we will see a   #
## vector of the groups. They are all individually subset. Compare that with   #
## the same operation on the help page example (i.e., "multcompLetters(dif3)").#
################################################################################
library(multcompView)
fit <- aov(units ~ hours)
model.tables(fit, "means")                      ## Table of means
summary(fit)                                    ## ANOVA
qtukey(0.95, 4, 24)                             ## q(0.95, 4, 24) = 3.90
TukeyHSD(fit)                                   
multcompTs(TukeyHSD(fit)$hours[, "p adj"])      ## Difference Truth Matrix
multcompLetters(TukeyHSD(fit)$hours[, "p adj"]) ## Homogeneous Subsets
by(units, hours, mean)[2:4] -                   ## Diminishing Marginal 
  by(units, hours, mean)[1:3]                   ## Returns. Point (2), (p 764)


################################################################################
## FIGURE 17.7                                                         (p 764) #
## Scatter plot and Fitted Quadratic Response                                  #
## Function--Piecework Trainees Example                                        #
################################################################################
fun <- function(x) -3.73571 + 9.175*x - 0.3125*x^2   ## (17.52) transformed
plot(as.numeric(as.vector(hours)), units, pch = 19,
  xlab = "Hour of Training", ylab = "Number of Acceptable Units")
curve(fun, 6, 12, add = TRUE)

rm(fun)


################################################################################
## TABLE 17.5                                                          (p 765) #
## Illustration of Data for Regression Analysis--Piecework Trainees Example    #
################################################################################
cbind(
  i   = as.numeric(hours),
  j   = employee,
  Y   = units,
  x   = as.vector(scale(as.numeric(as.vector(hours)), T, F)),
  xsq = as.vector(scale(as.numeric(as.vector(hours)), T, F))^2
);


################################################################################
## TABLE 17.6                                                          (p 765) #
## Analysis of Variance--Piecework Trainees Example                            #
################################################################################
fit <- lm(units ~ x + xsq, data = transform(data,
  x  = scale(as.numeric(as.vector(hours)), T, F),
  xsq = scale(as.numeric(as.vector(hours)), T, F)^2)
); 
anova(fit)                                      ## (a) Regression Model (17.51)
anova(lm(units ~ hours))                        ## (b) ANOVA (17.50)
fit <- anova(fit, lm(units ~ hours))            ## (c) ANOVA for Lack of Fit
fit[2,4] / (fit[2,2] / fit[2,1])                ## F* = MSLF / MSPE = 0.136
qf(0.95, fit[2, 3], fit[2,1])                   ## Decision on q(0.94, 1, 24)


################################################################################
## Clean up the R environment from this session                                #
################################################################################
detach(data)
rm(data, fit)
