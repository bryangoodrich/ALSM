################################################################################
## Input Teaching Effectiveness Data                                           #
################################################################################
data <- read.table("data/CH15FI01.txt", col.names = c("rating", "attend"))
attach(data)


################################################################################
## FIGURE 15.1                                                         (p 645) #
## Teaching Performance Comparison--Teaching Effectiveness Example             #
################################################################################
stripchart(split(rating, attend), method = "stack", pch = 16, xlab = "Rating")


################################################################################
## Input Quick Bread Volume Data                                               #
################################################################################
detach(data)
data <- read.table("data/CH15FI06.txt", col.names = c("volume", "temperature"))
data <- transform(data,
  temperature = factor(temperature, ordered = TRUE, 
    levels = c("Low", "Medium", "High", "Very_High")
  )
); 
attach(data)

  ## contr.treatment needs to be specified since temperature is ordered factor
fit <- lm(volume ~ temperature, 
  contrasts = list(temperature = "contr.treatment"))
summary(fit); anova(fit)


################################################################################
## FIGURE 15.6                                                         (p 659) #
## Summary Plot--Quick Bread Volume Example                                    #
##                                                                             #
## R defaults to plotting boxes for factors (try plot(volume ~ temperature)).  #
## One can work around this with a bit of work, but it is much more convenient #
## to just make use of xyplot (lattice).                                       #
################################################################################
library(lattice)
xyplot(volume ~ temperature, ylim = c(0, 1500), pch = 19)


################################################################################
## Input Blocked Quick Bread Volume Data                                       #
################################################################################
detach(data)
data <- read.table("data/CH15FI09.txt", 
  col.names = c("volume", "temperature", "plant"))
data <- transform(data,
  temperature = factor(temperature, ordered = TRUE,
    levels = c("Low", "Medium", "High", "Very_High")
  )
);
attach(data)

fit <- lm(volume ~ temperature + plant, 
  contrasts = list(temperature = "contr.treatment"))
summary(fit); anova(fit)


################################################################################
## FIGURE 15.9                                                         (p 662) #
## Summary Plot--Blocked Quick Bread Volume Optimization Example               #
################################################################################
library(lattice)
xyplot(volume ~ temperature, groups = plant,
  type = "b", pch = 19, auto.key = list(columns = 2))


################################################################################
## Input Skin Sensitivity Experiment Data                                      #
##                                                                             #
## The data is currently in a "wide" format and possibly can be reshaped using #
## the reshape (stats) function. However, the stack (utils) function will      #
## suffice. It returns the values with an indicator variable (ind) based on    #
## the selected columns to stack; in this case, ind will consist of factors    #
## for 'control' and 'experiment'. Thus, the model will just require this long #
## dataset combined with the subject numbers.                                  #
################################################################################
detach(data)
data <- read.table("data/CH15TA01.txt",
  col.names = c("subject", "control", "experiment", "difference"))
data <- transform(data, subject = factor(subject))
attach(data)

fit <- lm(values ~ ind + subject, 
  data = cbind(subject, stack(data, select = c(control, experiment)))
); ## end model fit


################################################################################
## TABLE 15.1                                                          (p 670) #
## Data and Descriptive Statistics--Skin Sensitivity Experiment                #
################################################################################
rbind(
  data, 
  Mean = c(NA, apply(data[-1], 2, mean)),
  SD = c(NA, round(apply(data[-1], 2, sd), 2))
);


################################################################################
## FIGURE 15.13                                                        (p 670) #
## Summary Plot--Skin Sensitivity Example                                      #
################################################################################
library(lattice)
xyplot(values ~ ind, data = model.frame(fit),
  group = subject, type = "b", pch = 19, auto.key = list(columns = 5))


################################################################################
## FIGURE 15.14                                                        (p 671) #
## Regression Results--Skin Sensitivity Experiment                             #
################################################################################
summary(fit); anova(fit)
anova(lm(values ~ ind, data = model.frame(fit)), fit)     ## (15.20)


################################################################################
## Clean up the R environment from this session                                #
################################################################################
detach(data)
rm(data, fit)


