# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input the Dwaine Studios Data
df <- get("CH06FI05", envir = env)
names(df) <- c("x1", "x2", "y")
fit <- lm(y ~ x1 + x2, df)



# FIGURE 6.4     (p 232)
# Scatter Plot Matrix and Correlation Matrix--Dwaine Studios Example
with(df, pairs(data.frame("SALES" = y, "TARGETPOP" = x1, "DISPOINC" = x2)))
with(df, cor  (data.frame("SALES" = y, "TARGETPOP" = x1, "DISPOINC" = x2)))



# FIGURE 6.5     (p 237)
# Multiple Regression Output and Basic Data--Dwaine Studios Example
cbind(
  "TARGETPOP (X1)" = df$x1,
  "DISPOINC (X2)"  = df$x2,
  "SALES (Y)"      = df$y,
  "FITTED"         = round(fitted(fit),2),
  "RESIDUAL"       = round(resid(fit), 2)
)
summary(fit)
anova(fit)



# FIGURE 6.6     (p 238)
# Plot of Point Clound before and after Spinning--Dwaine Studios Example
#
# R 3D plotting function requires library "lattice" to be loaded. There are
# other libaries such as "scatterplot3d" but lattice is popular and installed
# with R as a basic library to be used.
#
# To rotate the cloud we use the "screen" parameter that takes a list of the
# position angles of the coordinate axes. The default screen value is
#
#   screen = list(z = 40, x = -60, y = 0)
#
# So by eyeballing what the author has done, the choice is made to rotate
# along the z axis to z = -30. The "scales" parameter deals with the axes
# values and most commonly the "arrows = FALSE" option is utilized so that
# axis numeric values are listed instead of the default pointed arrows.
#
# Another approach used below is to make use of the "scatterplot3d", the
# "rgl" or the "Rcmdr" libraries. The latter two provide ways for drawing
# regression planes (see next example) as well as the ablity to rotate a plot
# by use of the mouse. These will have to be downloaded. Note that Rcmdr will
# require downloading a lot of other dependent packages for it to function.
# It runs as a separate GUI application for data processing.
#
# Note, scatterplot3d has advantages with regard to easy implementation,
# control, multiple graph plotting (par parameters) and additional plotting
# such as planes. For this example, though, it fails to be good for issuing
# rotations. That is better obsereved dynamically anyway for which rgl is
# much better suited. Rcmdr automatically plots a regression surface, too.
library(lattice)
cloud(y ~ x1 + x2, df, scales = list(arrows = FALSE), 
      xlab = "TARGETPOP", ylab = "DISPOINC", zlab = "SALES",
      main = "(a) Before Spinning" )

cloud(y ~ x1 + x2, df, scales = list(arrows = FALSE),
      xlab = "TARGETPOP", ylab = "DISPOINC", zlab = "SALES",
      main = "(b) After Spinning", screen = list(z = -40, x = -60, y = 0) )

library(scatterplot3d)
with(df, scatterplot3d(x1, x2, y))

library(rgl)
with(df, plot3d(x1, x2, y))

library(Rcmdr)
with(df, scatter3d(x1, y, x2))



# FIGURE 6.7     (p 240)
# Plot of Estimated Regression Surface--Dwaine Studios Example
#
# Per the above explication and examples we can either use the Rcmdr scatter
# plot that issues the regression surface automatically, or we can use the
# scatterplot3d to make a static plot issue its plane3d command. Since the
# former is a straightforward command done above, I will do the latter here.
library(scatterplot3d)
pplot <- with(df, 
  scatterplot3d(x2, x1, y,
                xlab = "DISPOINC",
                ylab = "TARGETPOP",
                zlab = "SALES")
);  # end with
pplot$plane3d(lm(y ~ x2 + x1, df))

rm(pplot)



# FIGURE 6.8     (p 242)
# Diagnostic Plots--Dwaine Studios Example
par(mfrow = c(2, 2), pch = 19)
plot(resid(fit) ~ fitted(fit), xlab = "Fitted", ylab = "Residual")
title("(a) Residual Plot against Y")

plot(resid(fit) ~ x1, df, xlab = "Targtpop", ylab = "Residual")
title("(b) Residual Plot against X1")

plot(resid(fit) ~ x2, df, xlab = "Dispoinc", ylab = "Residual")
title("(c) Residual Plot against X2")

plot(resid(fit) ~ I(x1 * x2), df, xlab = "X1X2", ylab = "Residual")
title("(d) Residual Plot against X1X2")



# FIGURE 6.9     (p 243)
# Additional Diagnotic Plots--Dwaine Studios Example
par(mfrow = c(2, 2), pch = 19)
plot(abs(resid(fit)) ~ fitted(fit), xlab = "Fitted", ylab = "Absresid")
title("(a) Plot of Absolute Residuals against Y")
qqnorm(resid(fit), main = "", xlab = "Expected", ylab = "Residual")
title("(b) Normal Probability Plot")
qqline(resid(fit))



# Estimation of Regression Parameters     (p 245)
# See Ch. 4 walk-through for why the default confint level = 0.95 is alright.
confint(fit)



# Estimation of Mean Response     (p 246)
predict(fit, data.frame(x1 = 65.4, x2 = 17.6), interval = "confidence")



# Prediction Limits for New Observations     (p 247)
# Scheffe and Bonferroni Procedures
#
# The ci.sim from Ch. 4 handles this multivariate case without change.
ci.sim <- function(model, newdata, type = c("B", "S"), alpha = 0.05) {
  g  <- nrow(newdata)
  CI <- predict(model, newdata, se.fit = TRUE)
  M  <- ifelse(match.arg(type) == "B",
          qt(1 - alpha / (2*g), model$df),              # B = (4.9a)
          sqrt( g * qf( 1 - alpha, g, model$df) )       # S = (4.8a)
        );  # end ifelse
  spred <- sqrt( CI$residual.scale^2 + (CI$se.fit)^2 )  #     (2.38) 
  x <- data.frame(
    "x"     = newdata,
    "spred" = spred,
    "fit"   = CI$fit,
    "lower" = CI$fit - M * spred,
    "upper" = CI$fit + M * spred
  );  # return value
  
  return(x)
} # end sim function

newdata <- data.frame( x1 = c(65.4, 53.1), x2 = c(17.6, 17.7) )
ci.sim(fit, newdata, "B", 0.1)  # Bonferroni Prediction
ci.sim(fit, newdata, "S", 0.1)  # Scheffe Prediction

rm(ci.sim, newdata)



# Clean up the R environment from this session
rm(df, fit)

rm(env)  # Ignore this if you're moving on to the next chapter.
