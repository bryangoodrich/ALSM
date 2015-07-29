# Load the data sets
env <- new.env()
load("data.rda", envir = env)



# Input the Toluca Company Data
df        <- get("CH01TA01", envir = env)
names(df) <- c("x", "y")



# TABLE 1.1     (p 19)
# Data on Lot Size and Work Hours and Needed Calculations for Least Squares
# Estimates--Toluca Company Example 
xx <- apply(df, 2, scale, TRUE, FALSE)[, 'x']
yy <- apply(df, 2, scale, TRUE, FALSE)[, 'y']

tab <- transform(df,
  "x"      = x,
  "y"      = y,
  "xdif"   = xx,
  "ydif"   = yy,
  "crp"    = round(xx * yy),
  "sqdevx" = xx^2,
  "sqdevy" = round(yy^2)
); ## End Table

tab <- rbind(
  tab,
  Total = round(apply(tab, 2, sum), 0),
  Mean  = c(colMeans(df), rep("", 5))
);

print(tab, quote = FALSE)
rm(tab, xx, yy)



# FIGURE 1.10 and 1.11     (p 20)
# Scatter Plot and Fitted Regression Line--Toluca Company Example
# Regression Output--Toluca Company Example
fit <- lm(y ~ x, data = df)
summary(fit)

par(mfrow = c(2, 2), pch = 20)
plot(y ~ x, df,  xlab = "Lot Size", ylab = "Hours")
title("(a) Scatter Plot")

# notice the different ways by which to issue plot commands
with(df, plot(x, y, xlab = "Lot Size", ylab = "Hours"))
title("(b) Fitted Regression Line")
abline(fit)



# TABLE 1.2     (p 22)
# Fitted Values, Residuals, and Squared Residuals--Toluca Company Example
#
# Notice the use of "accessor" methods fitted and resid. They are the
# recommended way of obtaining fitted values and residuals instead of direct
# named list element access such as model[["residuals"]]
tab <- cbind(
  "Lot Size (X)"       = df$x,
  "Work Hours (Y)"     = df$y,
  "Est. Mean Response" = round(fitted(fit),   2),
  "Residuals"          = round( resid(fit),   2),
  "Sq.Residuals"       = round( resid(fit)^2, 1)
); ## End Table
rbind(tab, Totals = colSums(tab))



# Clean up the R environment from this session
rm(df, fit, tab)
dev.off()

rm(env)  # Ignore this if you're moving on to the next chapter.
