## ----set-global-opts, include=FALSE--------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)

## ------------------------------------------------------------------------
data("CH15TA01", package = "ALSM")
data("CH15FI01", package = "ALSM")
data("CH15FI06", package = "ALSM")
data("CH15FI09", package = "ALSM")

## ------------------------------------------------------------------------
.data <- CH15FI01
names(.data) <- c("y", "x")

## ----stripchart, fig.width=12, fig.height=6------------------------------
stripchart(split(.data$y, .data$x), method = "stack", main = "Attendance", xlab = "Rating",
           group.names = c("Attend", "Not Attend"), pch = 19, col = "gray30")
abline(h = 2 - 0.05)  # Add a reference line  below the 2nd group

## ------------------------------------------------------------------------
.data <- CH15FI06
names(.data) <- c("y", "x")

# Treat the factor as ordered
.data <- transform(.data, x = factor(x, ordered = TRUE, levels = c("Low", "Medium", "High", "Very High")))


# Specify appropriate contrasts to great this ordered factor as a categorical factor
fit <- lm(y ~ x, .data, contrasts = list(x = "contr.treatment"))
summary(fit)
anova(fit)

## ------------------------------------------------------------------------
plot(y ~ unclass(x), .data, ylim = c(0, 1600), pch = 19, xaxt = 'n', xlab = "Oven Temperature", ylab = "Volume")
axis(1, at = seq(4), labels = levels(.data$x))
title("Summary Plot")

## ------------------------------------------------------------------------
.data <- CH15FI09
names(.data) <- c("y", "x1", "x2")
.data <- transform(.data, x1 = factor(x1, ordered = TRUE, levels = c("Low", "Medium", "High", "Very High")))

fit <- lm(y ~ x1 + x2, .data, contrasts = list(x1 = "contr.treatment"))
summary(fit)
anova(fit)

## ------------------------------------------------------------------------
plot(y ~ unclass(x1), .data, col = unclass(x2), pch = 20, xaxt = 'n', xlab = "Oven Temperature", ylab = "Volume")
lines(y ~ unclass(x1), .data, subset = x2 == 'Plant A', col = unclass(x2))
lines(y ~ unclass(x1), .data, subset = x2 == 'Plant B', col = unclass(x2))
axis(1, at = seq(4), labels = levels(.data$x1))
title("Summary Plot")
text(3.6, 1400, "Plant B", col = 2)
text(3.5, 1100, "Plant A")

## ------------------------------------------------------------------------
.data <- CH15TA01
names(.data) <- c("subject", "control", "experiment", "difference")

.data <- transform(.data,
                x1 = stack(.data, select = c(control, experiment))$ind,
                x2 = factor(subject),
                y  = stack(.data, select = c(control, experiment))$values)

fit <- lm(y ~ x1 + x2, .data)

## ------------------------------------------------------------------------
tab <- xtabs(y ~ x2 + x1, .data)  # Wide format table like before
tab <- addmargins(tab, 1, FUN = list(list('Mean' = mean, 'Std' = sd)))  # Add sample means and std. dev. to bottom margin
tab <- addmargins(tab, 2, FUN = diff)  # Add within-subject differences to side margin
round(tab, 4)

## ------------------------------------------------------------------------
x <- jitter(unclass(.data$x1), 0.12)
plot(y ~ x, .data, pch = 19, xaxt = 'n', xlab = "Treatment", ylab = "Diameter")
axis(1, at = seq(2), labels = levels(.data$x1))
for (id in .data$x2)
  lines(y ~ x, .data, subset = x2 == id)
title("Summary Plot")

## ------------------------------------------------------------------------
summary(fit)
anova(fit)
anova(lm(y ~ x1, .data), fit)  # (15.20)

