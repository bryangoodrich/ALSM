## ----set-global-opts, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)
options(width = 100, scipen = 7)

## -------------------------------------------------------------------------------------------------
data("CH01TA01", package = "ALSM")
data("CH02TA04", package = "ALSM")

## -------------------------------------------------------------------------------------------------
.data        <- CH01TA01
names(.data) <- c("x", "y")
fit       <- lm(y ~ x, .data)

## -------------------------------------------------------------------------------------------------
summary(fit)
anova(fit)   

tab <- c(
  'n'    = nrow(fit$model),              # 25
  'xbar' = mean(fit$model[[2]]),         # 70
  'b0'   = coef(fit)[[1]],               # 62.366
  'b1'   = coef(fit)[[2]],               # 3.570
  'SSE'  = anova(fit)["Residuals", 2],   # 54,825.460
  'MSE'  = anova(fit)["Residuals", 3],   # 2,383.716
  'SSx'  = sum(scale(.data$x, T, F)^2),  # 19,800
  'SSTO' = sum(anova(fit)[,2]))          # 307,203
  
round(tab, 2) 

## -------------------------------------------------------------------------------------------------
ci <- function(model, Xh = 0, alpha = 0.1, m = 1, interval = c("confidence", "prediction", "band")) 
{
  int   <- match.arg(interval)
  c     <- ifelse(int == "prediction", 1, 0)
  x     <- model$model[[2]]
  n     <- length(x)                                   
  df    <- summary(model)$df[2]                        
  SSx   <- sum( scale(model$model[[2]], T, F)^2 )                     
  MSE   <- anova(model)["Residuals", 3]                
  dev   <- Xh - mean(x)
  s     <- sqrt( MSE * (c/m + 1/n + dev^2 / SSx) )     
  
  if (int == "band") 
  {
    val <- sqrt( 2 * qf(1 - alpha, 2, df) )
  } else
    val <- qt(1 - alpha/2, df)

  pred  <- coef(fit)[[1]] + coef(fit)[[2]] * Xh
  lwr   <- pred - val * s
  upr   <- pred + val * s

  return(cbind(pred, lwr, upr))
}

## -------------------------------------------------------------------------------------------------
ci(fit)
confint(fit, parm = 1, level = 0.9)

## -------------------------------------------------------------------------------------------------
newdata = c(65, 100)
ci(fit, Xh = newdata)
predict(fit, data.frame(x = newdata), interval = "confidence", level = 0.9)

## -------------------------------------------------------------------------------------------------
newdata = 100
ci(fit, Xh = newdata, interval = "prediction")
predict(fit, data.frame(x = newdata), interval = "prediction", level = 0.9)

## -------------------------------------------------------------------------------------------------
ci(fit, newdata, interval = "prediction", m = 3)
predict(fit, data.frame(x = newdata), interval = "prediction", 
        level = 0.9, weights = 3)

## -------------------------------------------------------------------------------------------------
ci(fit, 100, interval = "band")  # Example case when Xh = 100

plot(y ~ x, .data, xlab = "Lot Size X", ylab = "Hours Y")

curve(ci(fit, x, interval = "band")[, "pred"], add = TRUE)
curve(ci(fit, x, interval = "band")[,  "lwr"], add = TRUE, col = "red")
curve(ci(fit, x, interval = "band")[,  "upr"], add = TRUE, col = "red")
curve(predict(fit, data.frame(x = x), int = "c")[, 2], add = TRUE, col = "blue")
curve(predict(fit, data.frame(x = x), int = "c")[, 3], add = TRUE, col = "blue")

## -------------------------------------------------------------------------------------------------
MSR   <- anova(fit)["x",         "Mean Sq"]   # 252,377.600
MSE   <- anova(fit)["Residuals", "Mean Sq"]   #   2,383.716
Fstat <- MSR / MSE                            #     105.876
fval  <- qf(0.95, 1, 23)                      #       4.279
round(c('MSR' = MSR, 'MSE' = MSE, 'F' = Fstat, 'fval' = fval), 4)

if (abs(Fstat) <= fval) {
  print("Conclude H0")                        # Fail to reject H0
} else
  print("Conclude Ha")                        # Reject H0

## -------------------------------------------------------------------------------------------------
tab <- c(
  'SSR'  = anova(fit)["x", "Mean Sq"],         # 252,377
  'SSTO' = sum( anova(fit)["Sum Sq"] ),        # 307,203
  'R2'   = summary(fit)$r.squared,             # 0.8215 (= SSR / SSTO)
  'r'    = with(.data, cor(x, y)))                # 0.9064 (= sqrt(R2) )

round(tab, 4)

## -------------------------------------------------------------------------------------------------
.data <- CH02TA04
y1 <- .data[, 1]
y2 <- .data[, 2]

## -------------------------------------------------------------------------------------------------
cbind(
  'Population'  = y1,
  'Expenditure' = y2,
  'R1'          = rank(y1),
  'R2'          = rank(y2))

# Calculate Spearman Correlation Information
c('rho' =  sum(scale(rank(y1), T, F) * scale(rank(y2), T, F)) /
           sqrt(sum(scale(rank(y1), T, F)^2) * sum(scale(rank(y2), T, F)^2)))

cor.test(y1, y2, method = "spearman")$estimate    # Spearman Rho = 0.895

c('rho' = cor(rank(y1), rank(y2)))                # Alternative R method

# t* calculation
c('tstar' = (cor(rank(y1), rank(y2)) * sqrt(length(y1) - 2)) / 
             sqrt(1 - cor(rank(y1), rank(y2))^2))

