## ----set-global-opts, include=FALSE-----------------------------------------------------------------------------------
knitr::opts_chunk$set(cache=FALSE, tidy=FALSE)
options(width = 120)

## ---------------------------------------------------------------------------------------------------------------------
# The data set (from chapter 6) to be used in this chapter
.data <- structure(list(x1 = c(68.5, 45.2, 91.3, 47.8, 46.9, 66.1, 49.5, 
52, 48.9, 38.4, 87.9, 72.8, 88.4, 42.9, 52.5, 85.7, 41.3, 51.7, 
89.6, 82.7, 52.3), x2 = c(16.7, 16.8, 18.2, 16.3, 17.3, 18.2, 
15.9, 17.2, 16.6, 16, 18.3, 17.1, 17.4, 15.8, 17.8, 18.4, 16.5, 
16.3, 18.1, 19.1, 16), y = c(174.4, 164.4, 244.2, 154.6, 181.6, 
207.5, 152.8, 163.2, 145.4, 137.2, 241.9, 191.1, 232, 145.3, 
161.1, 209.7, 146.4, 144, 232.6, 224.1, 166.5)), .Names = c("x1", 
"x2", "y"), class = "data.frame", row.names = c(NA, -21L))

## ---------------------------------------------------------------------------------------------------------------------
(A = c(4, 7, 10))          # Vector in R is not the type of entities we want.
as.matrix(A)               # Matrix of R vector is our basic column vector.
t(A)                       # Transpose of a matrix. Not the same as R vector! A is first coerced to column matrix.
                           
fit <- lm(y ~ x1 + x2, data = .data)

(X = model.matrix(fit))          # (5.6) Note: carries an attribute "assign"
t(X)                             # (5.7)
(Y = model.frame(fit)[1])        # (5.4) Note: this returns a data.frame object
t(Y)                             # (5.5) Note: Can still transpose frame.
(e = as.matrix(residuals(fit)))  # (5.10)

## ---------------------------------------------------------------------------------------------------------------------
(A = matrix(1:6, ncol = 2))             # 3x2 matrix
(B = matrix(c(1:3, 2:4), ncol = 2))     # 3x2 matrix
A + B
A - B                                   # 3x2 (5.8)

(A = matrix(c(2, 9, 7, 3), ncol = 2))   # 2x2 
4 * A                                   # 2x2 (5.11)

(A = matrix(c(2, 4, 5, 1), ncol = 2))   # 2x2
(B = matrix(c(4, 5, 6, 8), ncol = 2))   # 2x2
A %*% B                                 # 2x2

(A = matrix(c(1, 3, 4, 0, 5, 8), ncol = 3, byrow = TRUE)) ## 2x3
(B = matrix(c(3, 5, 2)))                # 3x1
A %*% B                                 # 2x1 (= 2x3 %*% 3x1)

## ---------------------------------------------------------------------------------------------------------------------
(beta = matrix(coef(fit)))        # coefficients from model
cbind (
  X %*% beta,                     # fitted values
  fitted(fit)                     # wrapper to get fitted values from object
); 

class(Y)
class(t(Y))                       # See difference
tryCatch(t(Y) %*% Y,              # Error produced, Y not a matrix
         error = function(e) {print(e)})
t(Y) %*% as.matrix(Y)             # (5.13) SSy = sum of squares of Y
t(X) %*% X                        # (5.14)
t(X) %*% as.matrix(Y)             # (5.15)
solve(t(X) %*% X)                 # (5.24) See Inverse below

## ---------------------------------------------------------------------------------------------------------------------
t( t(X) %*% X )                        # Symmetric matrix A = t(A)
(I = diag(3))                          # Diagonal matrix of size 3
I %*% (t(X) %*% X)                     # Also defines an identity matrix
2*I                                    # Scalar matrix
(J = matrix(1, 4, 4))                  # (5.18)
matrix(1, 4) %*% t(matrix(1, 4))       # Also equals J

(A = matrix(c(1, 2, 5, 1, 2, 2, 10, 6, 3, 4, 15, 1), ncol = 4, byrow = T))
qr(A)$rank                             # QR decomposition returns rank value
5*A[,1] + 0*A[,2] - 1*A[,3] + 0*A[,4]  # equals the 0 vector)

(A = matrix(c(2, 3, 4, 1), ncol = 2))
solve(A)                               # Inverse of A
round(solve(A) %*% A)                  # Round to remove unnecessary precision
A %*% solve(A)                         # Rounding unnecssary; first term 
                                       # in product determines precision?
det(A)                                 # (5.23b) 

## ---------------------------------------------------------------------------------------------------------------------
(A = matrix(ncol = 2, sample(seq(100), 4)))  # Randomly defined 2x2 matrices 
(B = matrix(ncol = 2, sample(seq(100), 4)))  # from 0 to 100 with equal prob.
(C = matrix(ncol = 2, sample(seq(100), 4)))
(k = sample(0:100, 1))

A + B == B + A                               # (5.25)
(A + B) + C == A + (B + C)                   # (5.26)
(A %*% B) %*% C == A %*% (B %*% C)           # (5.27)
C %*% (A + B) == C %*% A + C %*% B           # (5.28)
k * (A + B) == k*A + k*B                     # (5.29)
t(t(A)) == A                                 # (5.30)
t(A + B) == t(A) + t(B)                      # (5.31)
t(A %*% B) == t(B) %*% t(A)                  # (5.32)
t(A %*% B %*% C) == t(C) %*% t(B) %*% t(A)   # (5.33)
solve(A %*% B)
solve(B) %*% solve(A)                        # (5.34) Check manually, precision causes FALSE comparisons.
solve(A %*% B %*% C)
solve(C) %*% solve(B) %*% solve(A)           # (5.33) Same here.
solve(solve(A))
A                                            # (5.36) and here.
solve(t(A))
t(solve(A))                                  # (5.37) and here. 

## ---------------------------------------------------------------------------------------------------------------------
options(scipen = 3, digits = 4)              # Ease up on scientific notation
Y = as.matrix(Y)                             # Recall Y = model.frame(fit)[1]
vcov(fit)                                    # (5.40) see methods(class = "lm")
diag(vcov(fit))                              # Variances of coefficients
cbind(Y, X, e)
c(beta = beta)                               # (5.52) R indexes by 1, not 0
cbind(X %*% beta + e, Y)                     # (5.53) 
cbind(t(X) %*% X %*% beta, t(X) %*% Y)       # (5.59)
cbind(beta,solve(t(X) %*% X) %*% t(X) %*% Y) # (5.60)
fitted(fit)                                  # (5.70)
cbind(fitted(fit), X %*% beta)               # (5.71)

(H = X %*% solve(t(X) %*% X) %*% t(X))       # (5.73a) Big! 21x21 matrix
cbind(fitted(fit), H %*% Y)                  # (5.73)
round(H %*% H  - H, 12) == 0                 # (5.74) Difference to high precision is 0. H is idempotent
cbind(hatvalues(fit), diag(H))               # Useful in later chapters.

cbind(
  e,                                         # (5.75)
  Y - X %*% beta,                            # (5.76)
  Y - H %*% Y,
 (diag(nrow(Y)) - H) %*% Y                   # (5.78)
); 

## ---------------------------------------------------------------------------------------------------------------------
n = nrow(Y)
J = matrix(1, n, n)
I = diag(n)

anova(lm(y ~ 1, .data), fit)                      # SSTO = SSE for reduced model
(SSTO = t(Y) %*% Y - (1/n) * t(Y) %*% J %*% Y)   # (5.83)

                                                 
sum(anova(fit)[1:2, "Sum Sq"])                   # SSR in R is sum of SSR of all predictors in model
anova(fit)       
(SSE = t(e) %*% e)                               # (5.84)
(SSR = t(beta) %*% t(X) %*% Y - (1/n) * t(Y) %*% J %*% Y)  # (5.85)

c(SSTO = SSTO, SSE = SSE, SSR = SSR)             # In summary
anova(lm(y ~ 1, .data), fit)                      # All included here

(SSTO = t(Y) %*% (I - (1/n) * J) %*% Y)          # (5.89a)
(SSE  = t(Y) %*% (I - H) %*% Y)                  # (5.89b)
(SSR  = t(Y) %*% (H - (1/n) * J) %*% Y)          # (5.89b)

## ---------------------------------------------------------------------------------------------------------------------
(MSE  = as.vector(SSE / (n-3)))
(se   = MSE * solve(t(X) %*% X))                 # (5.93) 
vcov(fit)                                        # Convenient!
sqrt(diag(se))                                   # Standard error for coefs.
sqrt(diag(vcov(fit)))                            # Much easier!
summary(fit)                                     # See something familiar?

(h = matrix(c(1, 34, 13)))                       # (5.95) x1 = 34, x2 = 13
(t(h) %*% beta)                                  # (5.96)
(sh = MSE * (t(h) %*% solve(t(X) %*% X) %*% h))  # (5.98)
(sh = sqrt(sh))

# Compare with the below output for $fit and $se.fit
predict(fit, newdata = data.frame(x1 = 34, x2 = 13), se.fit = TRUE)

