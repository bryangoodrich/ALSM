#' Scheffe Multiple Comparison Procedure
#'
#' Computes post-hoc adjustment for multiple comparisons
#'
#' @param formula an object of class formula; a symbolic description of the model to be fitted.
#' @param data a data frame
#' @param cont a contrast matrix
#' @param conf.level confidence level
#' @param MSE value of mean squared error
#' @export
#' @rdname Scheffe
Scheffe <- function(formula, data, cont, conf.level = 0.90, MSE)
{
    f  <- function(contrast, mean)   {contrast %*% mean}                # (17.41)
    h  <- function(contrast, n, MSE) {sqrt(MSE * sum(contrast^2 / n))}  # (17.42)

    .data <- model.frame(formula, data)
    y  <- .data[, 1]
    x  <- .data[, 2]
    r  <- nrow(cont)

    means <- tapply(y, x, mean)
    L     <- apply(cont, 2, f, mean = means)
    se    <- apply(cont, 2, h, n = tapply(y, x, length), MSE = MSE)
    S     <- sqrt((r-1) * qf(conf.level, r-1, nrow(.data)-r))               # (17.43a)
    L + S*cbind('lwr' = -se, 'upr' = se)                                 # (17.43)
}
