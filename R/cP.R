#' Density for the continuous Poisson distribution with parameter lambda
#'
#' @param x vector of (non-negative) quantiles
#' @param lambda parameter value
#'
#' @return gives the density of the continuous Poisson distribution
#' @export
#'
#' @importFrom cubature adaptIntegrate
#'
#' @examples dcpois(x = 1, lambda = 1)

dcpois <- function(x, lambda){
  if(!is.numeric(x) | !is.numeric(lambda)) stop("inputs x and lambda must be numeric")
  if(lambda <= 0) stop("lambda must be positive")

  bef <- function(y){
    g <- function(x){
      exp(-(x[1] + x[2]))*(x[1]*x[2])^(y-1)*log(x[1]/x[2])
    }
    Int <- adaptIntegrate(g, lowerLimit = c(lambda, 0), upperLimit = c(Inf, lambda))
    return(Int$integral * (gamma(y)^(-2)))
  }
  suppressWarnings(pdf.val <- sapply(x, bef))
  pdf.val[which(is.na(pdf.val) == TRUE)] = 0
  pdf.val[which(pdf.val == Inf)] = 0
  return(pdf.val)
}

#' Distribution function for the continuous Poisson distribution with parameter lambda
#'
#' @param q vector of quantiles
#' @param lambda parameter value
#' @param lower.tail logical; if TRUE (default), probabilities are P(X <= x), otherwise, P(X > x)
#'
#' @return vector of CDFs (when lower.tail is TRUE), otherwise, vector of SFs
#' @export
#'
#' @examples pcpois(q = 1:5, lambda = 2)

pcpois <- function(q, lambda, lower.tail = TRUE){
  if(!is.numeric(q) | !is.numeric(lambda)) stop("inputs x and lambda must be numeric")
  if(lambda <= 0) stop("lambda must be positive")

  suppressWarnings(df <- pgamma(lambda, shape = q, lower.tail = !lower.tail))
  df[which(is.na(df) == TRUE)] = 0
  return(df)
}

#' Random number generation for the continuous Poisson distribution with parameter lambda
#'
#' @param n number of random values to generate
#' @param lambda parameter value
#'
#' @return n random values generated from the continuous Poisson distribution
#' @export
#'
#' @examples rcpois(n = 10, lambda = 5)

rcpois <- function(n, lambda){
  if(!is.numeric(n) | !is.numeric(lambda)) stop("inputs n and lambda must be numeric")
  if(n <= 0 | n%%1 != 0) stop("input n must be a positive integer")
  if(lambda <= 0) stop("input lambda must be a positive number")

  bb <- 1:n
  ans.roots <- array(NA, dim = length(bb))
  gen <- runif(n)  # n uniform values

  for(i in bb){
    func <- function(x){
      f <- (pgamma(lambda, x, 1, lower.tail = FALSE) - gen[i])^2
      return(f)
    }
    suppressWarnings(ans.roots[i] <- optim(par = 1, func)$par)
  }
  return(unlist(ans.roots))
}

#' quantile function for the continuous Poisson distribution with parameter lambda
#'
#' @param p vector of probabilities
#' @param lambda parameter value
#'
#' @return the quantile value(s) with given probabilities p
#' @export
#'
#' @examples qcpois(p = 0.3, lambda = 4)

qcpois <- function(p, lambda){
  if(!is.numeric(p) | !is.numeric(lambda)) stop("inputs p and lambda must be numeric")
  if(lambda <= 0) stop("input lambda must be a positive number")
  if(any(p <= 0) | any(p >= 1)) stop("input p should be in (0, 1)")

  bb <- 1:length(p)
  ans.roots <- array(NA, dim = length(bb))

  for(i in bb){
    func <- function(x){
      f <- (pgamma(lambda, x, 1, lower.tail = FALSE) - p[i])^2
      return(f)
    }
    suppressWarnings(ans.roots[i] <- optim(par = 0.5, func))
  }
  return(unlist(ans.roots))
}

#' MME estimate for lambda
#'
#' @param data our data needed for parameter estimation
#'
#' @return lambda estimate by method of moments
#' @export
#'
#' @examples cpoisMME(c(rcpois(n = 10, lambda = 3)))

cpoisMME <- function(data){
  if(!is.numeric(data)) stop("data and correction must contain numerical values")
  parfunc <- function(lambda){
    func <- function(x) 1 - pgamma(lambda, shape = x, lower.tail = FALSE)
    return((integrate(func, lower = 0, upper = Inf )$value - mean(data))^2)
  }
  suppressWarnings(val.lambda <- optim(par = 1, parfunc)$par)
  return(val.lambda)
}

#' Quantile-Quantile plot for continuous Poisson
#'
#' @param data required to check for model fitting
#' @param correction tail correction value
#'
#' @return a quantile-quantile plot
#' @export
#'
#' @examples qq.cpois(data = rcpois(n = 10, lambda = 2))

qq.cpois <- function(data, correction = 0.25){
  if(!is.numeric(data) | !is.numeric(correction)) stop("data and correction must contain numerical values")
  emp.cdf <- function(x = data, cor = correction){
    s <- sort(data); n <- length(data)
    uni <- unique(sort(data))
    cdf <- array(NA, dim = length(uni))
    for(i in uni){
      pos <- match(i, uni)
      cdf[pos] <- length(which(s <= i))/n
    }
    mon <- n/(n + cor)  # using correction for final (tail) value
    cdf[length(uni)] <- mon
    uni[length(uni)] <- quantile(s, probs = mon)
    return(data.frame(X = uni, Empirical_CDF = cdf))
  }
  jj <- emp.cdf( )
  lambda.est <- cpoisMME(data = data)
  mm <- qcpois(p = jj$Empirical_CDF, lambda = lambda.est)

  plot(mm, jj$X, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       main = "Continuous Poisson Q-Q Plot")
  abline(a = 0, b = 1, col = "red")
}
