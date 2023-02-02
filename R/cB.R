#'Density for the continuous binomial distribution with parameters size and probability (prob)
#'
#' @param x vector of quantiles
#' @param size number of trials
#' @param prob probability of success on each trial
#'
#' @return density value(s) of the vector x
#' @export
#'
#' @importFrom cubature adaptIntegrate
#' @examples dcbinom(x = 2, size = 5, prob = 0.5)

dcbinom <- function(x, size, prob){
  if(!is.numeric(x) | !is.numeric(size) | !is.numeric(prob)) stop("x, size, and prob must be numeric values")
  if(size <= 0) stop("size must be a positive integer")
  if(prob <= 0 | prob >= 1) stop("prob must be in the interval (0, 1)")

  bef <- function(y){
    g <- function(x){
      (x[1]*x[2])^(y - 1)*((1 - x[1])*(1 - x[2]))^(size - y)*(log(x[1]*(1 - x[2])) - log(x[2]*(1 - x[1])))
    }
    Int <- adaptIntegrate(g, lowerLimit = c(prob, 0), upperLimit = c(1, prob))
    return(Int$integral)
  }
  suppressWarnings(pdf.val <- sapply(x, bef))
  suppressWarnings(pdf.val <- pdf.val/((beta(x, size + 1 - x))^2))
  pdf.val[which(is.na(pdf.val) == TRUE)] = 0
  pdf.val[which(pdf.val == Inf)] = 0
  return(pdf.val)
}

#' Distribution function for the continuous binomial distribution
#'
#' @param x vector of quantiles
#' @param size number of trials
#' @param prob probability of success on each trial
#'
#' @return the CDF/SF values
#' @export
#'
#' @examples pcbinom(x = 4, size = 8, prob = 0.3)

pcbinom <- function(x, size, prob, lower.tail = TRUE){
  if(!is.numeric(x) | !is.numeric(size) | !is.numeric(prob)) stop("x, size, and prob must be numeric values")
  if(size <= 0) stop("size must be a positive integer")
  if(prob <= 0 | prob >= 1) stop("prob must be in the interval (0, 1)")

  suppressWarnings(int <- 1 - pbeta(prob, x, size - x + 1, lower.tail = lower.tail))

  int[which(is.na(int) == TRUE)] = 0

  int[which(x >= size + 1)] = 1
  return(int)
}

#' Random number generation for the continuous binomial distribution
#'
#' @param n number of observations
#' @param size number of trials
#' @param prob probability of success
#'
#' @return vector of random values from the continuous Binomial distribution
#' @export
#'
#' @examples rcbinom(n = 1, size = 5, prob = 0.3)

rcbinom <- function(n, size, prob){
  if(!is.numeric(size) | !is.numeric(prob) | !is.numeric(n)) stop("inputs n, size and prob must be numeric")
  if(n <= 0 | n%%1 != 0 ) stop("input n must be a positive integer greater than 1")
  if(size < 1) stop("size >= 1")

  bb <- 1:n
  ans.roots <- array(NA, dim = length(bb))
  gen <- runif(n)  # n uniform values

  for(i in bb){
    func <- function(x){
      f <- (1 - pbeta(prob, x, size - x + 1) - gen[i])^2
      return(f)
    }
    suppressWarnings(ans.roots[i] <- optim(par = size - 0.5, func))
  }
  return(unlist(ans.roots))
}

#' Quantile function for the continuous Binomial distribution
#'
#' @param p vector of probabilities
#' @param size number of trials
#' @param prob probability of success on each trial
#'
#' @return vector of quantiles
#' @export
#'
#' @examples qcbinom(p = 0.3, size = 10, prob = 0.5)

qcbinom <- function(p, size, prob){
  if(!is.numeric(size) | !is.numeric(prob) | !is.numeric(p)) stop("inputs n, size and prob must be numeric")
  if(prob <= 0 | prob >= 1) stop("prob must be in the interval (0, 1)")
  if(any(p <= 0) | any(p >= 1)) stop("prob must be in the interval (0, 1)")

  gen <- p
  bb <- 1:length(p)
  ans.roots <- array(NA, dim = length(bb))

  for(i in bb){
    func <- function(x){
      f <- (1 - pbeta(prob, x, size - x + 1) - gen[i])^2
      return(f)
    }
    suppressWarnings(ans.roots[i] <- optim(par = size - 0.5, func))
  }
  return(unlist(ans.roots))
}

#' MME estimates for r and beta of the Continuous Binomial distribution
#'
#' @param data required for parameter estimation
#'
#' @return size and probability estimates of the continuous Binomial distribution by MME
#' @export
#'
#' @examples cbinomMME(data = rcbinom(n = 1000, size = 5, prob = .5))

cbinomMME <- function(data){
  if(!is.numeric(data)) stop("data and correction must contain numerical values")
  do <- function(t){
    (integrate(function(x) pbeta(t[1], x, t[2] - x + 1), lower = 0, upper = t[2] + 1)$value - mean(data))^2 +
    (integrate(function(x) 2*x*pbeta(t[1], x, t[2] - x + 1), lower = 0, upper = t[2] + 1)$value - mean(data^2))^2
  }
  sec <- 1
  ini <- optim(c(0.3, sec), do, lower = c(1e-05, 1), upper = c(0.9999, Inf), method = "L-BFGS-B")$par
#  while(ini[1] > 1){
#  sec <- sec + 1
#  ini <- optim(c(0.3, sec), do)$par
#  }
  return(ini)
}

#' Quantile-Quantile plot for continuous Binomial distribution
#'
#' @param data required to check for model fitting
#' @param correction tail correction value
#'
#' @return a quantile-quantile plot
#' @export
#'
#' @examples qq.cbinom(data = rcbinom(n = 25, size = 5, prob = .4))

qq.cbinom <- function(data, correction = 0.25){
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
  para.est <- cbinomMME(data = data)
  mm <- qcbinom(p = jj$Empirical_CDF, size = para.est[2], prob = para.est[1])
  plot(mm, jj$X, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       main = "Continuous Binomial Q-Q Plot")
  abline(a = 0, b = 1, col = "red")
}
