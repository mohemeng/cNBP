#' Density for the continuous Negative Binomial distribution with parameters r and beta
#'
#' @param x vector of (non-negative) quantiles
#' @param r parameter from the hierarchical continuous gamma distribution
#' @param beta parameter from the hierarchical continuous beta distribution
#'
#' @return value of the probability density function at x
#' @export
#'
#' @examples dcnbinom(x = 1, r = 2, beta = 3)

dcnbinom <- function(x, r, beta){
  if(!is.numeric(x) | !is.numeric(r)| !is.numeric(beta)) stop("inputs x, r, and beta must be numeric")
  if(beta <= 0 | r <= 0) stop("r and beta must be positive")

  g.1 <- function(x, t){
    bb <- ((beta^r)/gamma(r))*(t^(r-1))*exp(-beta*t)*digamma(x)*pgamma(t, x)
    return(bb)}

  g.2 <- function(x, m){
    bb <- (m^(x - 1)/gamma(x))*exp(-m)*pgamma(beta*m, r, lower.tail = FALSE)*log(m)
    return(bb)}

  ar <- array(NA, dim = length(x))
  for(i in x){
    pos <- match(i, x)
    if(i > 0){
      suppressWarnings(ar[pos] <- integrate(g.1, lower = 0, upper = Inf, x = i)$value - integrate(g.2, lower = 0, upper = Inf, x = i)$value)}
    else{
      ar[pos] <- 0
    }
  }
  return(ar)}

#' Distribution function of the continuous Negative Binomial distribution with parameters r and beta
#'
#' @param x vector of (nonnegative) quantiles
#' @param r parameter from the hierarchical gamma distribution
#' @param beta parameter from the hierarchical gamma distribution
#' @param lower.tail if TRUE, returns CDF, otherwise computes the SF
#'
#' @return the SF/CDF of the continuous Negative Binomial distribution with parameters r and beta
#' @export
#'
#' @examples pcnbinom(x = 4, r = 3, beta = 5)

pcnbinom <- function(x, r, beta, lower.tail = TRUE){
  if(!is.numeric(x) | !is.numeric(r)| !is.numeric(beta)) stop("inputs x, r, and beta must be numeric")
  if(beta <= 0 | r <= 0) stop("r and beta must be positive")

  g <- function(x, t){
    bb <- ((beta^r)/gamma(r))*pgamma(t, x)*t^(r - 1)*exp(- beta*t)
    return(bb)
  }
  ar <- array(NA, dim = length(x))
  for(i in x){
    pos <- match(i, x)

    if(i >= 0){
      if(lower.tail){
        ar[pos] <- 1 - integrate(g, lower = 0, upper = Inf, x = i)$value}
      else{
        ar[pos] <- integrate(g, lower = 0, upper = Inf, x = i)$value
      }
    }
    else{ar[pos] = 0}
  }
  return(ar)
}

#' Random number generation for the continuous Negative Binomial distribution with parameters r and beta
#'
#' @param n number of observations
#' @param r parameter from the hierarchical gamma distribution
#' @param beta parameter from the hierarchical gamma distribution
#'
#' @return vector containing random values from the continuous Negative Binomial distribution with parameters r and beta
#' @export
#'
#' @examples rcnbinom(n = 100, r = 3, beta = 7)

rcnbinom <- function(n, r, beta){
  if(!is.numeric(n) | !is.numeric(r) | !is.numeric(beta)) stop("inputs n, r, and beta must be numeric")
  if(n <= 0 | n%%1 != 0) stop("input n must be a positive integer")
  if(r <= 0) stop("input r must be a positive number")
  if(beta <= 0) stop("input beta must be a positive number")

  suppressMessages(require(rootSolve))
  bb <- 1:n
  ans.roots <- array(NA, dim = length(bb))
  gen <- runif(n)  # n uniform values

  for(i in bb){
    lambda <- rgamma(1, shape = r, scale = 1/beta)
    func <- function(x){
      f <- pgamma(lambda, x, 1, lower.tail = FALSE) - gen[i]
      return(f)}

    suppressWarnings(ans.roots[i] <- uniroot(func, lower =  0.2, upper = 10, extendInt = "yes")$root)}
  return(ans.roots)
}

#' Quantile function for the continuous Negative Binomial distribution with parameters r and beta
#'
#' @param p vector of probabilities
#' @param r parameter from the hierarchical gamma distribution
#' @param beta parameter from the hierarchical gamma distribution
#'
#' @return vector of quantiles
#' @export
#'
#' @examples qcnbinom(p = 0.3, r = 2, beta = 5)

qcnbinom <- function(p, r, beta){
  if(!is.numeric(p) | !is.numeric(r) | !is.numeric(beta)) stop("inputs p, r, and beta must be numeric")
  if(any(p <= 0) | any(p > 1)) stop("input p must be in (0, 1)")
  if(r <= 0) stop("input r must be a positive number")
  if(beta <= 0) stop("input beta must be a positive number")

  suppressMessages(require(rootSolve))
  ans.roots <- array(NA, dim = length(p))
  for(i in p){
    pos <- match(i, p)
    func <- function(x, t){
      h <- function(t) (beta^r)*exp(-lgamma(r))*pgamma(t, x)*t^(r - 1)*exp(-beta*t)
      f <-  1 - integrate(h, lower = 0, upper = Inf)$value - i
      return(f)
    }
    suppressWarnings(ans.roots[pos] <- uniroot(func, lower =  0.0002, upper = 100, extendInt = "yes")$root)
  }
  return(ans.roots)
}

#' MME estimates for the parameters r and beta in the continuous Negative Binomial distribution
#'
#' @param data required to check for model fitting
#'
#' @return the estimates of r and beta for the continuous Negative Binomial distribution
#' @export
#'
#' @examples cnbinomMME(data = rcnbinom(n = 100, r = 2, beta = 3))

cnbinomMME <- function(data){
  est <- function(para){
    suppressMessages(require(nleqslv)); suppressMessages(require(cubature))
    (adaptIntegrate(function(x) (((para[2]/(1 - para[2]))^para[1])/gamma(para[1]))*pgamma(x[2], x[1])*x[2]^(para[1] - 1)*exp(- (para[2]/(1 - para[2]))*x[2]), lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))$integral - mean(data))^2 +
    (adaptIntegrate(function(x) 2*x[1]*(((para[2]/(1 - para[2]))^para[1])/gamma(para[1]))*pgamma(x[2], x[1])*x[2]^(para[1] - 1)*exp(- (para[2]/(1-para[2]))*x[2]), lowerLimit = c(0, 0), upperLimit = c(Inf, Inf))$integral - mean(data^2))^2
  }
  parastart <- c(2, .2)
  final <- optim(parastart, est)$par
  final <- list(r = final[1], beta = final[2]/(1 - final[2]))
  return(final)
}

#' Quantile-Quantile plot for Continuous Negative Binomail distribution
#'
#' @param data required for model fitting
#' @param correction tail correction value
#'
#' @return a quantile-quantile plot
#' @export
#'
#' @examples qq.cnbinom(data = rcnbinom(100, r = 2, beta = 3))

qq.cnbinom <- function(data, correction = 0.25){
  if(!is.numeric(data) | !is.numeric(correction)) stop("data and correction must contain numerical values")
  emp.cdf <- function(dat = data, cor = correction){
    s <- sort(dat); n <- length(dat)
    uni <- unique(sort(dat))
    cdf <- array(NA, dim = length(uni))
    for(i in uni){
      pos <- match(i, uni)
      cdf[pos] <- length(which(s <= i))/(n + cor)
    }
    return(data.frame(X = uni, Empirical_CDF = cdf))
  }

  jj <- emp.cdf()
  est.val <- cnbinomMME(data = data)
  mm <- qcnbinom(p = jj$Empirical_CDF, r = est.val$r, beta = est.val$beta)

  plot(mm, jj$X, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       main = "Continuous Negative Binomial Q-Q Plot")
  abline(a = 0, b = 1, col = "red")
}
