#' Illustration of fisher.pool()
#'
#' 
#' @param p.vals A vector of p-values
#'
#' @return Return a pooled p-value calculated by Fisher's pooling method.
#'
#' @examples
#' x = runif(10)
#' fisher.pool(x)
#'
#' @export
#'
fisher.pool = function(p.vals) {
  X2 = -2*sum(log(p.vals))
  pooled.p = pchisq(X2, 2*length(p.vals), lower.tail = F)
  return(pooled.p)
}

#' Illustration of stouffer.pool()
#'
#' @param p.vals A vector of p-values
#'
#' @return Return a pooled p-value calculated by Stouffer's pooling method.
#'
#' @examples
#' x = runif(10)
#' stouffer.pool(x)
#'
#' @export
#' 
stouffer.pool = function(p.vals) {
  k = length(p.vals)
  T_stouffer = sum(qnorm(p.vals)/sqrt(k))
  pooled.p = pnorm(T_stouffer, lower.tail = T)
  return(pooled.p)
}

#' Illustration of minp.pool()
#' test
#' 
#' test
#' @param p.vals A vector of p-values
#'
#' @return Return a pooled p-value calculated by minP pooling method.
#'
#' @examples
#' x = runif(10)
#' minp.pool(x)
#'
#' @export
#' 
minp.pool = function(p.vals) {
  k = length(p.vals)
  T_min = min(p.vals)
  pooled.p = pbeta(T_min, 1, k)
  return(pooled.p)
}

#' Illustration of maxp.pool()
#'
#' @param p.vals A vector of p-values
#'
#' @return Return a pooled p-value calculated by maxP pooling method.
#'
#' @examples
#' x = runif(10)
#' maxp.pool(x)
#'
#' @export
#' 
maxp.pool = function(p.vals) {
  k = length(p.vals)
  T_max = max(p.vals)
  pooled.p = pbeta(T_max, k, 1)
  return(pooled.p)
}
