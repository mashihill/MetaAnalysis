#' Illustration of fisher.pool()
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param method2order method to order colors (\code{"hsv"} or \code{"cluster"})
#' @param p.vals vector of p-values
#' @param mar margin parameters; vector of length 4 (see \code{\link[graphics]{par}})
#'
#' @return Return a vector of pooled p-values calculated by Fisher's method.
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

stouffer.pool = function(p.vals) {
  k = length(p.vals)
  T_stouffer = sum(qnorm(p.vals)/sqrt(k))
  pooled.p = pnorm(T_stouffer, lower.tail = T)
  return(pooled.p)
}

#stouffer.pool2 = function(p.vals) {
#  k = length(p.vals)
#  T_stouffer = sum(qnorm(1-p.vals)/sqrt(k))
#  pooled.p = pnorm(T_stouffer, lower.tail = F)
#  return(pooled.p)
#}

min.pool = function(p.vals) {
  k = length(p.vals)
  T_min = min(p.vals)
  pooled.p = pbeta(T_min, 1, k)
  return(pooled.p)
}

max.pool = function(p.vals) {
  k = length(p.vals)
  T_max = max(p.vals)
  pooled.p = pbeta(T_max, k, 1)
  return(pooled.p)
}
