#' Illustration of fisher.pool()
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param p.vals vector of p-values
#'
#' @return Return a vector of pooled p-values calculated by Fisher's pooling method.
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
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param p.vals vector of p-values
#'
#' @return Return a vector of pooled p-values calculated by Fisher's pooling method.
#'
#' @examples
#' x = runif(10)
#' fisher.pool(x)
#'
#' @export
#' 
stouffer.pool = function(p.vals) {
  k = length(p.vals)
  T_stouffer = sum(qnorm(p.vals)/sqrt(k))
  pooled.p = pnorm(T_stouffer, lower.tail = T)
  return(pooled.p)
}

#' Illustration of min.pool()
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param p.vals vector of p-values
#'
#' @return Return a vector of pooled p-values calculated by Fisher's pooling method.
#'
#' @examples
#' x = runif(10)
#' fisher.pool(x)
#'
#' @export
#' 
min.pool = function(p.vals) {
  k = length(p.vals)
  T_min = min(p.vals)
  pooled.p = pbeta(T_min, 1, k)
  return(pooled.p)
}

#' Illustration of max.pool()
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param p.vals vector of p-values
#'
#' @return Return a vector of pooled p-values calculated by Fisher's pooling method.
#'
#' @examples
#' x = runif(10)
#' fisher.pool(x)
#'
#' @export
#' 
max.pool = function(p.vals) {
  k = length(p.vals)
  T_max = max(p.vals)
  pooled.p = pbeta(T_max, k, 1)
  return(pooled.p)
}
