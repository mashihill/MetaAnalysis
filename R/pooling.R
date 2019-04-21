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
