set.seed(123)
p <- 6
data1 <- data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(p*200),ncol=p))
data2 <- data.frame(group=sample(1:2,150,replace=TRUE), matrix(rnorm(p*150),ncol=p))

library(dplyr)

num.data = 2

p.matrix = matrix(NA, nrow = num.data, ncol = p)
pooled.p.matrix = data.frame(matrix(NA, nrow = 4, ncol = p))
x <- c("Fisher", "Stouffer", "minP", "maxP")
rownames(pooled.p.matrix) <- x
colnames(pooled.p.matrix) <- names(data1[-1])

fisher.pool = function(p.vals) {
  X2 = -2*sum(log(p.vals))
  
  #cat("builtin: ", sumlog(p.vals)$p, "\n")
  pooled.p = pchisq(X2, 2*length(p.vals), lower.tail = F)
  #cat("own: ", pooled.p, "\n")
  return(pooled.p)
  ## sumlog(p.vals) ## metap library
}

stouffer.pool = function(p.vals) {
  k = length(p.vals)
  #print(sumz(p.vals))
  T_stouffer = sum(qnorm(p.vals)/sqrt(k))
  #cat("T: ", T_stouffer, "\n")
  #cat("builtin: ", sumz(p.vals)$p, "\n")
  pooled.p = pnorm(T_stouffer, lower.tail = T)
  #cat("own: ", pooled.p, "\n")
  return(pooled.p)
}

min.pool = function(p.vals) {
  k = length(p.vals)
  #print(sumz(p.vals))
  T_min = min(p.vals)
  #cat("T: ", T_stouffer, "\n")
  #cat("builtin: ", sumz(p.vals)$p, "\n")
  pooled.p = pbeta(T_min, 1, k)
  #cat("own: ", pooled.p, "\n")
  return(pooled.p)
}


max.pool = function(p.vals) {
  k = length(p.vals)
  #print(sumz(p.vals))
  T_max = max(p.vals)
  #cat("T: ", T_stouffer, "\n")
  #cat("builtin: ", sumz(p.vals)$p, "\n")
  pooled.p = pbeta(T_max, k, 1)
  #cat("own: ", pooled.p, "\n")
  return(pooled.p)
}


my.aov = function(data) {
  p.value.vec = rep(NA, p)
  for (k in 1:p) {
    listids <- list()
    bioname = paste("X", k, sep = '')
    for (ids in unique(data$group)){
      subdf <- subset(x=data, subset=group==ids)
      listids[[ids]] <- shapiro.test(get(bioname, subdf))
    }
    
    isNormal = T
    for (i in 1:length(listids)) {
        if (listids[[i]]$p.value < 0.05) {
          isNormal = F
        }
    }
    
    if (isNormal == F) {
      res.aov = kruskal.test(as.formula(paste(bioname, "~ factor(group)")), data = data)
      p.val = res.aov$p.value
    } else {
      res.aov = aov(as.formula(paste(bioname, "~ factor(group)")), data = data)
      p.val = summary(res.aov)[[1]][["Pr(>F)"]][[1]]
    }
    p.value.vec[k] = p.val
  }
  return(p.value.vec)
}

p.matrix[1,] = my.aov(data1)
p.matrix[2,] = my.aov(data2)



for (i in 1:dim(p.matrix)[2]) {
  pooled.p.matrix["Fisher", i] = fisher.pool(p.matrix[,i])
  pooled.p.matrix["Stouffer", i] = stouffer.pool(p.matrix[,i])
  pooled.p.matrix["minP", i] = min.pool(p.matrix[,i])
  pooled.p.matrix["maxP", i] = max.pool(p.matrix[,i])
}

pooled.p.matrix

rbind(ppool(t(p.matrix), method = "FI")@pvalue,
      ppool(t(p.matrix), method = "SI")@pvalue,
      ppool(t(p.matrix), method = "MI")@pvalue,
      ppool(t(p.matrix), method = "MA")@pvalue)