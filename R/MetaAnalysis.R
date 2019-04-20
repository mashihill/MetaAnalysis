library(dplyr)
library(JumpTest)

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

my.aov = function(data) {
  p = dim(data)[2] - 1
  p.value.vec = rep(NA, p)
  test.performed = c()
  
  for (k in 1:p) {
    shapiro.vec = c()
    bioname = paste("X", k, sep = '')
    p.val = NA
    
    ## check normality
    for (ids in unique(data$group)) {
      subdf <- subset(x=data, subset=group==ids)
      shapiro.vec <- c(shapiro.vec, shapiro.test(get(bioname, subdf))$p.value)
    }
    
    if (length(unique(data$group)) > 2) {
      if (min(shapiro.vec) < 0.05) { ## Not Normal
        #print('Kruskal')
        test.performed = c(test.performed, 'Kruskal')
        res.aov = kruskal.test(as.formula(paste(bioname, "~ factor(group)")), data = data)
        p.val = res.aov$p.value
      } else { ## Normal
        #print('Anova')
        test.performed = c(test.performed, 'Anova')
        res.aov = aov(as.formula(paste(bioname, "~ factor(group)")), data = data)
        p.val = summary(res.aov)[[1]][["Pr(>F)"]][[1]]
      }

    } else { ##  group size == 2
      subdf1 <- subset(x=data, subset=group==unique(data$group)[1])
      subdf2 <- subset(x=data, subset=group==unique(data$group)[2])

      if (min(shapiro.vec) >= 0.05) {  ## Normal
        f.pval = var.test(get(bioname, subdf1), get(bioname, subdf2))$p.value
        
        if(f.pval > 0.05) {  ## Equal variance
          test.performed = c(test.performed, 'Equal T')
          p.val <- t.test(get(bioname, subdf1),get(bioname, subdf2),var.equal=T)$p.value
        } else {  ## Unequal variance
          test.performed = c(test.performed, 'Unequal T')
          p.val <- t.test(get(bioname, subdf1),get(bioname, subdf2),var.equal=F)$p.value
        }
      
      } else {  ## Not normal
        test.performed = c(test.performed, 'Wilcox')
        p.val <- wilcox.test(get(bioname, subdf1),get(bioname, subdf2))$p.value
      }
    }
    
    p.value.vec[k] = p.val
    #print(p.val)
  }
 
   return(list(p.value.vec=p.value.vec, test.performed=test.performed))
}

meta.analysis = function(..., method) {
  
  data.list = list(...)
  num.data = length(data.list)
  
  p.matrix = matrix(NA, nrow = num.data, ncol = p)
  print(names(data.list))
  #rownames(p.matrix) = names(data.list)
  colnames(p.matrix) = names(data.list[[1]][-1])
  pooled.p.matrix = data.frame(matrix(NA, nrow = 4, ncol = p))
  x = c("Fisher", "Stouffer", "minP", "maxP")
  rownames(pooled.p.matrix) = x
  colnames(pooled.p.matrix) = names(data.list[[1]][-1])
  test.performed = c()

  
  for (i in 1:length(data.list)) {
    res = my.aov(data.list[[i]])
    p.matrix[i,] = res$p.value.vec
    test.performed = c(test.performed, res$test.performed)
  }
  
  for (i in 1:dim(p.matrix)[2]) {
    pooled.p.matrix["Fisher", i] = fisher.pool(p.matrix[,i])
    pooled.p.matrix["Stouffer", i] = stouffer.pool(p.matrix[,i])
    pooled.p.matrix["minP", i] = min.pool(p.matrix[,i])
    pooled.p.matrix["maxP", i] = max.pool(p.matrix[,i])
  }
  
  builtin.pooled.p = rbind(ppool(t(p.matrix), method = "FI")@pvalue,
          ppool(t(p.matrix), method = "SI")@pvalue,
          ppool(t(p.matrix), method = "MI")@pvalue,
          ppool(t(p.matrix), method = "MA")@pvalue)
  rownames(builtin.pooled.p) = x

  
  return(list(p.matrix=p.matrix, 
              pooled.p.matrix=pooled.p.matrix, 
              builtin.pooled.p=builtin.pooled.p,
              test.performed=test.performed))

}