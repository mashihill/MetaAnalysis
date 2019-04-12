set.seed(123)
p <- 5
n1 <- 20
n2 <- 15
data1 <- data.frame(group=sample(1:3,n1,replace=TRUE), matrix(rnorm(p*n1),ncol=p))
data2 <- data.frame(group=sample(1:2,n2,replace=TRUE), matrix(rnorm(p*n2),ncol=p))

library(dplyr)
library(JumpTest)


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


test.performed <- c()

my.aov = function(data) {
  p.value.vec = rep(NA, p)
  for (k in 1:p) {
    shapiro.vec = c()
    bioname = paste("X", k, sep = '')
    for (ids in unique(data$group)){
      subdf <- subset(x=data, subset=group==ids)
      shapiro.vec <- c(shapiro.vec, shapiro.test(get(bioname, subdf))$p.value)
    }
    
    isNormal = T
    if (min(shapiro.vec) < 0.05) {
      isNormal = F
    }
    
    if (length(unique(data$group)) > 2) {
      if (isNormal == F) {
        print('Kruskal')
        test.performed = c(test.performed, 'Kruskal')
        
        res.aov = kruskal.test(as.formula(paste(bioname, "~ factor(group)")), data = data)
        p.val = res.aov$p.value
      } else {
        print('Anova')
        test.performed = c(test.performed, 'Anova')
        print(test.performed)
        res.aov = aov(as.formula(paste(bioname, "~ factor(group)")), data = data)
        p.val = summary(res.aov)[[1]][["Pr(>F)"]][[1]]
      }
      
    
    } else { ##  #group == 2
      subdf1 <- subset(x=data, subset=group==unique(data$group[1]))
      subdf2 <- subset(x=data, subset=group==unique(data$group[2]))
      
      if (isNormal == T) {
        fp <- var.test(get(bioname, subdf1), get(bioname, subdf2))$p.value
        
        if(fp>0.05) {  #equal variance
          print('Equal T')
          test.performed = c(test.performed, 'Equal T')
          p.val <- t.test(get(bioname, subdf1),get(bioname, subdf2),var.equal=T)$p.value
        } else {
          print('Unequal T')
          test.performed = c(test.performed, 'Unequal T')
          p.val <- t.test(get(bioname, subdf1),get(bioname, subdf2),var.equal=F)$p.value
        }
      
      } else {
        print('Wilcox')
        test.performed = c(test.performed, 'Wilcox')
        p.val <- wilcox.test(get(bioname, subdf1),get(bioname, subdf2))$p.value
      }

      }
    p.value.vec[k] = p.val
  }
  return(list(p.value.vec=p.value.vec, test.performed=test.performed))
  
}


## Main
p.matrix[1,] = my.aov(data1)$p.value.vec
test.performed = c(test.performed, my.aov(data1)$test.performed)
p.matrix[2,] = my.aov(data2)$p.value.vec
test.performed = c(test.performed, my.aov(data2)$test.performed)



for (i in 1:dim(p.matrix)[2]) {
  pooled.p.matrix["Fisher", i] = fisher.pool(p.matrix[,i])
  pooled.p.matrix["Stouffer", i] = stouffer.pool(p.matrix[,i])
  pooled.p.matrix["minP", i] = min.pool(p.matrix[,i])
  pooled.p.matrix["maxP", i] = max.pool(p.matrix[,i])
}


p.matrix
pooled.p.matrix

unique(test.performed)
#
#rbind(ppool(t(p.matrix), method = "FI")@pvalue,
#      ppool(t(p.matrix), method = "SI")@pvalue,
#      ppool(t(p.matrix), method = "MI")@pvalue,
#      ppool(t(p.matrix), method = "MA")@pvalue)