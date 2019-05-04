my.aov = function(data) {
  p = dim(data)[2] - 1
  p.value.vec = rep(NA, p)
  test.performed = c()
  p.alpha = 0.05
  
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
      if (min(shapiro.vec) < p.alpha) { ## Not Normal
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
      
      if (min(shapiro.vec) < p.alpha) {  ## Not Normal
        test.performed = c(test.performed, 'Wilcox')
        p.val = wilcox.test(get(bioname, subdf1),get(bioname, subdf2))$p.value
        
      } else {  ## Normal
        f.pval = var.test(get(bioname, subdf1), get(bioname, subdf2))$p.value
        
        if(f.pval > 0.05) {  ## Equal variance
          test.performed = c(test.performed, 'Equal T')
          p.val <- t.test(get(bioname, subdf1),get(bioname, subdf2),var.equal=T)$p.value
        } else {  ## Unequal variance
          test.performed = c(test.performed, 'Unequal T')
          p.val <- t.test(get(bioname, subdf1),get(bioname, subdf2),var.equal=F)$p.value
        }
      }
    }
    p.value.vec[k] = p.val
  }
  
  return(list(p.value.vec=p.value.vec, test.performed=test.performed))
}
