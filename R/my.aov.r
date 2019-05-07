#' Illustration of my.aov()
#'
#' @param data Dataframe with the first column as group, other columns are biomarkers
#' @param alpha The significance level to use within statistical test, default = 0.05
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{p.value.vec}}{A vector of p-values for each biomarkers.}
#'   \item{\code{test.performed}}{A vector of method (\code{string}) of statistical test used.}
#' }
#'
#' @examples
#' data <- data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(6*200),ncol=6))
#' res = my.aov(data.list, 0.1)
#' res$p.value.vec
#'
#' @export
#' 
my.aov = function(data, alpha=0.05) {
  p = dim(data)[2] - 1  # Number of biomarkers
  p.value.vec = rep(NA, p)
  p.alpha = alpha
  test.performed = c()
  
  for (k in 1:p) {
    shapiro.vec = c()
    p.val = NA
    
    ## check normality
    for (ids in unique(data$group)) {
      subdf <- subset(x=data, subset=group==ids)
      shapiro.vec <- c(shapiro.vec, shapiro.test(subdf[,k+1])$p.value)
    }
    
    if (length(unique(data$group)) > 2) {
      if (min(shapiro.vec) < p.alpha) { ## Not Normal
        test.performed = c(test.performed, 'Kruskal')
        res.aov = kruskal.test(data[,k+1]~factor(data[,1]))
        p.val = res.aov$p.value
      } else { ## Normal
        test.performed = c(test.performed, 'Anova')
        res.aov = aov(data[,k+1]~factor(data[,1]))
        p.val = summary(res.aov)[[1]][["Pr(>F)"]][[1]]
      }
      
    } else { ##  group size == 2
      subdf1 <- subset(x=data, subset=group==unique(data$group)[1])
      subdf2 <- subset(x=data, subset=group==unique(data$group)[2])
      
      if (min(shapiro.vec) < p.alpha) {  ## Not Normal
        test.performed = c(test.performed, 'Wilcox')
        p.val = wilcox.test(subdf1[,k+1],subdf2[,k+1])$p.value
        
      } else {  ## Normal
        f.pval = var.test(subdf1[,k+1],subdf2[,k+1])$p.value
        
        if(f.pval > 0.05) {  ## Equal variance
          test.performed = c(test.performed, 'T(eq)')
          p.val <- t.test(subdf1[,k+1], subdf2[,k+1], var.equal=T)$p.value
        } else {  ## Unequal variance
          test.performed = c(test.performed, 'T(Uneq)')
          p.val <- t.test(subdf1[,k+1], subdf2[,k+1], var.equal=F)$p.value
        }
      }
    }
    p.value.vec[k] = p.val
  }
  
  return(list(p.value.vec=p.value.vec, test.performed=test.performed))
}
