library(dplyr)
library(JumpTest)
source('./pooling.R')
source('./my.aov.R')

meta.analysis = function(..., method) {
  
  ## TODO: method
  
  data.list = list(...)
  num.data = length(data.list)
  p = dim(data.list[[1]])[2] - 1
  
  
  ## Error handling
  if (num.data < 2) {
    stop("Number of data is less than 2.")
  }
  
  if (num.data > 5) {
    stop("Number of data is greater than 5.")
  }
  
  if (!all(sapply(data.list, function(x) {return(dim(x)[2]-1)}) == p)) {
    stop("Number of columns are different.")
  }
  
  
  ## Output initialization
  p.matrix = matrix(NA, nrow = num.data, ncol = p)
  colnames(p.matrix) = names(data.list[[1]][-1])
  pooled.p.matrix = data.frame(matrix(NA, nrow = 4, ncol = p))
  x = c("Fisher", "Stouffer", "minP", "maxP")
  rownames(pooled.p.matrix) = x
  colnames(pooled.p.matrix) = names(data.list[[1]][-1])
  test.performed = c()

  ## Calculating p-matrix
  for (i in 1:length(data.list)) {
    res = my.aov(data.list[[i]])
    p.matrix[i,] = res$p.value.vec
    test.performed = c(test.performed, res$test.performed)
  }
  
  ## Calculating pooled p-matrix
  for (i in 1:dim(p.matrix)[2]) {
    pooled.p.matrix["Fisher", i] = fisher.pool(p.matrix[,i])
    pooled.p.matrix["Stouffer", i] = stouffer.pool(p.matrix[,i])
    pooled.p.matrix["minP", i] = min.pool(p.matrix[,i])
    pooled.p.matrix["maxP", i] = max.pool(p.matrix[,i])
  }
  
  ## Getting pooled p-matrix by R package: JumpTest  
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