#library(JumpTest)
source('./R/pooling.R')
source('./R/my.aov.R')

#' Illustration of meta.analysis()
#'
#' @param ... Dataframe, at least 2 and up to 5 dataframes.
#' @param method String or vector of string, containing one of the following (case sensitive): \code{'Fisher', 'Stouffer', 'minP', 'maxP'}.
#' @param alpha The significance level to use within statistical test, default = 0.05
#' 
#' @return A list with components:
#' \describe{
#'   \item{\code{p.matrix}}{A p-values matrix of dimension (number of dataframes, number of biomarkers).}
#'   \item{\code{test.performed}}{A matrix of method (\code{string}) of statistical test used.}
#'   \item{\code{pooled.p.matrix}}{A matrix of pooled p-values for each biomarkers.}
#' }
#'
#' @examples
#' data1 <- data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(100*200),ncol=100))
#' data2 <- data.frame(group=sample(1:2,150,replace=TRUE), matrix(rnorm(100*150),ncol=100))
#' data3 <- data.frame(group=sample(1:4,400,replace=TRUE), matrix(rnorm(100*400),ncol=100))
#' 
#' 
#' res = meta.analysis(data1, data2, method='Fisher')
#' p.matrix = res$p.matrix
#' pooled.p.matrix = res$pooled.p.matrix
#' test.performed = res$test.performed
#'
#' res = meta.analysis(data1, data2, data3, method=c('Stouffer', 'minP', 'maxP'))
#' p.matrix = res$p.matrix
#' pooled.p.matrix = res$pooled.p.matrix
#' test.performed = res$test.performed
#'
#' @export
#' 
meta.analysis = function(..., method=c('Fisher', 'Stouffer', 'minP', 'maxP'), alpha=0.05) {

  method = c(method)
  data.list = list(...)
  num.data = length(data.list)
  p = dim(data.list[[1]])[2] - 1  # Number of biomarkers
  
  ### Error handling
  if (num.data < 2) {
    stop("Number of data is less than 2.")
  }
  
  if (num.data > 5) {
    stop("Number of data is greater than 5.")
  }
  
  if (!all(sapply(data.list, dim)[2,] == (p+1))) {
    stop("Number of columns are different.")
  }
  
  if (!all(method %in% c("Fisher", "Stouffer", "minP", "maxP"))) {
    stop("Invalid method provided.")
  }

  ### Output initialization
  # p value matrix
  p.matrix = matrix(NA, nrow = num.data, ncol = p)
  colnames(p.matrix) = names(data.list[[1]][-1])
  
  # test performed matrix
  test.performed = data.frame(p.matrix)
  
  # pooled p value matrix
  pooled.p.matrix = data.frame(matrix(NA, nrow = length(method), ncol = p))
  rownames(pooled.p.matrix) = method
  colnames(pooled.p.matrix) = names(data.list[[1]][-1])
  

  ### Calculating p-matrix
  for (i in 1:length(data.list)) {
    res = my.aov(data.list[[i]], alpha)
    p.matrix[i,] = res$p.value.vec
    test.performed[i,] = res$test.performed
  }
  
  ### Calculating pooled p-matrix
  for (i in 1:dim(p.matrix)[2]) {
    for (m in method) {
      if (m == 'Fisher') {
        pooled.p.matrix[m, i] = fisher.pool(p.matrix[,i])
      } else if (m == 'Stouffer') {
        pooled.p.matrix[m, i] = stouffer.pool(p.matrix[,i])
      } else if (m == 'minP') {
        pooled.p.matrix[m, i] = min.pool(p.matrix[,i])
      } else if (m == 'maxP') {
        pooled.p.matrix[m, i] = max.pool(p.matrix[,i])
      } else {
        stop('Unknown pooling method')
      }
    }
  }
  
  print('her')
  
  ## Getting pooled p-matrix by R package: JumpTest  
  builtin.pooled.p = rbind(ppool(t(p.matrix), method = "FI")@pvalue,
          ppool(t(p.matrix), method = "SI")@pvalue,
          ppool(t(p.matrix), method = "MI")@pvalue,
          ppool(t(p.matrix), method = "MA")@pvalue)
  #rownames(builtin.pooled.p) = x

  return(list(p.matrix=p.matrix, 
              pooled.p.matrix=pooled.p.matrix, 
              builtin.pooled.p=builtin.pooled.p,
              test.performed=test.performed))

}