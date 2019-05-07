source('./R/meta.analysis.R')
#browseVignettes("MetaAnalysis")
#document()

## Main
set.seed(123)
data1 <- data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(6*200),ncol=6))
data2 <- data.frame(group=sample(1:2,150,replace=TRUE), matrix(rnorm(6*150),ncol=6))
data3 <- data.frame(group=sample(1:4,400,replace=TRUE), matrix(rnorm(6*400),ncol=6))

res = meta.analysis(data1, data2, data3, method=c('Fisher', 'Stouffer', 'minP', 'maxP'), alpha = 0.01)
p.matrix = res$p.matrix
test.performed = res$test.performed
pooled.p.matrix = res$pooled.p.matrix

p.matrix
test.performed
pooled.p.matrix


smd <- data.frame(group=sample(1:4,20,replace=TRUE), matrix(rnorm(6*20),ncol=6))
