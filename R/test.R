source('./MetaAnalysis.R')

set.seed(123)
p <- 8
n1 <- 20
n2 <- 15
n3 <- 30
d0 <- data.frame(group=sample(6:8,n1,replace=TRUE), matrix(rnorm(p*n1),ncol=p))
data1 <- data.frame(group=sample(1:3,n1,replace=TRUE), matrix(rnorm(p*n1),ncol=p))
data2 <- data.frame(group=sample(1:2,n2,replace=TRUE), matrix(rnorm(p*n2),ncol=p))
data3 <- data.frame(group=sample(2:4,n3,replace=TRUE), matrix(rnorm(p*n3),ncol=p))

## Main
res = meta.analysis(d0,data1, data2, data3, method='test')
