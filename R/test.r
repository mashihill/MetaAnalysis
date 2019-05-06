#library(dplyr)
#source('./R/meta.analysis.R')
#library(MetaAnalysis)

#set.seed(123)
#p <- 100
#n1 <- 200
#n2 <- 150
#n3 <- 30
#d0 <- data.frame(group=sample(6:8,n1,replace=TRUE), matrix(rnorm(p*n1),ncol=p))
#data1 <- data.frame(group=sample(1:3,n1,replace=TRUE), matrix(rnorm(p*n1),ncol=p))
#data2 <- data.frame(group=sample(1:2,n2,replace=TRUE), matrix(rnorm(p*n2),ncol=p))
#data3 <- data.frame(group=sample(2:4,n3,replace=TRUE), matrix(rnorm(p*n3),ncol=p))
#res = meta.analysis(d0,data1, data2, data3, method='Fisher')

#browseVignettes("MetaAnalysis")

## Main
set.seed(123)
data1 <- data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(100*200),ncol=100))
data2 <- data.frame(group=sample(1:2,150,replace=TRUE), matrix(rnorm(100*150),ncol=100))
data3 <- data.frame(group=sample(1:4,400,replace=TRUE), matrix(rnorm(100*400),ncol=100))

res = meta.analysis(data1, data2, data3, method=c('maxP'))
#res = meta.analysis(data1, data2, data3, method=c('Fisher', 'Stouffer', 'minP', 'maxP'))
p.matrix = res$p.matrix
pooled.p.matrix = res$pooled.p.matrix
test.performed = res$test.performed
pooled.p.matrix
