#library(dplyr)
source('./R/meta.analysis.R')
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
#document()

## Main
set.seed(123)
data1 <- data.frame(group=sample(1:3,200,replace=TRUE), matrix(rnorm(6*200),ncol=6))
data2 <- data.frame(group=sample(1:2,150,replace=TRUE), matrix(rnorm(6*150),ncol=6))
data3 <- data.frame(group=sample(1:4,400,replace=TRUE), matrix(rnorm(6*400),ncol=6))

res = meta.analysis(data1, data2, data3, method=c('Fisher', 'Stouffer', 'minP', 'maxP'))
p.matrix = res$p.matrix
test.performed = res$test.performed
pooled.p.matrix = res$pooled.p.matrix

p.matrix
test.performed
pooled.p.matrix

