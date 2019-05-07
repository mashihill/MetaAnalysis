#source('./R/meta.analysis.R')
## Main
#set.seed(123)
p=6
data1 <- data.frame(group=sample(1:3,20,replace=TRUE), matrix(rnorm(p*20),ncol=p))
data2 <- data.frame(group=sample(1:2,15,replace=TRUE), matrix(rnorm(p*15),ncol=p))
data3 <- data.frame(group=sample(1:2,40,replace=TRUE), matrix(rnorm(p*40),ncol=p))
data4 <- data.frame(group=sample(1:4,40,replace=TRUE), matrix(rnorm(p*40),ncol=p))
data5 <- data.frame(group=sample(1:3,40,replace=TRUE), matrix(rnorm(p*40),ncol=p))

res = meta.analysis(data1, data2, data3, data4, data5,
                    method=c('Fisher',  'maxP', 'Stouffer', 'minP'), alpha = 0.01)
p.matrix = res$p.matrix
test.performed = res$test.performed
pooled.p.matrix = res$pooled.p.matrix

p.matrix
test.performed
pooled.p.matrix

matplot(t(pooled.p.matrix), type='l', col = c('black', 'green', 'red', 'blue'), lty = 'solid')
legend(0,0,inset=c(-10,-10), legend = rownames(pooled.p.matrix), xpd = T, col = c('black', 'green', 'red', 'blue'), lty = 'solid')
