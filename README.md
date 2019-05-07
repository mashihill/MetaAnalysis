MetaAnalysis
================
Yen-An Chen
2019-05-07

Introduction
------------

The purpose of this package is to take a series of data sets containing *p* biomarkers as represented in the columns of the data frame, and calculate a series of p-values for each column, as well as run a pooling method based on user input.

Installing
----------

There are two options to install MetaAnalysis package from source,

1.  Open terminal application and enter command:

    ``` bash
    R CMD INSTALL [PATH_TO_FILE]/MetaAnalysis_0.1.tar.gz
    ```

2.  In RStudio, enter command:

``` r
install.packages("[PATH_TO_FOLDER]/MetaAnalysis_0.1.tar", repos = NULL, type="source")
```

Workflow
--------

The package will not take less than two data frames and no more than five. In addition, the package will not run unless the number of columns between all the data frames is the same. The amount rows in each data set do not have to be the same, however, each row must be organized into a group in a column placed before all the biomarkers.

For each biomarker, the package will check the group difference (and gets a p-value). The statistical test this package choose to use will be described in next section.

The second function of this package is to take the collection of p-values and pool them by biomarker. Each biomarker will undergo pooling based on four tests (As highlighted in the “Pooling Methods” Section), the Fisher, Stouffer, Minimum p-value (minP), and Maximum p-value (maxP). If one of some of these tests is specified in the “method” parameter, the package will then provide a matrix of p-values based on these tests.

Statistical Tests
-----------------

-   If number of groups in dataframe is greater than 2, we will first check its normality, if yes we then perform one way ANOVA, otherwise we use Kruskal Wallis test

-   If number of groups in dataframe equal to 2, we will also first check its normality, if yes then we check the equal variance properties, and do the T-test accordingly, if the data is not normal we will perform Wilcox rank sum test.

As the final result, user can see the statistical test performed, one way anova is denoted as `Anova`, Kruskal Wallis test as `Kruskal`, equal variance T test as `T(eq)`, unequal variance T test as `T(uneq)`, and Wilcox rank sum test as `Wilcox`.

Pooling Methods
---------------

1.  Fisher's Method The Fisher's method sums up the log-transformed p-values obtained from individual studies. Fisher's statistic $\\chi\_{Fisher}^{2} = -2 \\sum\_{i=1}^K \\log(p\_i)$ follows a *χ*<sup>2</sup> distribution with 2*K* degrees of freedom under the null hypothesis. Smaller p-values contribute larger scores to the Fisher's statistic.

2.  Stouffer's Method Stouffer's methods sums the inverse normal trasformed p-values. Stouffer's statistic $T = \\sum\_{i=1}^K \\frac{\\Phi^{-1}(p\_i)}{\\sqrt{K}}$ ( *Φ* is the cdf of a standard normal distribution) follows a standard normal distribution under the null hypothesis. Similar to Fisher's method, smaller p-values contribute more to the Stouffer's score.

3.  Minimum p-value Method The minimun p-value method take the minimum p-value among the *K* studies as the test statistic. It follows a beta distribution with degrees of freedom *α* = 1 and *β* = *K* under the null hypothesis.

4.  Maximum p-value Method The maximum p-value method take the maximum p-value as the test statistic. It follows a beta distribution with degrees of freedom *α* = *K* and *β* = 1 under the null hypothesis.

You can use: `meta.analysis(...,method)` to compute pooled p value for different biomarkers by using four different methods: Fisher, Stouffer, Maxium p value and Minium P value.

`method` can be a single string or a vector of the collection of the following options:

**Pooling Commands**

-   Fisher: set method=`'Fisher'`
-   Stouffer: set method=`'Stouffer'`
-   MinP: set method=`'minP'`
-   MaxP: set method=`'maxP'`

*Please note that setting the type is case sensitive, and will default to all methods, i.e. `c('Fisher', 'Stouffer', 'minP', 'maxP')`*

Example
-------

1.  Input dataframes
    -   The number of input dataframes should between 2 and 5.
    -   Each dataframe should have the same number of columns.
    -   The function `mata.analysis()` will automatically check if all inputs meet requirements.
    -   Any row with missing values will be dropped.
    -   If not specified, the default level of confidence is *α* = 0.05.

Consider a scenario: We have collected data from 3 different institutions, each data belongs to a certain group, and have 5 numerical values of it's different biomarkers. We simulate data for illustration.

``` r
library(MetaAnalysis)
set.seed(123)  ## Fix seed for reproducible result
p <- 5
# dataframe1 has 5 columns and 3 groups
data1 <- data.frame(group=sample(1:3,200,replace=TRUE),matrix(rnorm(p*200),ncol=p))

head(data1)
#>   group          X1         X2          X3         X4         X5
#> 1     1 -0.71040656 -0.7152422 -0.60189285 -0.7282191 -1.0141142
#> 2     3  0.25688371 -0.7526890 -0.99369859 -1.5404424 -0.7913139
#> 3     2 -0.24669188 -0.9385387  1.02678506 -0.6930946  0.2995937
#> 4     3 -0.34754260 -1.0525133  0.75106130  0.1188494  1.6390519
#> 5     3 -0.95161857 -0.4371595 -1.50916654 -1.3647095  1.0846170
#> 6     1 -0.04502772  0.3311792 -0.09514745  0.5899827 -0.6245675
```

data2 and data3 are generated in the same fashion.

``` r
# dataframe1 has 5 columns and 2 groups
data2 <- data.frame(group=sample(1:2,150,replace=TRUE),matrix(rnorm(p*150),ncol=p))
# dataframe1 has 5 columns and 4 groups
data3 <- data.frame(group=sample(1:4,400,replace=TRUE),matrix(rnorm(p*400),ncol=p))
```

To do the meta analysis using this package, input dataframes sequentially in the first three parameters, and specify pooling method/methods to use if wanted (as introduced in previous paragraph). We can also customize significance level.

``` r
res = meta.analysis(data1, data2, data3, method=c('Fisher', 'Stouffer', 'minP', 'maxP'), alpha = 0.05)
```

*Input data error: if number of dataframes are out of limit, or number of columns are not equal.*

*Any row with missing values will be dropped in the package*

1.  Obtaining p-value matrix
    -   To view the p-values matrix, one can access the `p.matrix` component of resluting variable.

``` r
res$p.matrix
#>              X1        X2        X3        X4        X5
#> Data1 0.8391408 0.3286983 0.9805006 0.3796453 0.6302030
#> Data2 0.2160643 0.8853164 0.8557419 0.2849488 0.3575808
#> Data3 0.6228494 0.6347043 0.5886171 0.5202245 0.5911930
```

1.  Tests performed
    -   User can also get a record of the results from tests performed, which is corresponding to p-value matrix.

``` r
res$test.performed
#>            X1    X2      X3      X4    X5
#> Data1 Kruskal Anova   Anova Kruskal Anova
#> Data2  Wilcox T(eq)   T(eq)   T(eq) T(eq)
#> Data3   Anova Anova Kruskal Kruskal Anova
```

1.  Pooling p-values
    -   Finally, show the pooled p-value by accessing `pooled.p.matrix` element.

``` r
res$pooled.p.matrix
#>                 X1        X2        X3        X4        X5
#> Fisher   0.6278110 0.7601185 0.9651910 0.4511969 0.6724216
#> Stouffer 0.6176331 0.7378441 0.9734358 0.3171520 0.5455228
#> minP     0.5182282 0.6974805 0.9303792 0.6343955 0.7348721
#> maxP     0.5908870 0.6938979 0.9426351 0.1407902 0.2502888
```

Example (Comparison)
--------------------

Suppose we want to see the visualized "trend" between each pooling method, we can plot the matrix using `matplot` function

``` r
set.seed(123)
p <- 20
data1 <- data.frame(group=sample(1:3,20,replace=TRUE), matrix(rnorm(p*20),ncol=p))
data2 <- data.frame(group=sample(1:2,15,replace=TRUE), matrix(rnorm(p*15),ncol=p))
data3 <- data.frame(group=sample(1:4,40,replace=TRUE), matrix(rnorm(p*40),ncol=p))
data4 <- data.frame(group=sample(1:3,50,replace=TRUE), matrix(rnorm(p*50),ncol=p))
data5 <- data.frame(group=sample(1:2,40,replace=TRUE), matrix(rnorm(p*40),ncol=p))

res = meta.analysis(data1, data2, data3, data4, data5,
                    method=c('Fisher',  'maxP', 'Stouffer', 'minP'), 
                    alpha = 0.05)
pooled.p.matrix = res$pooled.p.matrix

matplot(t(pooled.p.matrix), type='l', col = c('black', 'green', 'red', 'blue'), lty = 'solid')
legend(0, 0, inset=c(-10,-10), legend = rownames(pooled.p.matrix), xpd = T, col = c('black', 'green', 'red', 'blue'), lty = 'solid')
```

![](Introduction_files/figure-markdown_github/unnamed-chunk-8-1.png)
