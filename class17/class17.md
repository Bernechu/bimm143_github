# class 17
Berne Chu (A18608434)

## Section 1. Proportion og G/G in a population

Read CSV file downloaded from emsembl

``` r
mxl <- read.csv("class17data.csv")
head(mxl)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  NA19648 (F)                       A|A ALL, AMR, MXL      -
    2                  NA19649 (M)                       G|G ALL, AMR, MXL      -
    3                  NA19651 (F)                       A|A ALL, AMR, MXL      -
    4                  NA19652 (M)                       G|G ALL, AMR, MXL      -
    5                  NA19654 (F)                       G|G ALL, AMR, MXL      -
    6                  NA19655 (M)                       A|G ALL, AMR, MXL      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
mxl$Genotype <- as.factor(mxl$Genotype)
table(mxl$Genotype)
```


    A|A A|G G|A G|G 
     22  21  12   9 

``` r
table(mxl$Genotype) / nrow(mxl) * 100
```


        A|A     A|G     G|A     G|G 
    34.3750 32.8125 18.7500 14.0625 
