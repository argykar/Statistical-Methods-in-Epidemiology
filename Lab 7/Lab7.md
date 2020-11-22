Statistical Methods in Epidemiology - Lab 7
================

## Table of Contents

1.  [Matched studies](#matched)
    -   [Endometrial cancer, 1:1 matched sets](#1:1)
    -   [Endometrial cancer, 1:4 matched sets](#1:4)
2.  [Conditional logistic regression](#conditional)
3.  [Conditional vs. unconditional logistic](#vs)

# Matched case-control studies

In this lab, we are going to display some ways to analyse data from
matched case control studies in R. As usual, we’ll begin by importing
the libraries that are going to be used.

``` r
library(haven)
library(kableExtra)
library(survival)
library(tidyverse)
library(lmtest)
library(Epi)
```

## I. Matched studies <a name='matched'></a>

We are going to use two datasets, *bdendo* and *bdendo11*. The first one
comes from the study of endometrial cancer and the medical use of
estrogens, described in Breslow and Day (1980). Each case was matched to
4 control women who were alive and living in the retirement community at
the time the case was diagnosed, who were born within one year of the
case, had the same marital status and who had entered the community at
approximately the same time. The variable **age3** codes subjects in
three board groups (55-64, 65-74 and 75+) using the age of the case in
each set, so subjects in the same matched set always have the same value
for **age3**. The variable **agegrp** does the same thing for 5 year age
groups. Variable **d** indicates the case-control status (1:case,
0:control), while variable **est** denotes the use of estrogens. The
second dataset is a subset of the first one, and in particular for each
case one control woman has been chosen (1:1 format). We’ll load both of
them and create factors for the categorical variables.

``` r
data1 <- read_dta("bdendo11.dta")
# factoring
data1$gall <- factor(data1$gall, labels = c("No", "Yes"))
data1$hyp <- factor(data1$hyp, labels = c("No", "Yes"))
data1$ob <- factor(data1$ob, labels = c("No", "Yes"))
data1$est <- factor(data1$est, labels = c("No", "Yes"))
data1$dur <- factor(data1$dur, ordered = TRUE)
data1$non <- factor(data1$non, labels = c("No", "Yes"))
data1$cest <- factor(data1$cest, ordered = TRUE)
data1$agegrp <- factor(data1$agegrp, labels = c("55-59", "60-64", "65-69", "70-74", "75-79", "80-84"))
data1$age3 <- factor(data1$age3, labels = c("<64", "65-75", "75+"))


data2 <- read_dta("bdendo.dta")
# factoring
data2$gall <- factor(data2$gall, labels = c("No", "Yes"))
data2$hyp <- factor(data2$hyp, labels = c("No", "Yes"))
data2$ob <- factor(data2$ob, labels = c("No", "Yes"))
data2$est <- factor(data2$est, labels = c("No", "Yes"))
data2$dur <- factor(data2$dur, ordered = TRUE)
data2$non <- factor(data2$non, labels = c("No", "Yes"))
data2$cest <- factor(data2$cest, ordered = TRUE)
data2$agegrp <- factor(data2$agegrp, labels = c("55-59", "60-64", "65-69", "70-74", "75-79", "80-84"))
data2$age3 <- factor(data2$age3, labels = c("<64", "65-75", "75+"))


# simple function for not writing all the time kable functions
# default for digits is 3, but you can modify if you want more or less
display <- function(x, digits = 3) {
  x %>%
    kable(digits = digits) %>%
    kable_styling()
}
```

### i) Endometrial cancer, 1:1 matched sets <a name='1:1'></a>

In this section we’ll work with **bdendo11** dataset. First, we’ll
examine the effect of estrogen with Mantel - Haenszel. We’ll work with
`mantelhaen.test()` function. In the previous lab, we mainly used
`mhor()`, as it displays the odds ratios of every stratum. Here, where
strata are defined by the sets, this information has no use and
`mantelhaen.test()` does not provide such information, so it works as we
want. In the code below, we are extracting the pieces we want from an
object created by `mantelhaen.test()` in order to present them in a
table.

``` r
test1 <- with(data1, mantelhaen.test(table(d, est, set), correct = FALSE))

bind_cols("OR" = test1$estimate, 
           "chi2(1)" = test1$statistic, 
           "p>chi2" = test1$p.value,
           "Lower_l" = test1$conf.int[1],
          "Upper_l" = test1$conf.int[2]) %>%
  display(digits = 4)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
chi2(1)
</th>
<th style="text-align:right;">
p&gt;chi2
</th>
<th style="text-align:right;">
Lower\_l
</th>
<th style="text-align:right;">
Upper\_l
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
9.6667
</td>
<td style="text-align:right;">
21.125
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.9447
</td>
<td style="text-align:right;">
31.7331
</td>
</tr>
</tbody>
</table>

The Mantel - Haenszel estimate is the ratio of the number of matched
pairs where only the case is exposed to the number of matched pairs
where only the control is exposed. We can verify it by the code below.
Here, we are grouping by set (and **age3**, but since **age3** remains
constant within set it is like grouping only by set, we just want to
have the age available in the dataframe below) and we are creating two
variables, one that counts the total number of exposed people within
sets (can be 0, 1, or 2) and one that counts the total number of exposed
cases within sets (can be 0 or 1). By keeping only the sets that contain
exactly one exposed person, we are essentially keeping the sets where
either have exactly one exposed case or exactly one exposed control.
Thus, in the final dataset, **Est\_ca** is equal to one in sets where
only the case is exposed and equal to zero in sets where only the
control is exposed.

``` r
m_pairs <- data1 %>%
  group_by(set, age3) %>%
  summarise(Est = sum((est == "Yes")), 
            Est_Ca = sum((est == "Yes") & (d == 1))) %>% 
  filter(Est == 1)

sum(m_pairs$Est_Ca == 1) / sum(m_pairs$Est_Ca == 0)
```

    ## [1] 9.666667

There is no use controlling for **age3** as people are matched on age,
but we can check whether age modifies the effect of estrogen. So, we’ll
calculate Mantel - Haenszel estimate of odds ratio controlled for set in
every stratum of **age3**. We’ll accomplish it by applying
`mantelhaen.test()` to every level of **age3** through `group_map()`.
Note the use of `sapply()` that helps us extract the same information as
before for every level of **age3**.

``` r
test2 <- data1 %>%
  group_by(age3) %>%
  group_map(~ with(.x, mantelhaen.test(table(d, est, set), correct = FALSE)))

bind_cols("Age" = levels(data1$age3),
           "OR" = sapply(1:3, function(i) test2[[i]]$estimate),
           "chi2(1)" = sapply(1:3, function(i) test2[[i]]$statistic),
           "p>chi2" = sapply(1:3, function(i) test2[[i]]$p.value),
           "Lower_l" = sapply(1:3, function(i) test2[[i]]$conf.int[1]),
           "Upper_l" = sapply(1:3, function(i) test2[[i]]$conf.int[2])) %>%
  display(digits = 4)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Age
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
chi2(1)
</th>
<th style="text-align:right;">
p&gt;chi2
</th>
<th style="text-align:right;">
Lower\_l
</th>
<th style="text-align:right;">
Upper\_l
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
&lt;64
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
3.5714
</td>
<td style="text-align:right;">
0.0588
</td>
<td style="text-align:right;">
0.7224
</td>
<td style="text-align:right;">
49.8372
</td>
</tr>
<tr>
<td style="text-align:left;">
65-75
</td>
<td style="text-align:right;">
15
</td>
<td style="text-align:right;">
12.2500
</td>
<td style="text-align:right;">
0.0005
</td>
<td style="text-align:right;">
1.9814
</td>
<td style="text-align:right;">
113.5556
</td>
</tr>
<tr>
<td style="text-align:left;">
75+
</td>
<td style="text-align:right;">
8
</td>
<td style="text-align:right;">
5.4444
</td>
<td style="text-align:right;">
0.0196
</td>
<td style="text-align:right;">
1.0006
</td>
<td style="text-align:right;">
63.9625
</td>
</tr>
</tbody>
</table>

The chi-square test of heterogeneity of odds ratios between the strata
of **age3** is in fact, in this case, a Pearson Chi-square performed
only in the sets which have exactly one exposed case or exactly one
exposed control, by tabulating the age and the variable that shows
whether the case or the control is exposed. So, we are going to use the
**m\_pairs** dataframe that we created before.

``` r
with(m_pairs, chisq.test(Est_Ca, age3))
```

    ## 
    ##  Pearson's Chi-squared test
    ## 
    ## data:  Est_Ca and age3
    ## X-squared = 0.41452, df = 2, p-value = 0.8128

Afterwards, we’ll try to control for hypertension. As
`mantelhaen.test()` does not support stratification by more than one
variable, we have to manually create a variable from the levels of
**set** and **hyp**. We will use `interaction()`, which automatically
does that work, and specifying the argument *drop = TRUE*, unused
combinations of levels are dropped. Also, it’s necessary to drop sets in
which case and control have different hypertension, as it will result in
strata of set and hypertension with a single observation and
`mantelhaen.test()` supports only strata with more than one observation.
The forementioned tasks will be executed with the aid of **Tidyverse**
tools.

``` r
help_data <- data1 %>% 
  group_by(set, hyp) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  mutate(set_hyp = interaction(factor(set), hyp, drop = TRUE))
  
test3 <- with(help_data, mantelhaen.test(d, est, set_hyp, correct = FALSE))

bind_cols("OR" = test3$estimate, 
           "chi2(1)" = test3$statistic, 
           "p>chi2" = test3$p.value,
           "Lower_l" = test3$conf.int[1],
          "Upper_l" = test3$conf.int[2]) %>% 
  display(digits = 5)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
chi2(1)
</th>
<th style="text-align:right;">
p&gt;chi2
</th>
<th style="text-align:right;">
Lower\_l
</th>
<th style="text-align:right;">
Upper\_l
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
17
</td>
<td style="text-align:right;">
14.22222
</td>
<td style="text-align:right;">
0.00016
</td>
<td style="text-align:right;">
2.2624
</td>
<td style="text-align:right;">
127.7403
</td>
</tr>
</tbody>
</table>

### ii) Endometrial cancer, 1:4 matched sets <a name='1:4'></a>

Here, we’ll work with the full dataset (i.e 1:4 matched sets) and we’ll
answer the same questions as before. With the code below, we are
obtaining the Mantel - Haenszel estimate of the Odds Ratio for the use
of estrogens.

``` r
test4 <- with(data2, mantelhaen.test(table(d, est, set), correct = FALSE))

bind_cols("OR" = test4$estimate, 
           "chi2(1)" = test4$statistic, 
           "p>chi2" = test4$p.value,
           "Lower_l" = test4$conf.int[1],
          "Upper_l" = test4$conf.int[2]) %>% 
  display(digits = 5)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
chi2(1)
</th>
<th style="text-align:right;">
p&gt;chi2
</th>
<th style="text-align:right;">
Lower\_l
</th>
<th style="text-align:right;">
Upper\_l
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
8.46154
</td>
<td style="text-align:right;">
31.15563
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
3.41153
</td>
<td style="text-align:right;">
20.98696
</td>
</tr>
</tbody>
</table>

Next, we are calculating the Mantel - Haenszel estimate of the Odds
Ratio for the use of estrogens in each category of **age3**.

``` r
test5 <- data2 %>%
  group_by(age3) %>%
  group_map(~ with(.x, mantelhaen.test(table(d, est, set), correct = FALSE)))

bind_cols("Age" = levels(data2$age3),
           "OR" = sapply(1:3, function(i) test5[[i]]$estimate),
           "chi2(1)" = sapply(1:3, function(i) test5[[i]]$statistic),
           "p>chi2" = sapply(1:3, function(i) test5[[i]]$p.value),
           "Lower_l" = sapply(1:3, function(i) test5[[i]]$conf.int[1]),
           "Upper_l" = sapply(1:3, function(i) test5[[i]]$conf.int[2])) %>%
  display(digits = 4)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Age
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
chi2(1)
</th>
<th style="text-align:right;">
p&gt;chi2
</th>
<th style="text-align:right;">
Lower\_l
</th>
<th style="text-align:right;">
Upper\_l
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
&lt;64
</td>
<td style="text-align:right;">
3.8000
</td>
<td style="text-align:right;">
3.3793
</td>
<td style="text-align:right;">
0.0660
</td>
<td style="text-align:right;">
0.7307
</td>
<td style="text-align:right;">
19.7606
</td>
</tr>
<tr>
<td style="text-align:left;">
65-75
</td>
<td style="text-align:right;">
10.6667
</td>
<td style="text-align:right;">
18.6889
</td>
<td style="text-align:right;">
0.0000
</td>
<td style="text-align:right;">
2.8733
</td>
<td style="text-align:right;">
39.5981
</td>
</tr>
<tr>
<td style="text-align:left;">
75+
</td>
<td style="text-align:right;">
13.5000
</td>
<td style="text-align:right;">
9.7656
</td>
<td style="text-align:right;">
0.0018
</td>
<td style="text-align:right;">
2.1099
</td>
<td style="text-align:right;">
86.3788
</td>
</tr>
</tbody>
</table>

Finally, we are going to obtain the Mantel - Haenszel estimate of the
Odds Ratio for the use of estrogens adjusted for hypertension.

``` r
help_data_II <- data2 %>% 
  group_by(set, hyp) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  mutate(set_hyp = interaction(factor(set), hyp, drop = TRUE))
  
test6 <- with(help_data_II, mantelhaen.test(d, est, set_hyp, correct = FALSE))


bind_cols("OR" = test6$estimate, 
           "chi2(1)" = test6$statistic, 
           "p>chi2" = test6$p.value,
           "Lower_l" = test6$conf.int[1],
          "Upper_l" = test6$conf.int[2]) %>% 
  display(digits = 5)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
chi2(1)
</th>
<th style="text-align:right;">
p&gt;chi2
</th>
<th style="text-align:right;">
Lower\_l
</th>
<th style="text-align:right;">
Upper\_l
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
29.22222
</td>
<td style="text-align:right;">
33.47134
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
7.24748
</td>
<td style="text-align:right;">
117.8255
</td>
</tr>
</tbody>
</table>

## II. Conditional logistic regression <a name='conditional'></a>

In this section, we are going to use conditional logistic regression for
the analysis. One function which we can use for that purpose is
`clogit()` from **survival** package, which has similar syntax with the
rest of the functions we have seen that fit models. The *formula*
argument has to be in the form of **y \~ x + strata(matched.set)**. The
*matched.set* variable has to be numeric. First, we’ll work with the 1:1
matched dataset. We are going to examine the effect of estrogen.

For not repeating the same code chunk, below we created a function which
make the results appear in a nice table format like STATA. We added an
option for exponentiating the coefficients if you’d like. Default option
is **exp = TRUE**. We use `na.omit()`, so as not to include lines that
contain *NAs*, because of the matching variables coefficients.

``` r
table_stata <- function(model, digits = 3, exp = TRUE) {
  df <- as.data.frame(summary(model)$coefficients) %>%
    rownames_to_column(var = "Predictors") %>%
    bind_cols(as.data.frame(summary(model)$conf.int) %>% select(-`exp(coef)`)) %>%
    mutate(SE = `exp(coef)` * `se(coef)`,
           lower = `coef` - 1.96 * `se(coef)`,
           upper = `coef` + 1.96 * `se(coef)`) %>%
    na.omit()
  
    {if (exp)
      df %>% select(Predictors,
             'OR' = `exp(coef)`,
             SE,
             `z`,
             `Pr(>|z|)`,
             `lower .95`,
             `upper .95`)
    else
      df %>% select(
        Predictors,
        'Coef' = `coef`,
        'se' = `se(coef)`,
        `z`,
        `Pr(>|z|)`,
        lower,
        upper
      )
  } %>%
    display(digits = digits)
}
```

``` r
fit1 <- clogit(d ~ est + strata(set), data1)

table_stata(fit1)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
9.667
</td>
<td style="text-align:right;">
5.863
</td>
<td style="text-align:right;">
3.741
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
2.945
</td>
<td style="text-align:right;">
31.733
</td>
</tr>
</tbody>
</table>

Next, we are going to examine whether the effect of age is modified by
**age3**. To do this, first we have to fit a model with an interaction
of **est** and **age3**. However, age is the matching variable, so there
will be no calculations for the main effects of **age3**. Subsequently,
we are performing a likelihood ratio test to test the interaction
parameters.

``` r
fit2 <- clogit(d ~ est * age3 + strata(set), data1)

table_stata(fit2)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
6.000
</td>
<td style="text-align:right;">
6.481
</td>
<td style="text-align:right;">
1.659
</td>
<td style="text-align:right;">
0.097
</td>
<td style="text-align:right;">
0.722
</td>
<td style="text-align:right;">
49.837
</td>
</tr>
<tr>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
estYes:age365-75
</td>
<td style="text-align:right;">
2.500
</td>
<td style="text-align:right;">
3.736
</td>
<td style="text-align:right;">
0.613
</td>
<td style="text-align:right;">
0.540
</td>
<td style="text-align:right;">
0.134
</td>
<td style="text-align:right;">
46.774
</td>
</tr>
<tr>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
estYes:age375+
</td>
<td style="text-align:right;">
1.333
</td>
<td style="text-align:right;">
2.018
</td>
<td style="text-align:right;">
0.190
</td>
<td style="text-align:right;">
0.849
</td>
<td style="text-align:right;">
0.069
</td>
<td style="text-align:right;">
25.912
</td>
</tr>
</tbody>
</table>

``` r
lrtest(fit1, fit2) %>%
  display()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
\#Df
</th>
<th style="text-align:right;">
LogLik
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
Chisq
</th>
<th style="text-align:right;">
Pr(&gt;Chisq)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-31.444
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-31.239
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.815
</td>
</tr>
</tbody>
</table>

Same for variable **age**, which is continuous (age is the matching
variable, but in its continuous form there are a few differences in its
values among people of one set, so we have main-effects estimation).

``` r
fit3 <- clogit(d ~ est + age + strata(set), data1)

table_stata(fit3)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
10.249
</td>
<td style="text-align:right;">
6.392
</td>
<td style="text-align:right;">
3.732
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
3.019
</td>
<td style="text-align:right;">
34.796
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0.741
</td>
<td style="text-align:right;">
0.299
</td>
<td style="text-align:right;">
-0.744
</td>
<td style="text-align:right;">
0.457
</td>
<td style="text-align:right;">
0.336
</td>
<td style="text-align:right;">
1.633
</td>
</tr>
</tbody>
</table>

``` r
fit4 <- clogit(d ~ est*age + strata(set), data1)

table_stata(fit4)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
2.348
</td>
<td style="text-align:right;">
15.058
</td>
<td style="text-align:right;">
0.133
</td>
<td style="text-align:right;">
0.894
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
675908.171
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0.712
</td>
<td style="text-align:right;">
0.313
</td>
<td style="text-align:right;">
-0.771
</td>
<td style="text-align:right;">
0.441
</td>
<td style="text-align:right;">
0.301
</td>
<td style="text-align:right;">
1.688
</td>
</tr>
<tr>
<td style="text-align:left;">
estYes:age
</td>
<td style="text-align:right;">
1.021
</td>
<td style="text-align:right;">
0.092
</td>
<td style="text-align:right;">
0.230
</td>
<td style="text-align:right;">
0.818
</td>
<td style="text-align:right;">
0.855
</td>
<td style="text-align:right;">
1.218
</td>
</tr>
</tbody>
</table>

``` r
lrtest(fit3, fit4) %>%
  display()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
\#Df
</th>
<th style="text-align:right;">
LogLik
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
Chisq
</th>
<th style="text-align:right;">
Pr(&gt;Chisq)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-31.170
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-31.143
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.053
</td>
<td style="text-align:right;">
0.818
</td>
</tr>
</tbody>
</table>

We are now going to adjust for hypertension. Below we can see that the
effect is increased, however not a statistically significant change on
hypertension.

``` r
fit5 <- clogit(d ~ est + hyp + strata(set), data1)

table_stata(fit5)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
10.049
</td>
<td style="text-align:right;">
6.231
</td>
<td style="text-align:right;">
3.721
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
2.981
</td>
<td style="text-align:right;">
33.878
</td>
</tr>
<tr>
<td style="text-align:left;">
hypYes
</td>
<td style="text-align:right;">
0.867
</td>
<td style="text-align:right;">
0.379
</td>
<td style="text-align:right;">
-0.326
</td>
<td style="text-align:right;">
0.744
</td>
<td style="text-align:right;">
0.368
</td>
<td style="text-align:right;">
2.043
</td>
</tr>
</tbody>
</table>

Check for interaction between hypertension and estrogen. We are going to
test for the significance of the interaction using LR test with the
models **fit5** and **fit6**

``` r
fit6 <- clogit(d ~ est * hyp + strata(set), data1)

table_stata(fit6)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
11.539
</td>
<td style="text-align:right;">
8.259
</td>
<td style="text-align:right;">
3.417
</td>
<td style="text-align:right;">
0.001
</td>
<td style="text-align:right;">
2.837
</td>
<td style="text-align:right;">
46.927
</td>
</tr>
<tr>
<td style="text-align:left;">
hypYes
</td>
<td style="text-align:right;">
1.250
</td>
<td style="text-align:right;">
1.219
</td>
<td style="text-align:right;">
0.229
</td>
<td style="text-align:right;">
0.819
</td>
<td style="text-align:right;">
0.185
</td>
<td style="text-align:right;">
8.451
</td>
</tr>
<tr>
<td style="text-align:left;">
estYes:hypYes
</td>
<td style="text-align:right;">
0.641
</td>
<td style="text-align:right;">
0.684
</td>
<td style="text-align:right;">
-0.417
</td>
<td style="text-align:right;">
0.677
</td>
<td style="text-align:right;">
0.079
</td>
<td style="text-align:right;">
5.197
</td>
</tr>
</tbody>
</table>

``` r
lrtest(fit5, fit6) %>%
  display()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
\#Df
</th>
<th style="text-align:right;">
LogLik
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
Chisq
</th>
<th style="text-align:right;">
Pr(&gt;Chisq)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-31.390
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-31.305
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.171
</td>
<td style="text-align:right;">
0.679
</td>
</tr>
</tbody>
</table>

Results from the LR test show that the interaction is not significant
and hence dropped of the model.

From now on, we will use the **bdendo.dta** file. We will use matching
1:4 but this time with models and not with mantel-haenzel tests. At the
beginning, we factored every variable needed to be factored and now we
are ready to proceed with the analysis.

First let’s take a look at our data.

``` r
data2 %>%
  select(set, d, est, age) %>%
  arrange(set) %>%
  head(10) %>%
  display()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
set
</th>
<th style="text-align:right;">
d
</th>
<th style="text-align:left;">
est
</th>
<th style="text-align:right;">
age
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:right;">
75
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:right;">
74
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:right;">
75
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:right;">
74
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:right;">
74
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:right;">
67
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:right;">
67
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:right;">
67
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
Yes
</td>
<td style="text-align:right;">
68
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:left;">
No
</td>
<td style="text-align:right;">
67
</td>
</tr>
</tbody>
</table>

As before, we create a fit object conducting a conditional logistic with
variables estrogen and groups the set variable.

``` r
fit11 <- clogit(d ~ est + strata(set), data2)

table_stata(fit11, exp = FALSE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
Coef
</th>
<th style="text-align:right;">
se
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower
</th>
<th style="text-align:right;">
upper
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
2.074
</td>
<td style="text-align:right;">
0.421
</td>
<td style="text-align:right;">
4.928
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.249
</td>
<td style="text-align:right;">
2.899
</td>
</tr>
</tbody>
</table>

``` r
# not meaningful since set is the matching variable and we can't get the estimates
fit12 <- clogit(d ~ age3 + strata(set), data2)
```

    ## Warning in fitter(X, Y, istrat, offset, init, control, weights = weights, :
    ## Loglik converged before variable 1,2 ; beta may be infinite.

We are now creating models with interaction effects and without between
the variables specified in the models as shown below. Then we are
testing the interaction significance with LR tests.

``` r
fit13 <- clogit(d ~ est * age3 + strata(set), data2)

table_stata(fit13, exp = FALSE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
Coef
</th>
<th style="text-align:right;">
se
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower
</th>
<th style="text-align:right;">
upper
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1
</td>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
1.431
</td>
<td style="text-align:right;">
0.826
</td>
<td style="text-align:right;">
1.733
</td>
<td style="text-align:right;">
0.083
</td>
<td style="text-align:right;">
-0.188
</td>
<td style="text-align:right;">
3.049
</td>
</tr>
<tr>
<td style="text-align:left;">
4
</td>
<td style="text-align:left;">
estYes:age365-75
</td>
<td style="text-align:right;">
0.847
</td>
<td style="text-align:right;">
1.034
</td>
<td style="text-align:right;">
0.820
</td>
<td style="text-align:right;">
0.412
</td>
<td style="text-align:right;">
-1.179
</td>
<td style="text-align:right;">
2.874
</td>
</tr>
<tr>
<td style="text-align:left;">
5
</td>
<td style="text-align:left;">
estYes:age375+
</td>
<td style="text-align:right;">
0.780
</td>
<td style="text-align:right;">
1.154
</td>
<td style="text-align:right;">
0.676
</td>
<td style="text-align:right;">
0.499
</td>
<td style="text-align:right;">
-1.482
</td>
<td style="text-align:right;">
3.042
</td>
</tr>
</tbody>
</table>

``` r
fit14 <- clogit(d ~ est + age3 + strata(set), data2)

table_stata(fit14, exp = FALSE)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
Coef
</th>
<th style="text-align:right;">
se
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower
</th>
<th style="text-align:right;">
upper
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
2.074
</td>
<td style="text-align:right;">
0.421
</td>
<td style="text-align:right;">
4.928
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.249
</td>
<td style="text-align:right;">
2.899
</td>
</tr>
</tbody>
</table>

``` r
lrtest(fit14, fit13) %>%
  display()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
\#Df
</th>
<th style="text-align:right;">
LogLik
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
Chisq
</th>
<th style="text-align:right;">
Pr(&gt;Chisq)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
-83.722
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-83.380
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.683
</td>
<td style="text-align:right;">
0.711
</td>
</tr>
</tbody>
</table>

``` r
fit15 <- clogit(d ~ est * age + strata(set), data2)

table_stata(fit15)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
2.077
</td>
<td style="text-align:right;">
9.452
</td>
<td style="text-align:right;">
0.161
</td>
<td style="text-align:right;">
0.872
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
15548.544
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0.734
</td>
<td style="text-align:right;">
0.185
</td>
<td style="text-align:right;">
-1.224
</td>
<td style="text-align:right;">
0.221
</td>
<td style="text-align:right;">
0.448
</td>
<td style="text-align:right;">
1.204
</td>
</tr>
<tr>
<td style="text-align:left;">
estYes:age
</td>
<td style="text-align:right;">
1.019
</td>
<td style="text-align:right;">
0.065
</td>
<td style="text-align:right;">
0.300
</td>
<td style="text-align:right;">
0.764
</td>
<td style="text-align:right;">
0.900
</td>
<td style="text-align:right;">
1.154
</td>
</tr>
</tbody>
</table>

``` r
fit16 <- clogit(d ~ est + age + strata(set), data2)

table_stata(fit16)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
8.128
</td>
<td style="text-align:right;">
3.462
</td>
<td style="text-align:right;">
4.919
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
3.527
</td>
<td style="text-align:right;">
18.730
</td>
</tr>
<tr>
<td style="text-align:left;">
age
</td>
<td style="text-align:right;">
0.750
</td>
<td style="text-align:right;">
0.181
</td>
<td style="text-align:right;">
-1.189
</td>
<td style="text-align:right;">
0.234
</td>
<td style="text-align:right;">
0.467
</td>
<td style="text-align:right;">
1.205
</td>
</tr>
</tbody>
</table>

``` r
lrtest(fit15, fit16) %>%
  display()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
\#Df
</th>
<th style="text-align:right;">
LogLik
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
Chisq
</th>
<th style="text-align:right;">
Pr(&gt;Chisq)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-82.963
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-83.008
</td>
<td style="text-align:right;">
-1
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.764
</td>
</tr>
</tbody>
</table>

``` r
fit17 <- clogit(d ~ est + hyp + strata(set), data2)

table_stata(fit17)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
7.961
</td>
<td style="text-align:right;">
3.399
</td>
<td style="text-align:right;">
4.858
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
3.447
</td>
<td style="text-align:right;">
18.384
</td>
</tr>
<tr>
<td style="text-align:left;">
hypYes
</td>
<td style="text-align:right;">
0.996
</td>
<td style="text-align:right;">
0.329
</td>
<td style="text-align:right;">
-0.011
</td>
<td style="text-align:right;">
0.991
</td>
<td style="text-align:right;">
0.521
</td>
<td style="text-align:right;">
1.905
</td>
</tr>
</tbody>
</table>

``` r
fit18 <- clogit(d ~ est * hyp + strata(set), data2)

table_stata(fit18)
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
SE
</th>
<th style="text-align:right;">
z
</th>
<th style="text-align:right;">
Pr(&gt;\|z\|)
</th>
<th style="text-align:right;">
lower .95
</th>
<th style="text-align:right;">
upper .95
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
estYes
</td>
<td style="text-align:right;">
9.496
</td>
<td style="text-align:right;">
4.795
</td>
<td style="text-align:right;">
4.458
</td>
<td style="text-align:right;">
0.000
</td>
<td style="text-align:right;">
3.530
</td>
<td style="text-align:right;">
25.546
</td>
</tr>
<tr>
<td style="text-align:left;">
hypYes
</td>
<td style="text-align:right;">
1.858
</td>
<td style="text-align:right;">
1.632
</td>
<td style="text-align:right;">
0.705
</td>
<td style="text-align:right;">
0.481
</td>
<td style="text-align:right;">
0.332
</td>
<td style="text-align:right;">
10.393
</td>
</tr>
<tr>
<td style="text-align:left;">
estYes:hypYes
</td>
<td style="text-align:right;">
0.502
</td>
<td style="text-align:right;">
0.459
</td>
<td style="text-align:right;">
-0.754
</td>
<td style="text-align:right;">
0.451
</td>
<td style="text-align:right;">
0.084
</td>
<td style="text-align:right;">
3.007
</td>
</tr>
</tbody>
</table>

``` r
lrtest(fit17, fit18) %>%
  display()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
\#Df
</th>
<th style="text-align:right;">
LogLik
</th>
<th style="text-align:right;">
Df
</th>
<th style="text-align:right;">
Chisq
</th>
<th style="text-align:right;">
Pr(&gt;Chisq)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
-83.722
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
-83.458
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.526
</td>
<td style="text-align:right;">
0.468
</td>
</tr>
</tbody>
</table>

## III. Conditional vs. unconditional logistic <a name='vs'></a>

``` r
data3 <- read_dta('lep.dta')
data3$d <- factor(data3$d, levels = c("1", "0"), labels = c("Cases", "Controls"))
data3$sex <- factor(data3$sex, labels = c("Male", "Female"))
data3$bcg <- factor(data3$bcg, levels = c("1", "0"), labels = c("Present", "Absent"))

data3$cage <- factor(data3$age, labels = c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34"))

data3$d <- fct_relevel(data3$d, "Controls")
data3$bcg <- fct_relevel(data3$bcg, "Absent")

fit31 <- glm(d ~ bcg + cage, data3, family = binomial)
# function from previous labs
r_table(fit31) %>% 
  display()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Predictors
</th>
<th style="text-align:right;">
OR
</th>
<th style="text-align:right;">
Std\_Err
</th>
<th style="text-align:right;">
z\_value
</th>
<th style="text-align:left;">
P-Value
</th>
<th style="text-align:right;">
LL
</th>
<th style="text-align:right;">
UL
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
bcgPresent
</td>
<td style="text-align:right;">
0.331
</td>
<td style="text-align:right;">
0.064
</td>
<td style="text-align:right;">
-5.732
</td>
<td style="text-align:left;">
&lt;0.001
</td>
<td style="text-align:right;">
0.227
</td>
<td style="text-align:right;">
0.483
</td>
</tr>
<tr>
<td style="text-align:left;">
cage5-9
</td>
<td style="text-align:right;">
1.339
</td>
<td style="text-align:right;">
0.351
</td>
<td style="text-align:right;">
1.113
</td>
<td style="text-align:left;">
0.266
</td>
<td style="text-align:right;">
0.801
</td>
<td style="text-align:right;">
2.240
</td>
</tr>
<tr>
<td style="text-align:left;">
cage10-14
</td>
<td style="text-align:right;">
1.237
</td>
<td style="text-align:right;">
0.374
</td>
<td style="text-align:right;">
0.704
</td>
<td style="text-align:left;">
0.481
</td>
<td style="text-align:right;">
0.684
</td>
<td style="text-align:right;">
2.238
</td>
</tr>
<tr>
<td style="text-align:left;">
cage15-19
</td>
<td style="text-align:right;">
1.625
</td>
<td style="text-align:right;">
0.523
</td>
<td style="text-align:right;">
1.508
</td>
<td style="text-align:left;">
0.131
</td>
<td style="text-align:right;">
0.865
</td>
<td style="text-align:right;">
3.052
</td>
</tr>
<tr>
<td style="text-align:left;">
cage20-24
</td>
<td style="text-align:right;">
3.402
</td>
<td style="text-align:right;">
1.003
</td>
<td style="text-align:right;">
4.153
</td>
<td style="text-align:left;">
&lt;0.001
</td>
<td style="text-align:right;">
1.909
</td>
<td style="text-align:right;">
6.063
</td>
</tr>
<tr>
<td style="text-align:left;">
cage25-29
</td>
<td style="text-align:right;">
2.422
</td>
<td style="text-align:right;">
0.596
</td>
<td style="text-align:right;">
3.596
</td>
<td style="text-align:left;">
&lt;0.001
</td>
<td style="text-align:right;">
1.496
</td>
<td style="text-align:right;">
3.923
</td>
</tr>
<tr>
<td style="text-align:left;">
cage30-34
</td>
<td style="text-align:right;">
2.085
</td>
<td style="text-align:right;">
0.500
</td>
<td style="text-align:right;">
3.064
</td>
<td style="text-align:left;">
0.002
</td>
<td style="text-align:right;">
1.303
</td>
<td style="text-align:right;">
3.335
</td>
</tr>
</tbody>
</table>

Now, we create a variable that is 1 for every subject and we are going
to declare it the matching variable. Next, we proceed with conditional
logistic regression similar with the model **fit31** but with matching.

``` r
data3$grp <- 1 

model <- clogistic(d ~ bcg + cage, strata = grp, data3) 

model
```

    ## 
    ## Call: 
    ## clogistic(formula = d ~ bcg + cage, strata = grp, data = data3)
    ## 
    ## 
    ## 
    ## 
    ##              coef exp(coef) se(coef)      z       p
    ## bcgPresent -1.104     0.332    0.193 -5.730 1.0e-08
    ## cage5-9     0.292     1.339    0.262  1.113 2.7e-01
    ## cage10-14   0.213     1.237    0.302  0.704 4.8e-01
    ## cage15-19   0.485     1.624    0.322  1.508 1.3e-01
    ## cage20-24   1.223     3.398    0.295  4.151 3.3e-05
    ## cage25-29   0.884     2.420    0.246  3.594 3.3e-04
    ## cage30-34   0.734     2.084    0.240  3.063 2.2e-03
    ## 
    ## Likelihood ratio test=107  on 7 df, p=0, n=1370

Here the matching variable **grp** is 1 for everyone, so the effect of
matching is neutral. In such case, both unconditional and conditional
logistic regression lead to approximately same results and matching can
be ignored.
