Statistical Methods in Epidemiology
================
13 Νοέμβριος, 2020

<p style="text-align:center";>
  <b> 
    <font size="6">Lab 3</font>
  </b>
</p>

## Table of Contents

1.  [Analyzing Rates](#rates)
2.  [Exposures with More Than Two Levels](#exposure)
3.  [Controlling for a Confounding Variable](#confounding)
4.  [Likelihood Ratio and Wald Tests](#wald)
5.  [Metric Exposures and Linear Effects](#explin)
6.  [Testing for Linearity](#linearity)
7.  [Interaction](#interaction)

In the practical that follows, *diet* data set, which contains data from
a pilot study of 337 men who kept a record of their fully weighted diet
over two weeks, will be analyzed. The particular data set lies in
**biostat3** package, so we will load it along with the rest of packages
that are going to be used.

``` r
library(biostat3)
library(tidyverse)
library(kableExtra)
library(lmtest)
```

A brief investigation of our data reveals that variable **energy** is
probably measured in calories, so we will convert it into kcals in order
to obtain the same results with the respective **stata** practical.
Also, we will create the categorical variables **eng3** and **htgrp**
that were also created in the previous lab and will be used later on.

``` r
data(diet)

diet <- diet %>%
  mutate(energy = energy / 1000, 
         eng3 = cut(energy, breaks = c(1.5, 2.5, 3, 4.5), labels = c("0", "1", "2")), 
         htgrp = cut(height, breaks = c(150, 170, 175, 180, 195), labels = c("0", "1", "2", "3")))
```

Preview of the variables created below:

<table>

<thead>

<tr>

<th style="text-align:right;">

energy

</th>

<th style="text-align:left;">

eng3

</th>

<th style="text-align:left;">

htgrp

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

2.02

</td>

<td style="text-align:left;">

0

</td>

<td style="text-align:left;">

1

</td>

</tr>

<tr>

<td style="text-align:right;">

2.45

</td>

<td style="text-align:left;">

0

</td>

<td style="text-align:left;">

2

</td>

</tr>

<tr>

<td style="text-align:right;">

2.28

</td>

<td style="text-align:left;">

0

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:right;">

2.47

</td>

<td style="text-align:left;">

0

</td>

<td style="text-align:left;">

0

</td>

</tr>

<tr>

<td style="text-align:right;">

2.36

</td>

<td style="text-align:left;">

0

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:right;">

2.38

</td>

<td style="text-align:left;">

0

</td>

<td style="text-align:left;">

0

</td>

</tr>

</tbody>

</table>

In what follows, we will mainly work with Poisson regression models. As,
on most occasions, we prefer the <u>incidence rate ratios</u> to effects
on a log scale, we will create a function that accepts a Poisson model
as argument and returns the incidence rate ratio for each of the
independent variables along with their standard error (calculated using
delta method), p-value and 95% confidence intervals (see also **mfx**
package for something similar).

<p style="color:red";><b>Don't be discouraged if you don't understand. The following requires enough practice with R and programming concepts in general. This function is for the authors' ease to continue with the Lab. The main idea for the creation of a function is to not repeat the same code chunks across the script.</b></p>

``` r
r_table <- function(model) {
  if (is(model, "glm") == TRUE) {
    if (family(model)$family == "poisson") {
      modification <- 
        model %>% 
        summary() %>% 
        coefficients() %>% 
        as.data.frame() %>% 
        slice(-1) %>% 
        rownames_to_column(var = "Predictors") %>% 
        mutate(IRR = exp(Estimate), 
               Std_Err = IRR * `Std. Error`, 
               pval = sprintf("%1.3f", `Pr(>|z|)`),
               pvalue = ifelse(pval < 0.001, '<0.001', pval)) %>%
        select(Predictors, 
               IRR, Std_Err, 
               "z_value" = `z value`, 
               "P-Value" = pvalue)
      
      if (nrow(modification) > 1) {
        final_table <- bind_cols(modification, 
                                 LL = exp(confint.default(model)[-1, ])[, 1], 
                                 UL = exp(confint.default(model)[-1, ])[, 2])
        } else {
          final_table <- bind_cols(modification, 
                                   LL = exp(confint.default(model)[-1, ])[1], 
                                   UL = exp(confint.default(model)[-1, ])[2])
        }
      } else if (family(model)$family == "binomial") {
              modification <- 
        model %>% 
        summary() %>% 
        coefficients() %>% 
        as.data.frame() %>% 
        slice(-1) %>% 
        rownames_to_column(var = "Predictors") %>% 
        mutate(OR = exp(Estimate), 
               Std_Err = OR * `Std. Error`, 
               pval = sprintf("%1.3f", `Pr(>|z|)`),
               pvalue = ifelse(pval < 0.001, '<0.001', pval)) %>%
        select(Predictors, 
               OR, Std_Err, 
               "z_value" = `z value`, 
               "P-Value" = pvalue)
      
      if (nrow(modification) > 1) {
        final_table <- bind_cols(modification, 
                                 LL = exp(confint.default(model)[-1, ])[, 1], 
                                 UL = exp(confint.default(model)[-1, ])[, 2])
        } else {
          final_table <- bind_cols(modification, 
                                   LL = exp(confint.default(model)[-1, ])[1], 
                                   UL = exp(confint.default(model)[-1, ])[2]) 
        } 
      } else {
        stop("Family of the model is not Poisson or Binomial")
        }
    } else {
      stop("Object is not of class glm")
      } 
    return(final_table) 
}
```

## I. Analyzing Rates <a name="rates"></a>

We will investigate the effect of energy intake on the rate of cases
with CHD fitting a poisson regression model with `glm()` function,
specifying the argument **family = “poisson”** (by default, when
**family = “poisson”** is specified, link function is the natural
logarithm). The version of energy intake that we will originally use is
the categorical one (variable **hieng**). Note, that person-years at
risk for each subject are already in our data set (variable **y**) and
we will introduce them as an offset into our model (argument
**offset**). Estimated coefficients along with further useful
information about the model can be provided with Base R function
`summary()` or `summ()` from **jtools** package. We will also use the
function that we created above in order to obtain the incidence rate
ratio of CHD comparing men with high energy intake to those with low.

``` r
fit1 <- glm(chd ~ hieng, data = diet, family = "poisson", offset = log(y))

r_table(fit1) %>%
  kbl(digits = 2) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Predictors

</th>

<th style="text-align:right;">

IRR

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

hienghigh

</td>

<td style="text-align:right;">

0.52

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:right;">

\-2.16

</td>

<td style="text-align:left;">

0.031

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:right;">

0.94

</td>

</tr>

</tbody>

</table>

## II. Exposures with More Than Two Levels <a name="exposure"></a>

As before, we ’ll investigate the effect of energy intake on the rate of
cases with CHD. This time, however, we ’ll use the second categorical
version of energy, **eng3**.

``` r
fit2 <- glm(chd ~ eng3, data = diet, family = "poisson", offset = log(y))

r_table(fit2) %>%
  kbl(digits = 2) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Predictors

</th>

<th style="text-align:right;">

IRR

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

eng31

</td>

<td style="text-align:right;">

0.65

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:right;">

\-1.33

</td>

<td style="text-align:left;">

0.182

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:right;">

1.23

</td>

</tr>

<tr>

<td style="text-align:left;">

eng32

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

\-2.87

</td>

<td style="text-align:left;">

0.004

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

0.67

</td>

</tr>

</tbody>

</table>

## III. Controlling for a Confounding Variable <a name="confounding"></a>

Now, we’ll find the effect of high energy (variable **hieng**)
controlled for **job**, as well as the effect of **eng3** controlled
for:

1.  **job**
2.  **job** and **htgrp**.

<!-- end list -->

``` r
fit3 <- glm(chd ~ hieng + job, data = diet, family = "poisson", offset = log(y))

r_table(fit3) %>%
  kbl(digits = 2) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Predictors

</th>

<th style="text-align:right;">

IRR

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

hienghigh

</td>

<td style="text-align:right;">

0.52

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:right;">

\-2.13

</td>

<td style="text-align:left;">

0.033

</td>

<td style="text-align:right;">

0.29

</td>

<td style="text-align:right;">

0.95

</td>

</tr>

<tr>

<td style="text-align:left;">

jobconductor

</td>

<td style="text-align:right;">

1.36

</td>

<td style="text-align:right;">

0.53

</td>

<td style="text-align:right;">

0.78

</td>

<td style="text-align:left;">

0.436

</td>

<td style="text-align:right;">

0.63

</td>

<td style="text-align:right;">

2.94

</td>

</tr>

<tr>

<td style="text-align:left;">

jobbank

</td>

<td style="text-align:right;">

0.88

</td>

<td style="text-align:right;">

0.32

</td>

<td style="text-align:right;">

\-0.34

</td>

<td style="text-align:left;">

0.736

</td>

<td style="text-align:right;">

0.43

</td>

<td style="text-align:right;">

1.81

</td>

</tr>

</tbody>

</table>

``` r
fit4 <- glm(chd ~ eng3 + job, data = diet, family = "poisson", offset = log(y))

r_table(fit4) %>%
  kbl(digits = 2) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Predictors

</th>

<th style="text-align:right;">

IRR

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

eng31

</td>

<td style="text-align:right;">

0.63

</td>

<td style="text-align:right;">

0.21

</td>

<td style="text-align:right;">

\-1.42

</td>

<td style="text-align:left;">

0.157

</td>

<td style="text-align:right;">

0.33

</td>

<td style="text-align:right;">

1.20

</td>

</tr>

<tr>

<td style="text-align:left;">

eng32

</td>

<td style="text-align:right;">

0.28

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

\-2.90

</td>

<td style="text-align:left;">

0.004

</td>

<td style="text-align:right;">

0.12

</td>

<td style="text-align:right;">

0.67

</td>

</tr>

<tr>

<td style="text-align:left;">

jobconductor

</td>

<td style="text-align:right;">

1.46

</td>

<td style="text-align:right;">

0.58

</td>

<td style="text-align:right;">

0.96

</td>

<td style="text-align:left;">

0.336

</td>

<td style="text-align:right;">

0.67

</td>

<td style="text-align:right;">

3.17

</td>

</tr>

<tr>

<td style="text-align:left;">

jobbank

</td>

<td style="text-align:right;">

0.92

</td>

<td style="text-align:right;">

0.34

</td>

<td style="text-align:right;">

\-0.22

</td>

<td style="text-align:left;">

0.828

</td>

<td style="text-align:right;">

0.45

</td>

<td style="text-align:right;">

1.89

</td>

</tr>

</tbody>

</table>

``` r
fit5 <- glm(chd ~ eng3 + job + htgrp, data = diet, family = "poisson", offset = log(y))

r_table(fit5) %>%
  kbl(digits = 2) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Predictors

</th>

<th style="text-align:right;">

IRR

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

eng31

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

0.23

</td>

<td style="text-align:right;">

\-1.08

</td>

<td style="text-align:left;">

0.282

</td>

<td style="text-align:right;">

0.36

</td>

<td style="text-align:right;">

1.35

</td>

</tr>

<tr>

<td style="text-align:left;">

eng32

</td>

<td style="text-align:right;">

0.35

</td>

<td style="text-align:right;">

0.16

</td>

<td style="text-align:right;">

\-2.37

</td>

<td style="text-align:left;">

0.018

</td>

<td style="text-align:right;">

0.15

</td>

<td style="text-align:right;">

0.84

</td>

</tr>

<tr>

<td style="text-align:left;">

jobconductor

</td>

<td style="text-align:right;">

1.25

</td>

<td style="text-align:right;">

0.50

</td>

<td style="text-align:right;">

0.55

</td>

<td style="text-align:left;">

0.584

</td>

<td style="text-align:right;">

0.57

</td>

<td style="text-align:right;">

2.73

</td>

</tr>

<tr>

<td style="text-align:left;">

jobbank

</td>

<td style="text-align:right;">

1.40

</td>

<td style="text-align:right;">

0.53

</td>

<td style="text-align:right;">

0.89

</td>

<td style="text-align:left;">

0.371

</td>

<td style="text-align:right;">

0.67

</td>

<td style="text-align:right;">

2.92

</td>

</tr>

<tr>

<td style="text-align:left;">

htgrp1

</td>

<td style="text-align:right;">

0.70

</td>

<td style="text-align:right;">

0.25

</td>

<td style="text-align:right;">

\-1.02

</td>

<td style="text-align:left;">

0.308

</td>

<td style="text-align:right;">

0.35

</td>

<td style="text-align:right;">

1.40

</td>

</tr>

<tr>

<td style="text-align:left;">

htgrp2

</td>

<td style="text-align:right;">

0.59

</td>

<td style="text-align:right;">

0.24

</td>

<td style="text-align:right;">

\-1.29

</td>

<td style="text-align:left;">

0.197

</td>

<td style="text-align:right;">

0.27

</td>

<td style="text-align:right;">

1.32

</td>

</tr>

<tr>

<td style="text-align:left;">

htgrp3

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

\-0.02

</td>

<td style="text-align:left;">

0.988

</td>

<td style="text-align:right;">

0.00

</td>

<td style="text-align:right;">

Inf

</td>

</tr>

</tbody>

</table>

## IV. Likelihood Ratio and Wald Tests <a name="wald"></a>

We ’ll perform a likelihood ratio test to test the hypothesis that
coefficients of **eng3** from model *fit2* are zero. We ’ll present two
ways:

1.  the manual using `logLik()` function, from Base R, that calculates
    the log-likelihood of a model
2.  the automatic using `lrtest()` function from **lmtest** package.

In both cases, initially, we need to create the reduced model (i.e. the
model without the variables whose coefficients are tested to be equal to
zero). We ’ll accomplish it using `update()` function.

``` r
fit6 <- update(fit2, .~.-eng3)

pchisq(q = -2 * (as.numeric(logLik(fit6)) - as.numeric(logLik(fit2))), df = 2, lower.tail = FALSE)
```

    ## [1] 0.01

``` r
lrtest(fit2, fit6) %>%
  kbl(digits = 3) %>%
  kable_styling()
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

Pr(\>Chisq)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

\-173

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

1

</td>

<td style="text-align:right;">

\-177

</td>

<td style="text-align:right;">

\-2

</td>

<td style="text-align:right;">

9.2

</td>

<td style="text-align:right;">

0.01

</td>

</tr>

</tbody>

</table>

An approximation of likelihood ratio test, the Wald test, can be
implemented with `waldtest()` function, again from **lmtest** package,
specifying the argument **test = “Chisq”**.

``` r
waldtest(fit2, fit6, test = "Chisq") %>%
  kbl(digits = 4) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:right;">

Res.Df

</th>

<th style="text-align:right;">

Df

</th>

<th style="text-align:right;">

Chisq

</th>

<th style="text-align:right;">

Pr(\>Chisq)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

334

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

336

</td>

<td style="text-align:right;">

\-2

</td>

<td style="text-align:right;">

8.24

</td>

<td style="text-align:right;">

0.0162

</td>

</tr>

</tbody>

</table>

5.  Same as before, but this time we create a model with **hieng**.

<!-- end list -->

``` r
hiengfit <- glm(chd ~ hieng + job, data = diet, family = "poisson", offset = log(y))

nohiengfit <- glm(chd ~  job, data = diet, family = "poisson", offset = log(y))

lrtest(hiengfit, nohiengfit) %>%
  kbl(digits = 3) %>%
  kable_styling()
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

Pr(\>Chisq)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

\-174

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

\-177

</td>

<td style="text-align:right;">

\-1

</td>

<td style="text-align:right;">

4.69

</td>

<td style="text-align:right;">

0.03

</td>

</tr>

</tbody>

</table>

6.  Testing the addition of **eng3** in the model that contains already
    the **job**.

<!-- end list -->

``` r
eng3fit <- glm(chd ~ eng3 + job, data = diet, family = "poisson", offset = log(y))

noeng3fit <- glm(chd ~ job, data = diet, family = "poisson", offset = log(y))

lrtest(eng3fit, noeng3fit) %>%
  kbl(digits = 3) %>%
  kable_styling()
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

Pr(\>Chisq)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

5

</td>

<td style="text-align:right;">

\-172

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

\-177

</td>

<td style="text-align:right;">

\-2

</td>

<td style="text-align:right;">

9.33

</td>

<td style="text-align:right;">

0.009

</td>

</tr>

</tbody>

</table>

## V. Metric Exposures and Linear Effects <a name="explin"></a>

``` r
fit8 <- glm(chd ~ height, data = diet, family = "poisson", offset = log(y))

fit8 %>%
  r_table() %>%
  kbl(digits = 3) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Predictors

</th>

<th style="text-align:right;">

IRR

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

height

</td>

<td style="text-align:right;">

0.912

</td>

<td style="text-align:right;">

0.02

</td>

<td style="text-align:right;">

\-4.21

</td>

<td style="text-align:left;">

\<0.001

</td>

<td style="text-align:right;">

0.873

</td>

<td style="text-align:right;">

0.952

</td>

</tr>

</tbody>

</table>

For the model with **fat** as predictor, we will use the *diet.dta*
dataset that we used in class. The `*biostat3::diet*` does not include
fat as a variable in the dataset.

``` r
diet2 <- haven::read_dta('diet.dta')

# The diet2 dataset does not contain the person-years, so we have to create it.
diet2 <- diet2 %>%
  mutate(y = lubridate::time_length(difftime(dox, doe), 'years'))

fatfit <- glm(chd ~ fat, data = diet2, family = "poisson", offset = log(y))

fatfit %>%
  r_table() %>%
  kbl(digits = 4) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Predictors

</th>

<th style="text-align:right;">

IRR

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

fat

</td>

<td style="text-align:right;">

0.98

</td>

<td style="text-align:right;">

0.0069

</td>

<td style="text-align:right;">

\-2.89

</td>

<td style="text-align:left;">

0.004

</td>

<td style="text-align:right;">

0.966

</td>

<td style="text-align:right;">

0.994

</td>

</tr>

</tbody>

</table>

For every 1 g/day increase in fat we expect on average 2% decrease in
the CHD mortality rate. But what if we want the 10 g/day increase? All
we have to do is isolate the IRR and raise it to the tenth power.

``` r
fatfit %>%
  r_table() %>%
  .[['IRR']] %>%
  . ^ 10
```

    ## [1] 0.816

So, for every 10 g/day increase in fat we expect on average 18.4%
decrease in CHD mortality rate.

## VI. Testing for Linearity <a name="linearity"></a>

First of, we are going to create the **htsq** variable, that is defined
as ![htsq =
height^2](https://latex.codecogs.com/png.latex?%5Cdpi{150}htsq%20%3D%20height%5E2
"htsq = height^2"). Next, we will proceed with the model that contains
**htsq** and **height** as predictors.

``` r
diet2 <- diet2 %>%
  mutate(htsq = height ^ 2)
```

``` r
htsqfit <- glm(chd ~ height + htsq, data = diet2, family = "poisson", offset = log(y))

# Extract the coefficients of the model above
res <- htsqfit %>% 
  summary() %>% 
  .$coef

# Combine it with the 95% CIs
cbind(res, confint.default(htsqfit)) %>% 
  kbl(digits = 3) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

Estimate

</th>

<th style="text-align:right;">

Std. Error

</th>

<th style="text-align:right;">

z value

</th>

<th style="text-align:right;">

Pr(\>|z|)

</th>

<th style="text-align:right;">

2.5 %

</th>

<th style="text-align:right;">

97.5 %

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

(Intercept)

</td>

<td style="text-align:right;">

\-27.832

</td>

<td style="text-align:right;">

73.954

</td>

<td style="text-align:right;">

\-0.376

</td>

<td style="text-align:right;">

0.707

</td>

<td style="text-align:right;">

\-172.779

</td>

<td style="text-align:right;">

117.115

</td>

</tr>

<tr>

<td style="text-align:left;">

height

</td>

<td style="text-align:right;">

0.370

</td>

<td style="text-align:right;">

0.873

</td>

<td style="text-align:right;">

0.424

</td>

<td style="text-align:right;">

0.672

</td>

<td style="text-align:right;">

\-1.341

</td>

<td style="text-align:right;">

2.080

</td>

</tr>

<tr>

<td style="text-align:left;">

htsq

</td>

<td style="text-align:right;">

\-0.001

</td>

<td style="text-align:right;">

0.003

</td>

<td style="text-align:right;">

\-0.530

</td>

<td style="text-align:right;">

0.596

</td>

<td style="text-align:right;">

\-0.006

</td>

<td style="text-align:right;">

0.004

</td>

</tr>

</tbody>

</table>

Notice that from the Wald test we get ![z ^ 2 = 0.281, p = 0.596
\> 0.05](https://latex.codecogs.com/png.latex?%5Cdpi{150}z%20%5E%202%20%3D%200.281%2C%20p%20%3D%200.596%20%3E%200.05
"z ^ 2 = 0.281, p = 0.596 \> 0.05") for the **htsq** term.

Another way of doing the same as above:

Combining *fit8* and *htsqfit*, we want to evaluate the addition of the
squared term with the Likelihood ratio test.

``` r
lrtest(fit8, htsqfit) %>%
  kbl(digits = 3) %>%
  kable_styling()
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

Pr(\>Chisq)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

\-164

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

\-164

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0.301

</td>

<td style="text-align:right;">

0.583

</td>

</tr>

</tbody>

</table>

From the table above, we observe that ![p = 0.583
\> 0.05](https://latex.codecogs.com/png.latex?%5Cdpi{150}p%20%3D%200.583%20%3E%200.05
"p = 0.583 \> 0.05"), hence the squared term for the height is not
needed.

## VII. Interaction <a name="interaction"></a>

We are going to investigate the interaction between **hieng** and
**job** and its statistical significance.

``` r
inter1 <- glm(chd ~ job * hieng, data = diet, family = "poisson", offset = log(y))

nointer1 <- glm(chd ~ job + hieng, data = diet, family = "poisson", offset = log(y))

lrtest(inter1, nointer1) %>%
  kbl(digits = 3) %>%
  kable_styling()
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

Pr(\>Chisq)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

\-174

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

4

</td>

<td style="text-align:right;">

\-174

</td>

<td style="text-align:right;">

\-2

</td>

<td style="text-align:right;">

0.333

</td>

<td style="text-align:right;">

0.847

</td>

</tr>

</tbody>

</table>

We conclude that, the addition of the interaction terms is not
statistically significant in the 5% confidence level since ![X^2 = 0.33
(df = 2), p
= 0.847](https://latex.codecogs.com/png.latex?%5Cdpi{150}X%5E2%20%3D%200.33%20%28df%20%3D%202%29%2C%20p%20%3D%200.847
"X^2 = 0.33 (df = 2), p = 0.847").

Hence, the effect of high energy on the rates of CHD mortality are not
different within the different levels of job.

Next, we will investigate the interaction between **hieng** and
**htgrp**.

``` r
inter2 <- glm(chd ~ hieng * htgrp, data = diet, family = "poisson", offset = log(y)) 

nointer2 <-  glm(chd ~ hieng + htgrp, data = diet, family = "poisson", offset = log(y))

lrtest(inter2, nointer2) %>%
  kbl(digits = 3) %>%
  kable_styling()
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

Pr(\>Chisq)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

8

</td>

<td style="text-align:right;">

\-159

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

5

</td>

<td style="text-align:right;">

\-159

</td>

<td style="text-align:right;">

\-3

</td>

<td style="text-align:right;">

0.947

</td>

<td style="text-align:right;">

0.814

</td>

</tr>

</tbody>

</table>

From the likelihood ratio test, we infer that ![X^2 = 0.947(df = 3), p
= 0.814](https://latex.codecogs.com/png.latex?%5Cdpi{150}X%5E2%20%3D%200.947%28df%20%3D%203%29%2C%20p%20%3D%200.814
"X^2 = 0.947(df = 3), p = 0.814")

The above results, show that the effect of high energy on the rates of
CHD mortality are not different within the different levels of height
groups.

Finally, we are going to test the interaction between **hieng** and
**height**.

``` r
inter3 <- glm(chd ~ hieng * height, data = diet, family = "poisson", offset = log(y)) 

nointer3 <- glm(chd ~ hieng + height, data = diet, family = "poisson", offset = log(y)) 

lrtest(inter3, nointer3) %>%
  kbl(digits = 3) %>%
  kable_styling()
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

Pr(\>Chisq)

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

4

</td>

<td style="text-align:right;">

\-162

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

\-163

</td>

<td style="text-align:right;">

\-1

</td>

<td style="text-align:right;">

0.388

</td>

<td style="text-align:right;">

0.533

</td>

</tr>

</tbody>

</table>

Same as before, we do not have enough enough evidence to conclude that
interaction between the variables mentioned above is statistically
significant in the 5% confidence level, since ![X^2 = 0.388(df = 1), p
= 0.533](https://latex.codecogs.com/png.latex?%5Cdpi{150}X%5E2%20%3D%200.388%28df%20%3D%201%29%2C%20p%20%3D%200.533
"X^2 = 0.388(df = 1), p = 0.533").
