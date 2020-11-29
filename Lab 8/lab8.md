Statistical Methods in Epidemiology - Lab8
================

# Matched Case - Control Studies: The Salmonella Typhimurium dataset

On today’s practical, *typhi* dataset is going to be examined. It’s
about a matched case - control study, set up by the Danish Zoonosis
Centre in order to find the sources that led to a large number of
Salmonella Typhimurium cases in the fall of 1996 in Fyn County in
Denmark. Each case was matched to two controls of the same age, sex and
residency. The rest of the variables include participants’ food intake
during the last two weeks. The dataset lies in **SME** package, which is
created by us for the labs’ needs. It currently contains all the labs’
datasets, along with a few functions that help us present the results of
the analysis. You can download it by typing:
`devtools::install_github("argykar/SME")`(package **devtools** has to be
installed first). Now, we are going to load the libraries and the
dataset that will be used throughout the lab and create factors for the
categorical variables.

``` r
library(kableExtra)
library(SME)
library(survival)
library(tidyverse)
```

``` r
data("typhi")

typhi$sex <- factor(typhi$sex, labels = c("Male", "Female"))
typhi$beef <- factor(typhi$beef, labels = c("No", "Yes"))
typhi$pork <- factor(typhi$pork, labels = c("No", "Yes"))
typhi$poultry <- factor(typhi$poultry, labels = c("No", "Yes"))
typhi$liverp <- factor(typhi$liverp, labels = c("No", "Yes"))
typhi$meat <- factor(typhi$meat, labels = c("No", "Yes"))
typhi$veg <- factor(typhi$veg, labels = c("No", "Yes"))
typhi$fruit <- factor(typhi$fruit, labels = c("No", "Yes"))
typhi$egg <- factor(typhi$egg, labels = c("No", "Yes"))
typhi$plant7 <- factor(typhi$plant7, labels = c("No", "Yes"))
```

Here, we are combining the functions `table_stata()` and `display()`
from **SME** and creating the function *tbl*, so as to reduce the lines
of code in the code chunks that follow.

``` r
tbl <- function(x, exp = FALSE, digits = 3) {
  x %>%
    table_stata(exp = exp) %>%
    display(digits = digits)
}
```

First, we’ll examine the effect of pork using conditional logistic
regression.

``` r
fit_pork <- clogit(case ~ pork + strata(set), data = typhi)

fit_pork %>%
  tbl()
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

Pr(\>|z|)

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

porkYes

</td>

<td style="text-align:right;">

0.266

</td>

<td style="text-align:right;">

0.454

</td>

<td style="text-align:right;">

0.585

</td>

<td style="text-align:right;">

0.558

</td>

<td style="text-align:right;">

\-0.624

</td>

<td style="text-align:right;">

1.155

</td>

</tr>

</tbody>

</table>

From the table above, we can see that pork consumption is associated
with an average increase of the odds of salmonella by exp(0.2657) - 1 =
1.3043 - 1 = 30.43%.

Next, we’ll perform the same univariate analysis for the rest of the
variables. In order to avoid fitting manually the different models, we
are employing an automatic approach. Let’s analyze our steps. First, we
extract the names of the variables that will be used as explanatory in
the models and save them in a vector called *vars* (character vector).
Next, we are creating a formula for each variable of *vars*, using
`sapply()`. As *vars* is a character vector we need to use the function
`as.formula()` to convert in into a formula and as we also need both the
dependent (i.e *case*) and the strata variables (i.e *set*) to construct
the formula, we need to use quotation marks around them and use
`paste0()` to create a single character object from the arguments that
we pass into it. After the construction of the formulae, we have to fit
the model for each of them, so we use `lapply()`. What `bquote()` and
`eval()` functions do in this case, is not obvious. We could have
avoided them, as, in this case, our results would not be
distinguishable. However, let’s describe our idea. Every fitted model
has a *Call* component (e.g fit1$call). Had we not use `bquote()`, the
clogit call for each of the above variables would be
*clogit(as.formula(paste0(“case \~”, var, “+ strata(set)”)), data =
typhi)* and that is at least uninformative for someone that would notice
it. By using `bquote` the call, e.g for *veg* variable, becomes
*clogit(case \~ veg + strata(set), data = typhi)*. The `eval()`
function, essentially, executes the call (i.e fits the model), as
without it, we would have obtained only a list with calls.

``` r
vars <- names(typhi[c(4, 6:length(typhi))])

frms <- sapply(vars, function(var) {as.formula(paste0("case ~", var,  "+ strata(set)"))})
               
fits <- lapply(1:8, function(i) {eval(bquote(clogit(.(frms[[i]]), data = typhi)))})
```

Let’s extract the piece of information we want from the fitted models.

``` r
fits[[1]] %>%
  tbl()
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

Pr(\>|z|)

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

beefYes

</td>

<td style="text-align:right;">

\-0.119

</td>

<td style="text-align:right;">

0.422

</td>

<td style="text-align:right;">

\-0.281

</td>

<td style="text-align:right;">

0.778

</td>

<td style="text-align:right;">

\-0.947

</td>

<td style="text-align:right;">

0.709

</td>

</tr>

</tbody>

</table>

``` r
fits[[2]] %>%
  tbl()
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

Pr(\>|z|)

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

poultryYes

</td>

<td style="text-align:right;">

0.128

</td>

<td style="text-align:right;">

0.39

</td>

<td style="text-align:right;">

0.327

</td>

<td style="text-align:right;">

0.743

</td>

<td style="text-align:right;">

\-0.636

</td>

<td style="text-align:right;">

0.891

</td>

</tr>

</tbody>

</table>

``` r
fits[[3]] %>%
  tbl()
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

Pr(\>|z|)

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

liverpYes

</td>

<td style="text-align:right;">

\-1.792

</td>

<td style="text-align:right;">

1.155

</td>

<td style="text-align:right;">

\-1.552

</td>

<td style="text-align:right;">

0.121

</td>

<td style="text-align:right;">

\-4.055

</td>

<td style="text-align:right;">

0.471

</td>

</tr>

</tbody>

</table>

``` r
fits[[4]] %>%
  tbl()
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

Pr(\>|z|)

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

meatYes

</td>

<td style="text-align:right;">

0.613

</td>

<td style="text-align:right;">

0.463

</td>

<td style="text-align:right;">

1.324

</td>

<td style="text-align:right;">

0.185

</td>

<td style="text-align:right;">

\-0.295

</td>

<td style="text-align:right;">

1.521

</td>

</tr>

</tbody>

</table>

``` r
fits[[5]] %>%
  tbl()
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

Pr(\>|z|)

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

vegYes

</td>

<td style="text-align:right;">

\-0.63

</td>

<td style="text-align:right;">

0.687

</td>

<td style="text-align:right;">

\-0.917

</td>

<td style="text-align:right;">

0.359

</td>

<td style="text-align:right;">

\-1.977

</td>

<td style="text-align:right;">

0.717

</td>

</tr>

</tbody>

</table>

``` r
fits[[6]] %>%
  tbl()
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

Pr(\>|z|)

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

fruitYes

</td>

<td style="text-align:right;">

\-1.812

</td>

<td style="text-align:right;">

0.659

</td>

<td style="text-align:right;">

\-2.749

</td>

<td style="text-align:right;">

0.006

</td>

<td style="text-align:right;">

\-3.104

</td>

<td style="text-align:right;">

\-0.52

</td>

</tr>

</tbody>

</table>

``` r
fits[[7]] %>%
  tbl()
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

Pr(\>|z|)

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

eggYes

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0.53

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

\-1.039

</td>

<td style="text-align:right;">

1.039

</td>

</tr>

</tbody>

</table>

``` r
fits[[8]] %>%
  tbl()
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

Pr(\>|z|)

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

plant7Yes

</td>

<td style="text-align:right;">

1.497

</td>

<td style="text-align:right;">

0.519

</td>

<td style="text-align:right;">

2.883

</td>

<td style="text-align:right;">

0.004

</td>

<td style="text-align:right;">

0.479

</td>

<td style="text-align:right;">

2.515

</td>

</tr>

</tbody>

</table>

The univariate analysis indicate that only fruit consumption and meat
consumption from plant7 are potentially significant predictors. As
regards fruit consumption, we found out that it has a protective effect.
More specifically, the odds of salmonella is 1 - exp(-1.8119) =
-15.3337% lower for people who eat fruits compared with people who do
not include fruits in their daily nutrition.

Next, we will investigate whether an interaction between *fruit* and
*plant7* exists. We will fit a model that includes the main effects of
them, as well as one more including an interaction between them.

``` r
fit_me <- clogit(case ~ fruit + plant7 + strata(set), data = typhi)

fit_me %>%
  tbl()
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

Pr(\>|z|)

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

fruitYes

</td>

<td style="text-align:right;">

\-1.419

</td>

<td style="text-align:right;">

0.753

</td>

<td style="text-align:right;">

\-1.885

</td>

<td style="text-align:right;">

0.059

</td>

<td style="text-align:right;">

\-2.894

</td>

<td style="text-align:right;">

0.056

</td>

</tr>

<tr>

<td style="text-align:left;">

plant7Yes

</td>

<td style="text-align:right;">

1.496

</td>

<td style="text-align:right;">

0.586

</td>

<td style="text-align:right;">

2.554

</td>

<td style="text-align:right;">

0.011

</td>

<td style="text-align:right;">

0.348

</td>

<td style="text-align:right;">

2.645

</td>

</tr>

</tbody>

</table>

Adjusting for plant7, the protective effect of fruit consumption
decreases both in strength and in significance.

``` r
fit_int <- clogit(case ~ fruit*plant7 + strata(set), data = typhi)

fit_int %>%
  tbl()
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

Pr(\>|z|)

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

fruitYes

</td>

<td style="text-align:right;">

\-1.118

</td>

<td style="text-align:right;">

1.116

</td>

<td style="text-align:right;">

\-1.002

</td>

<td style="text-align:right;">

0.317

</td>

<td style="text-align:right;">

\-3.306

</td>

<td style="text-align:right;">

1.070

</td>

</tr>

<tr>

<td style="text-align:left;">

plant7Yes

</td>

<td style="text-align:right;">

2.012

</td>

<td style="text-align:right;">

1.558

</td>

<td style="text-align:right;">

1.291

</td>

<td style="text-align:right;">

0.197

</td>

<td style="text-align:right;">

\-1.042

</td>

<td style="text-align:right;">

5.066

</td>

</tr>

<tr>

<td style="text-align:left;">

fruitYes:plant7Yes

</td>

<td style="text-align:right;">

\-0.588

</td>

<td style="text-align:right;">

1.632

</td>

<td style="text-align:right;">

\-0.360

</td>

<td style="text-align:right;">

0.719

</td>

<td style="text-align:right;">

\-3.787

</td>

<td style="text-align:right;">

2.612

</td>

</tr>

</tbody>

</table>

Exponentiating the summation of all coefficients, we obtain the odds
ratio comparing people with fruit and meat consumption from plant 7 with
those that neither eat fruits nor consume meat from plant7. As *fruit*
is a protective predictor and *plant7* is a harmful one, this odds ratio
is meaningless. We’d like to compare people that eat fruits and do not
consume meat from plant 7 to those that do not eat fruit and consume
meat from plant7. Consequently, we’ll change the reference category of
*plant7* and fit the interaction model as before.

``` r
typhi$plant7 <- fct_relevel(typhi$plant7, "Yes")

fit_inter <- clogit(case ~ fruit*plant7 + strata(set), data = typhi)

fit_inter %>%
  tbl()
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

Pr(\>|z|)

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

fruitYes

</td>

<td style="text-align:right;">

\-1.706

</td>

<td style="text-align:right;">

1.138

</td>

<td style="text-align:right;">

\-1.499

</td>

<td style="text-align:right;">

0.134

</td>

<td style="text-align:right;">

\-3.936

</td>

<td style="text-align:right;">

0.525

</td>

</tr>

<tr>

<td style="text-align:left;">

plant7No

</td>

<td style="text-align:right;">

\-2.012

</td>

<td style="text-align:right;">

1.558

</td>

<td style="text-align:right;">

\-1.291

</td>

<td style="text-align:right;">

0.197

</td>

<td style="text-align:right;">

\-5.066

</td>

<td style="text-align:right;">

1.042

</td>

</tr>

<tr>

<td style="text-align:left;">

fruitYes:plant7No

</td>

<td style="text-align:right;">

0.588

</td>

<td style="text-align:right;">

1.632

</td>

<td style="text-align:right;">

0.360

</td>

<td style="text-align:right;">

0.719

</td>

<td style="text-align:right;">

\-2.612

</td>

<td style="text-align:right;">

3.787

</td>

</tr>

</tbody>

</table>

From the table above, we conclude that the odds of salmonella is on
average 1 - exp(-1.7059 -2.0115 + 0.5878) = 96% lower for fruit
consumers who do not eat meat originated from plant7 in comparison to
non-fruit consumers who eat meat from plant7 (the interaction term,
however, is not significant).

Lastly, we are considering a model with *fruit* and *plant7* where
reference category is the people that do not consume meat from plant7
independently of the fruit consumption. In this way, we can investigate
whether the effect of *fruit* is protective only among subjects that
have eaten meat from plant 7. So, we are creating a variable with the
following coding: 1 for those that have not eaten meat from plant7
(independently of fruit consumption), 2 for those that have eaten meat
from plant7 and are not fruit consumers and 3 for those that have eaten
meat from plant7 and are also fruit consumers. The model with this
variable as explanatory is not an interaction model. Before we proceed,
we are relevelling *plant7* variable once more.

``` r
typhi$plant7 <- fct_relevel(typhi$plant7, "No")

typhi <- typhi %>%
  mutate(fp = factor(case_when(plant7 == "No" ~ 1,
                               plant7 == "Yes" & fruit == "No" ~ 2,
                               plant7 == "Yes" & fruit == "Yes" ~ 3)))

model <- clogit(case ~ fp + strata(set), data = typhi)

model %>%
  tbl()
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

Pr(\>|z|)

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

fp2

</td>

<td style="text-align:right;">

3.055

</td>

<td style="text-align:right;">

1.170

</td>

<td style="text-align:right;">

2.612

</td>

<td style="text-align:right;">

0.009

</td>

<td style="text-align:right;">

0.762

</td>

<td style="text-align:right;">

5.347

</td>

</tr>

<tr>

<td style="text-align:left;">

fp3

</td>

<td style="text-align:right;">

1.304

</td>

<td style="text-align:right;">

0.593

</td>

<td style="text-align:right;">

2.200

</td>

<td style="text-align:right;">

0.028

</td>

<td style="text-align:right;">

0.142

</td>

<td style="text-align:right;">

2.466

</td>

</tr>

</tbody>

</table>

We conclude that those who eat meat from plant7 and are not fruit
consumers have on average exp(3.055) = 21.21 times higher odds of
salmonella compared with people who have not eaten meat from plant7
(independently of whether they eat fruits or not), while the odds ratio
of disease comparing those that eat fruits and have eaten meat from
plant7 to those who do not eat meat from plant7 (independently of fruit
consumption) is exp(1.304) = 3.68. It’s a strong indication that the
aggravating effect of plant7 is remarkably limited for fruit consumers.
