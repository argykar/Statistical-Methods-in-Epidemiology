Statistical Methods in Epidemiology - Lab9
================

# Interactions

In the *hyper.dta* file, data of 328 records of individuals are stored.
The outcome of interest is hypertension (variable *hyper*). We want to
examine whether the development of hypertension is associated with a
person’s BMI and a score indicating the dedication to Mediterranean
diet, categorised in two groups. Most importantly, we are interested in
the examination of the interaction of these two characteristics. We’ll
start by importing the libraries we need, reading the dataset and
creating factors for the categorical variables.

``` r
library(haven)
library(interactionR)
library(kableExtra)
library(nlWaldTest)
library(tidyverse)
library(SME)
```

``` r
data <- read_dta("hyper.dta")

data$hyper <- factor(data$hyper, labels = c("No", "Yes"))
data$bmi <- factor(data$bmi, labels = c("<=30", ">30"))
data$nomd <- factor(data$nomd, labels = c(">=5", "<5"))
```

## Interactions - Approach 1

We’ll investigate the interaction between BMI and Mediterranean diet by
creating a variable from the combination of the levels of these two
variables and introducing it in a logistic regression model.

``` r
data$bmi_nomd <- with(data, interaction(bmi, nomd))

fit1 <- glm(hyper ~ bmi_nomd, data, family = binomial)

r_table(fit1) %>%
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

bmi\_nomd\>30.\>=5

</td>

<td style="text-align:right;">

2.611

</td>

<td style="text-align:right;">

1.982

</td>

<td style="text-align:right;">

1.264

</td>

<td style="text-align:left;">

0.206

</td>

<td style="text-align:right;">

0.589

</td>

<td style="text-align:right;">

11.565

</td>

</tr>

<tr>

<td style="text-align:left;">

bmi\_nomd\<=30.\<5

</td>

<td style="text-align:right;">

1.547

</td>

<td style="text-align:right;">

1.463

</td>

<td style="text-align:right;">

0.461

</td>

<td style="text-align:left;">

0.645

</td>

<td style="text-align:right;">

0.242

</td>

<td style="text-align:right;">

9.878

</td>

</tr>

<tr>

<td style="text-align:left;">

bmi\_nomd\>30.\<5

</td>

<td style="text-align:right;">

3.787

</td>

<td style="text-align:right;">

2.977

</td>

<td style="text-align:right;">

1.694

</td>

<td style="text-align:left;">

0.090

</td>

<td style="text-align:right;">

0.811

</td>

<td style="text-align:right;">

17.681

</td>

</tr>

</tbody>

</table>

We’ll calculate and examine three indices for interaction. More
specifically:

  - Relative Excess Risk for Interaction and Synergy Index (interaction
    on additive scale)
  - Ratio of Odds Ratios (interaction on multiplicative scale)

They are calculated as:

  
![&#10;RERI = OR\_{A=1, B=1} - OR\_{A=1, B=0} - OR\_{A=0, B=1}
+ 1&#10;SI = \\frac{OR\_{A=1, B=1} - 1}{OR\_{A=1, B=0} + OR\_{A=0, B=1}
- 2}&#10;ROR = \\frac{OR\_{A=1, B=1}}{OR\_{A=1, B=0}\*OR\_{A=0,
B=1}}&#10;](https://latex.codecogs.com/png.latex?%5Cdpi{100}%0ARERI%20%3D%20OR_%7BA%3D1%2C%20B%3D1%7D%20-%20OR_%7BA%3D1%2C%20B%3D0%7D%20-%20OR_%7BA%3D0%2C%20B%3D1%7D%20%2B%201%0ASI%20%3D%20%5Cfrac%7BOR_%7BA%3D1%2C%20B%3D1%7D%20-%201%7D%7BOR_%7BA%3D1%2C%20B%3D0%7D%20%2B%20OR_%7BA%3D0%2C%20B%3D1%7D%20-%202%7D%0AROR%20%3D%20%5Cfrac%7BOR_%7BA%3D1%2C%20B%3D1%7D%7D%7BOR_%7BA%3D1%2C%20B%3D0%7D%2AOR_%7BA%3D0%2C%20B%3D1%7D%7D%0A
"
RERI = OR_{A=1, B=1} - OR_{A=1, B=0} - OR_{A=0, B=1} + 1
SI = \\frac{OR_{A=1, B=1} - 1}{OR_{A=1, B=0} + OR_{A=0, B=1} - 2}
ROR = \\frac{OR_{A=1, B=1}}{OR_{A=1, B=0}*OR_{A=0, B=1}}
")  

They are, essentialy, nonlinear transformations of the coefficients of
*fit1*, so we are going to use `nlConfint()` function from
**nlWaldTest** package to calculate 95% approximate confidence intervals
for them. As *SI* and *ROR* are ratios, we have to ensure that their
bounds are positive, so we’ll first calculate 95% approximate confidence
intervals for their natural logarithm and then we’ll exponentiate them
in order to obtain the appropriate confidence intervals.

``` r
reri <- as.data.frame(nlConfint(fit1, c("exp(b[4]) - exp(b[3]) - exp(b[2]) + 1")))
si <- as.data.frame(nlConfint(fit1, c("log(exp(b[4]) - 1) - log(exp(b[3]) + exp(b[2]) - 2)")))
ror <- as.data.frame(nlConfint(fit1, c("log(exp(b[4]) / (exp(b[3])*exp(b[2])))")))

data.frame(Index = c("RERI", "SI", "ROR"),
           Estimate = c(reri$value, exp(si$value), exp(ror$value)),
           Lower_L = c(reri$`2.5 %`, exp(si$`2.5 %`), exp(ror$`2.5 %`)),
           Upper_L = c(reri$`97.5 %`, exp(si$`97.5 %`), exp(ror$`97.5 %`))) %>%
  kbl(digits = 4) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Index

</th>

<th style="text-align:right;">

Estimate

</th>

<th style="text-align:right;">

Lower\_L

</th>

<th style="text-align:right;">

Upper\_L

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

RERI

</td>

<td style="text-align:right;">

0.6293

</td>

<td style="text-align:right;">

\-2.4824

</td>

<td style="text-align:right;">

3.7409

</td>

</tr>

<tr>

<td style="text-align:left;">

SI

</td>

<td style="text-align:right;">

1.2916

</td>

<td style="text-align:right;">

0.2872

</td>

<td style="text-align:right;">

5.8096

</td>

</tr>

<tr>

<td style="text-align:left;">

ROR

</td>

<td style="text-align:right;">

0.9377

</td>

<td style="text-align:right;">

0.1284

</td>

<td style="text-align:right;">

6.8474

</td>

</tr>

</tbody>

</table>

We conclude that:

  - The combined effect (relative excess risk) of low MDScale and high
    BMI is 0.6 more from the sum of the relative excess risks due to a)
    low MDScale but low BMI and b) the presence of high BMI but high
    MDScale (individual effects). So, that is an indication for
    super-additive interaction, but the difference is not significant as
    RERI’s CI contains 0.

  - The combined presence of low MDScale and high BMI results in a 29%
    increased excess relative risk, as compared to the excess relative
    risk due to low MDScale but low BMI and the presence of high BMI but
    high MDScale (this is also not significant, as SI’s CI contains 1
    and it’s something we expect as *SI* and *RERI* are equivalent).

  - Obesity and low dedication to Mediterranean diet have
    sub-multiplicative interaction effect on hypertension (as 0.94 \<
    1). However, the 95% CI is (0.13 , 6.85), which includes 1. Thus, we
    cannot reject the null hypothesis of no multiplicative interaction.
    In this case, we conclude that the effects do not deviate from
    multiplicativity.

## Interactions - Approach 2

Here, we are going to examine the same indices. The difference is that
we are going to follow the classical approach, i.e. fitting a logistic
regression model that contains the main effects of obesity and
Mediterranean diet, as well as an interaction term for them.

``` r
fit2 <- glm(hyper ~ bmi*nomd, family = binomial, data)

r_table(fit2) %>%
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

bmi\>30

</td>

<td style="text-align:right;">

2.611

</td>

<td style="text-align:right;">

1.982

</td>

<td style="text-align:right;">

1.264

</td>

<td style="text-align:left;">

0.206

</td>

<td style="text-align:right;">

0.589

</td>

<td style="text-align:right;">

11.565

</td>

</tr>

<tr>

<td style="text-align:left;">

nomd\<5

</td>

<td style="text-align:right;">

1.547

</td>

<td style="text-align:right;">

1.463

</td>

<td style="text-align:right;">

0.461

</td>

<td style="text-align:left;">

0.645

</td>

<td style="text-align:right;">

0.242

</td>

<td style="text-align:right;">

9.878

</td>

</tr>

<tr>

<td style="text-align:left;">

bmi\>30:nomd\<5

</td>

<td style="text-align:right;">

0.938

</td>

<td style="text-align:right;">

0.951

</td>

<td style="text-align:right;">

\-0.063

</td>

<td style="text-align:left;">

0.949

</td>

<td style="text-align:right;">

0.128

</td>

<td style="text-align:right;">

6.847

</td>

</tr>

</tbody>

</table>

In this case, instead of using `nlConfint()` as before, we can use
`interactionR()` from **interactionR** package which automatically
calculates the quantities we want.

``` r
int <- interactionR(fit2, exposure_names = c("bmi", "nomd"), em = FALSE)

int$dframe %>%
  filter(Measures %in% c("RERI", "SI", "Multiplicative scale")) %>%
  kbl(digits = 4) %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Measures

</th>

<th style="text-align:right;">

Estimates

</th>

<th style="text-align:right;">

CI.ll

</th>

<th style="text-align:right;">

CI.ul

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Multiplicative scale

</td>

<td style="text-align:right;">

0.9377

</td>

<td style="text-align:right;">

0.1284

</td>

<td style="text-align:right;">

6.8474

</td>

</tr>

<tr>

<td style="text-align:left;">

RERI

</td>

<td style="text-align:right;">

0.6293

</td>

<td style="text-align:right;">

\-2.4824

</td>

<td style="text-align:right;">

3.7409

</td>

</tr>

<tr>

<td style="text-align:left;">

SI

</td>

<td style="text-align:right;">

1.2916

</td>

<td style="text-align:right;">

0.2872

</td>

<td style="text-align:right;">

5.8096

</td>

</tr>

</tbody>

</table>

The results are identical with the previous approach. However, in both
approaches, we have to make sure that the variables that are included in
interaction terms are coded in such way that they can can be perceived
as risk factors.
