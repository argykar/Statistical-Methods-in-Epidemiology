Statistical Methods in Epidemiology
================
14 October 2020

Παρακαλόυμε κάθε lab να διαβάζεται σε συνδυασμό με το βιβλίο [R 4 Data
Science](https://r4ds.had.co.nz/) (Hadley Wickham 2017)

-----

Το παρακάτω θα εγκαταστήσει και θα φορτώσει αυτόματα τις βιβλιοθήκες
στην R σε περίπτωση που δεν τις έχετε. Αν τις έχετε, αγνοήστε το και
φορτώστε κανονικά με το `library()`.

``` r
   packages <- c('tidyverse', 'haven', 'lubridate', 'plotly', 'knitr',
                 'epiR')
   if (require(packages) == FALSE) {
     install.packages(packages, repos = "http://cran.us.r-project.org")
   }

   lapply(packages, require, character.only = TRUE)
```

-----

Θέτουμε αρχικά το working directory μας.

``` r
#setwd('path/to/working_directory')
```

Το αρχείο είναι .dta (stata). Μέσω της βιβλιοθήκης `haven` η οποία έχει
την εντολή `read_dta()` φορτώνουμε το dataset.

``` r
diet <- read_dta('diet.dta')
```

## Ερώτημα 1

Η εντολή `str()` λειτουργεί ακριβώς όπως η **describe** του STATA.

``` r
str(diet)
```

    ## tibble [337 x 13] (S3: tbl_df/tbl/data.frame)
    ##  $ id    : num [1:337] 1 2 3 4 5 6 7 8 9 10 ...
    ##   ..- attr(*, "label")= chr "Subject identity number"
    ##   ..- attr(*, "format.stata")= chr "%9.0g"
    ##  $ doe   : Date[1:337], format: "1964-08-16" "1964-12-16" ...
    ##  $ dox   : Date[1:337], format: "1976-12-01" "1976-12-01" ...
    ##  $ chd   : num [1:337] 0 0 0 0 0 0 0 0 0 0 ...
    ##   ..- attr(*, "label")= chr "Outcome: 1= chd, 0 otherwise"
    ##   ..- attr(*, "format.stata")= chr "%9.0g"
    ##  $ dob   : Date[1:337], format: "1915-01-04" "1914-06-03" ...
    ##  $ job   : num [1:337] 0 0 0 0 0 0 0 0 0 0 ...
    ##   ..- attr(*, "label")= chr "Occupation"
    ##   ..- attr(*, "format.stata")= chr "%8.0g"
    ##  $ month : num [1:337] 8 12 11 9 9 3 11 5 2 7 ...
    ##   ..- attr(*, "label")= chr "month of survey"
    ##   ..- attr(*, "format.stata")= chr "%8.0g"
    ##  $ energy: num [1:337] 2.87 1.98 2.67 2.84 2.94 ...
    ##   ..- attr(*, "label")= chr "Total energy (1000kcals/day)"
    ##   ..- attr(*, "format.stata")= chr "%9.0g"
    ##  $ height: num [1:337] 175 164 169 167 174 ...
    ##   ..- attr(*, "label")= chr "Height (cm)"
    ##   ..- attr(*, "format.stata")= chr "%9.0g"
    ##  $ weight: num [1:337] 71.5 70.1 71.9 74.9 78.4 ...
    ##   ..- attr(*, "label")= chr "Weight (kg)"
    ##   ..- attr(*, "format.stata")= chr "%9.0g"
    ##  $ fat   : num [1:337] 141.7 85.8 107.7 132.2 126.3 ...
    ##   ..- attr(*, "label")= chr "Total fat (g/day)"
    ##   ..- attr(*, "format.stata")= chr "%9.0g"
    ##  $ fibre : num [1:337] 17.83 9.49 15.99 17.04 14.54 ...
    ##   ..- attr(*, "label")= chr "Total fibre (g/day)"
    ##   ..- attr(*, "format.stata")= chr "%9.0g"
    ##  $ hieng : num [1:337] 1 0 0 1 1 0 0 1 0 1 ...
    ##   ..- attr(*, "label")= chr "Indicator for energy > 2.75"
    ##   ..- attr(*, "format.stata")= chr "%9.0g"
    ##  - attr(*, "label")= chr "Diet data with dates"

Η μεταβλητή *chd* του dataset είναι ακέραιος (0, 1), επομένως θα την
μετατρέψουμε σε κατηγορική. Με την `factor()` παίρνουμε το
ζητούμενο. Στη συνέχεια, θα ορίσουμε τις κατηγορίες της, με την
`levels()`.

``` r
diet$chd <- factor(diet$chd, levels = c('1', '0'), labels = c('chd', 'otherwise'))
diet$hieng <- factor(diet$hieng, levels = c('1', '0'), labels = c('high', 'low'))
```

Θα παρατηρήσετε πως συχνά γίνεται χρήση του pipe (`%>%`). Το pipe
βρίσκεται στην βιβλιοθήκη **magrittr** και χρησιμοποιήτε για
λόγους ευκολίας, καθώς η R χρησιμοποιείται από κάθε είδους
επιστήμονα και όχι μόνο από προγραμματιστές. Το pipe
περιλαμβάνεται στο πακέτο **tidyverse** και σε όλα τα συναφή.

Στο STATA, με την **tab chd** λαμβάνουμε περιγραφικά για την κατηγορική
*chd*. Παρακάτω, με την χρήση του pipe, παίρνουμε ακριβώς τον ίδιο
πίνακα.

``` r
diet %>%
  group_by(chd) %>%
  summarise(count = n()) %>%
  mutate(percent = count * 100 / sum(count),
         cum = cumsum(percent))
```

    ## # A tibble: 2 x 4
    ##   chd       count percent   cum
    ##   <fct>     <int>   <dbl> <dbl>
    ## 1 chd          46    13.6  13.6
    ## 2 otherwise   291    86.4 100

Για να γίνει πιο κατανοητή η χρήση του pipe:

  - από το dataset μας (*diet*), γκρούπαρε με βάση το *chd*
    (`group_by()`)
  - συνοψίζουμε την πληροφορία (`summarise()`) και στη συνέχεια ζητάμε
    τις συχνότητες της κατηγορικής (`n()`), ενώ **count** είναι το
    όνομα της μεταβλητής που θα εμφανιστεί στο output.
  - στη συνέχεια, με την `mutate()` δημιουργούμε δύο νέες μεταβλητές
    (ποσοστό κατηγοριών από το σύνολο, αθροιστική συχνότητα,
    `cumsum()`)

## Ερώτημα 2

Σε περίπτωση που θέλουμε να δούμε τα δεδομένα μας, είτε για λόγους
δημιουργίας νέων μεταβλητών, είτε για λόγους τροποποίησης των ήδη
υπάρχοντων, μπορεί να γίνει με την `View()` για όλο το dataset ή την
`head()` για 5 γραμμές ή περισσότερες.

Στο εργαστήριο χρησιμοποιείται το **list** του STATA για τις 20 πρώτες
γραμμές.

``` r
head(diet, n = 20)
```

    ## # A tibble: 20 x 13
    ##       id doe        dox        chd   dob          job month energy height weight
    ##    <dbl> <date>     <date>     <fct> <date>     <dbl> <dbl>  <dbl>  <dbl>  <dbl>
    ##  1     1 1964-08-16 1976-12-01 othe~ 1915-01-04     0     8   2.87   175.   71.5
    ##  2     2 1964-12-16 1976-12-01 othe~ 1914-06-03     0    12   1.98   164.   70.1
    ##  3     3 1965-11-16 1976-12-01 othe~ 1907-02-03     0    11   2.67   169.   71.9
    ##  4     4 1965-09-16 1976-12-01 othe~ 1906-12-25     0     9   2.84   167.   74.9
    ##  5     5 1965-09-16 1976-03-31 othe~ 1906-04-01     0     9   2.94   174.   78.4
    ##  6     6 1965-03-16 1968-08-31 othe~ 1914-03-23     0     3   2.47   177.   72.4
    ##  7     7 1958-11-16 1976-12-01 othe~ 1913-09-26     0    11   2.56   169.   64.2
    ##  8     8 1965-05-16 1976-12-01 othe~ 1914-12-11     0     5   2.99   166.   73.8
    ##  9     9 1959-02-16 1962-01-10 othe~ 1892-01-10     0     2   2.31   166.   49.1
    ## 10    10 1964-07-16 1974-05-16 othe~ 1904-05-16     0     7   3.12   181.   78.3
    ## 11    11 1964-10-16 1974-04-08 othe~ 1904-04-08     0    10   2.16   174.   63.9
    ## 12    12 1964-07-16 1974-08-03 othe~ 1904-08-03     0     7   3.77   175.   77.5
    ## 13    13 1964-09-16 1974-02-16 othe~ 1904-02-17     0     9   2.87   165.   80.2
    ## 14    14 1959-12-16 1976-12-01 othe~ 1916-01-03     0    12   2.34   172.   89.1
    ## 15    15 1962-05-16 1976-08-20 othe~ 1906-08-21     0     5   2.81   164.   59.6
    ## 16    16 1959-05-16 1959-12-31 chd   1896-09-15     0     5   2.22   171.   89.4
    ## 17    17 1959-02-16 1965-01-14 othe~ 1899-05-06     0     2   3.18   177.   85.7
    ## 18    18 1959-02-16 1968-03-08 othe~ 1898-03-08     0     2   3.20   178.   94.8
    ## 19    19 1959-02-16 1966-03-12 othe~ 1898-03-01     0     2   2.77   169.   83.7
    ## 20    20 1959-02-16 1969-12-26 othe~ 1899-12-26     0     2   3.05   184.   85.5
    ## # ... with 3 more variables: fat <dbl>, fibre <dbl>, hieng <fct>

Αν θέλουμε συγκεκριμένες μεταβλητές του data.frame μπορούμε να τις
προσδιορίσουμε μέσω της `select()`

``` r
diet %>%
  select(id, doe, dox, chd) %>%
  head(20)
```

    ## # A tibble: 20 x 4
    ##       id doe        dox        chd      
    ##    <dbl> <date>     <date>     <fct>    
    ##  1     1 1964-08-16 1976-12-01 otherwise
    ##  2     2 1964-12-16 1976-12-01 otherwise
    ##  3     3 1965-11-16 1976-12-01 otherwise
    ##  4     4 1965-09-16 1976-12-01 otherwise
    ##  5     5 1965-09-16 1976-03-31 otherwise
    ##  6     6 1965-03-16 1968-08-31 otherwise
    ##  7     7 1958-11-16 1976-12-01 otherwise
    ##  8     8 1965-05-16 1976-12-01 otherwise
    ##  9     9 1959-02-16 1962-01-10 otherwise
    ## 10    10 1964-07-16 1974-05-16 otherwise
    ## 11    11 1964-10-16 1974-04-08 otherwise
    ## 12    12 1964-07-16 1974-08-03 otherwise
    ## 13    13 1964-09-16 1974-02-16 otherwise
    ## 14    14 1959-12-16 1976-12-01 otherwise
    ## 15    15 1962-05-16 1976-08-20 otherwise
    ## 16    16 1959-05-16 1959-12-31 chd      
    ## 17    17 1959-02-16 1965-01-14 otherwise
    ## 18    18 1959-02-16 1968-03-08 otherwise
    ## 19    19 1959-02-16 1966-03-12 otherwise
    ## 20    20 1959-02-16 1969-12-26 otherwise

## Ερώτημα 3

Δημιουργούμε την μεταβλητή time-since-entry (**tse**) \[in years
(365.25)\] και την **failure**, η οποία δεν έχει καμία διαφορά με την
**chd**. Γίνεται χρήση της `ifelse()` η οποία έχει ως πρώτο όρισμα ένα
λογικό επιχείρημα, ως δεύτερο μία τιμή σε περίπτωση που το λογικό
επιχείρημα είναι *YES* και ως τρίτη μία τιμή σε περίπτωση που είναι
*NO*.

Πιο συγκεκριμένα αν **chd** είναι 1, βάλε στην **failure** 1,
διαφορετικά 0.

Στη συνέχεια, αθροίζουμε τις αποτυχίες με την `summarise(sum())`

``` r
diet %>%
  mutate(tse = (dox - doe) / 365.25,
         failure = ifelse(chd == 'otherwise', 0, 1)) %>% 
  summarise(failures = sum(failure))
```

    ## # A tibble: 1 x 1
    ##   failures
    ##      <dbl>
    ## 1       46

``` r
dim(diet) # μας δίνει γραμμές και στήλες του data.frame
```

    ## [1] 337  13

Σε αντίθεση με το STATA, στην R δεν χρείαζεται να προσδιορίσουμε από
πριν τι είδους δεδομένα έχουμε. Οτιδήποτε χρειαστεί το δηλώνουμε
μέσα στην συνάρτηση που θα χρησιμοποιήσουμε για το μοντέλο μας ή
ό,τι είναι αυτό που θέλουμε να τρέξουμε, αρκεί βέβαια να έχουμε
καθορίσει τον σωστό τύπο μεταβλητών (integer, factor, …)

*Total time at risk*

``` r
diet %>%
  mutate(tse = (dox - doe) / 365.25) %>% 
  summarise(time_at_risk = sum(tse))
```

    ## # A tibble: 1 x 1
    ##   time_at_risk 
    ##   <drtn>       
    ## 1 4603.669 days

ΣΗΜΕΙΩΣΗ: δεν χρειάζεται κάθε φορά να δημιουργούμε την νέα μεταβλητή,
μπορούμε να την κάνουμε assign π.χ.

``` r
diet_tse <- diet %>%
  mutate(tse = lubridate::time_length(difftime(diet$dox, diet$doe), 'years'))

#Time difference (entry - exit) in years
diet_tse$tse
```

    ##   [1] 12.2929500 11.9589322 11.0417522 11.2087611 10.5379877  3.4606434
    ##   [7] 18.0424367 11.5455168  2.8993840  9.8316222  9.4757016 10.0479124
    ##  [13]  9.4182067 16.9609856 14.2642026  0.6269678  5.9110198  9.0568104
    ##  [19]  7.0663929 10.8583162  9.0759754  4.7665982  4.9609856  4.4709103
    ##  [25]  0.2874743 17.3798768 17.3798768  9.3032170 15.4305270 17.3798768
    ##  [31] 17.2950034 11.5099247  9.9000684  7.7097878  0.6214921  5.5770021
    ##  [37] 17.4620123 17.4620123  4.4572211  3.5044490 12.2792608 12.6324435
    ##  [43]  9.4647502 12.2464066 12.3449692 11.1895962 17.7138946 17.7905544
    ##  [49] 12.7802875 17.4620123 16.9582478 10.7022587 12.4818617 13.9055441
    ##  [55]  8.8186174  7.1567420 16.5448323 11.7125257  8.9691992  5.3333333
    ##  [61] 17.4620123  4.3723477 17.1279945 17.1279945 16.7118412 16.8761123
    ##  [67] 17.1279945 17.1279945 17.1279945 17.0431211 17.0431211 16.0410678
    ##  [73] 16.1259411 16.2929500 16.1259411 16.0410678  1.5824778  5.6180698
    ##  [79]  9.0896646 16.1259411 15.9589322 15.9589322 16.1259411 15.9589322
    ##  [85] 15.8740589 15.9589322 15.7125257 11.2936345  9.9548255 11.1594798
    ##  [91] 11.9561944 15.7125257 14.8884326 15.7891855 15.7125257 11.8740589
    ##  [97] 11.8740589 11.7125257 12.0410678 10.7132101 10.8747433 10.8747433
    ## [103]  8.6789870 16.8761123  9.2073922 14.4613279 14.0424367 12.8761123
    ## [109]  9.8726899 16.6269678 11.6030116 16.9609856 16.8761123  9.2895277
    ## [115] 11.7809719 14.7679671 13.6317591  4.1943874 10.2121834 11.9890486
    ## [121] 13.2484600 11.6358658  9.2019165  2.8446270 13.9739904 14.0095825
    ## [127] 16.7912389 16.7912389 16.7912389 11.1266256  4.6187543 11.5455168
    ## [133] 11.4606434 11.2826831  6.2915811 11.5865845 10.0041068 16.7118412
    ## [139] 10.7132101  6.7268994 10.6283368 11.1266256  6.2368241 10.0588638
    ## [145] 12.2491444 17.2950034 10.4722793 17.2950034 16.3778234 12.2628337
    ## [151] 17.2101300 17.0431211 17.0431211 17.0431211 17.0431211 16.4599589
    ## [157] 16.5448323  0.5612594 16.5448323  9.2046543  9.5605749 15.8740589
    ## [163] 13.3689254  6.7843943  4.7036277 16.4599589 16.4599589 16.2080767
    ## [169] 16.2080767 15.9589322  8.6844627 16.0410678 15.9589322 14.6338125
    ## [175] 10.3791923 15.7891855 10.6283368 10.5462012 15.7125257  8.4407940
    ## [181] 11.0609172  5.3634497 14.6584531 12.0410678 11.7125257 10.2094456
    ## [187] 20.0410678 20.0410678 19.4524298 20.0410678 18.3189596 19.3374401
    ## [193] 15.3401780  9.4537988 20.0410678  2.4339493  3.1101985 19.9589322
    ## [199] 10.8473648 19.9589322  1.7002053 19.7891855 19.9589322 19.8740589
    ## [205] 16.1478439 19.8740589 19.8740589 17.8288843 19.8740589 18.6858316
    ## [211] 16.5092402 17.3004791 19.7891855 19.0417522 19.0417522 17.5934292
    ## [217] 19.3785079 19.3785079 19.8740589 19.8740589 16.3285421 15.3949350
    ## [223] 19.7891855  4.7994524 19.7125257 19.6112252 19.6276523 18.1136208
    ## [229] 12.3860370 19.5455168 19.5455168 19.5455168 15.9397673 14.6365503
    ## [235]  7.1238877 12.1232033 19.2087611 10.1218344 19.5455168  7.7453799
    ## [241] 19.3785079 10.0369610 14.5708419 16.1013005 18.9596167 17.8316222
    ## [247] 11.2744695  3.6194387  7.8740589 18.9596167 18.7898700 11.4798084
    ## [253] 18.8747433 18.9596167 18.9596167 18.9596167 18.7898700 17.8206708
    ## [259] 16.5886379 18.7132101  3.4524298 13.9822040 17.7905544 17.2320329
    ## [265] 16.4982888  3.0855578 14.9267625 17.7905544 17.7138946 17.7138946
    ## [271] 17.7138946  4.7446954 14.7843943 17.7138946 17.7138946 17.7138946
    ## [277] 17.7138946  4.3531828 17.6290212 15.2060233 17.2073922 10.9705681
    ## [283] 16.7118412 16.7118412 16.7118412 16.1259411 16.1259411 16.1259411
    ## [289]  0.9007529 16.0410678 16.0410678 16.0410678 16.9609856 16.9609856
    ## [295] 16.9609856 16.9609856 16.8761123 16.9609856 16.9609856 15.4195756
    ## [301] 16.8761123 16.8761123 16.8761123 15.4606434  1.4948665 16.8761123
    ## [307] 16.7912389 16.1533196 10.7761807 16.5831622  3.5181383 16.1259411
    ## [313] 15.5455168 15.5455168 15.5455168 15.5455168 15.1238877 15.1266256
    ## [319] 15.5455168 15.5455168 15.5455168 15.4606434 15.4606434  4.4763860
    ## [325] 15.4606434 15.4606434 15.4606434 15.4606434 13.6536619  9.0376454
    ## [331] 15.0417522 15.4606434 15.4606434 15.3785079 15.3785079 15.2936345
    ## [337] 15.2936345

Όσο περισσότερο γράφετε τόσο συνηθίζετε στο σκεπτικό του pipe.
Εναλλατικά, μπορεί κανείς να χρησιμοποιήσει base R.

## Ερώτημα 4

``` r
diet_tse %>%
  mutate(failure = ifelse(chd == 'otherwise', 0, 1),
         '_t0' = 0,
         '_t' = tse,
         '_d' = failure,
         '_st' = 1) %>%
  select(id, '_t0', '_t', '_d', '_st') %>%
  head(20)
```

    ## # A tibble: 20 x 5
    ##       id `_t0`   `_t`  `_d` `_st`
    ##    <dbl> <dbl>  <dbl> <dbl> <dbl>
    ##  1     1     0 12.3       0     1
    ##  2     2     0 12.0       0     1
    ##  3     3     0 11.0       0     1
    ##  4     4     0 11.2       0     1
    ##  5     5     0 10.5       0     1
    ##  6     6     0  3.46      0     1
    ##  7     7     0 18.0       0     1
    ##  8     8     0 11.5       0     1
    ##  9     9     0  2.90      0     1
    ## 10    10     0  9.83      0     1
    ## 11    11     0  9.48      0     1
    ## 12    12     0 10.0       0     1
    ## 13    13     0  9.42      0     1
    ## 14    14     0 17.0       0     1
    ## 15    15     0 14.3       0     1
    ## 16    16     0  0.627     1     1
    ## 17    17     0  5.91      0     1
    ## 18    18     0  9.06      0     1
    ## 19    19     0  7.07      0     1
    ## 20    20     0 10.9       0     1

## Ερώτημα 6

``` r
diet %>%
  group_by(hieng) %>%
  summarise(count = n()) %>%
  mutate(percent = count * 100 / sum(count),
         cum = cumsum(percent))
```

    ## # A tibble: 2 x 4
    ##   hieng count percent   cum
    ##   <fct> <int>   <dbl> <dbl>
    ## 1 high    182    54.0  54.0
    ## 2 low     155    46.0 100

## Ερώτημα 7

Με την `filter()` επιλέγουμε μόνο τους cases, αφήνοντας απέξω τους chd
== ‘otherwise’

``` r
diet_tse %>%
  filter(chd == 'chd') %>%
  group_by(hieng) %>%
  summarise(D = n(), Y = sum(diet_tse$tse) / 1000, rate = D / Y)
```

    ## # A tibble: 2 x 4
    ##   hieng     D     Y  rate
    ##   <fct> <int> <dbl> <dbl>
    ## 1 high     18  4.60  3.91
    ## 2 low      28  4.60  6.08

## Ερώτημα 8

Δημιουργούμε την μεταβλητή **htgrp**.

Ο συλλογισμός για το παρακάτω:

  - `filter()` : αφήνουμε απέξω τις ΝΑ τιμές
  - `mutate()` : δημιουργούμε την νέα μεταβλητή μας
      - `case_when()` : μετατρέπει την συνεχή σε κατηγορική μεταβλητή
      - `between()` : μας επιστρέφει TRUE/FALSE, αν η τιμή που είναι στο
        πρώτο όρισμα βρίσκεται μέσα στο διάστημα (δεύτερο, τρίτο όρισμα)
  - `group_by()` : ομαδοποίηση με βάση την κατηγορική που φτιάξαμε
  - `summarise()` : συνοψίζουμε με βάση τις συχνότητες της κάθε
    κατηγορίας
  - `mutate()` : (προαιρετικό) φτιάχνουμε και τις υπόλοιπες στήλες (το
    δίνει το output του STATA, γι’αυτό και το κάνουμε)

<!-- end list -->

``` r
diet_htgrp <- diet_tse %>%
  filter(!is.na(height)) %>%
  mutate(htgrp = factor(case_when(
    between(height, min(diet_tse$height, na.rm = TRUE), 169.999) ~ 0,
    between(height, 170,174.999) ~ 1,
    between(height, 175, 179.999) ~ 2,
    between(height, 180, 195) ~ 3)
    ))

diet_htgrp
```

    ## # A tibble: 332 x 15
    ##       id doe        dox        chd   dob          job month energy height weight
    ##    <dbl> <date>     <date>     <fct> <date>     <dbl> <dbl>  <dbl>  <dbl>  <dbl>
    ##  1     1 1964-08-16 1976-12-01 othe~ 1915-01-04     0     8   2.87   175.   71.5
    ##  2     2 1964-12-16 1976-12-01 othe~ 1914-06-03     0    12   1.98   164.   70.1
    ##  3     3 1965-11-16 1976-12-01 othe~ 1907-02-03     0    11   2.67   169.   71.9
    ##  4     4 1965-09-16 1976-12-01 othe~ 1906-12-25     0     9   2.84   167.   74.9
    ##  5     5 1965-09-16 1976-03-31 othe~ 1906-04-01     0     9   2.94   174.   78.4
    ##  6     6 1965-03-16 1968-08-31 othe~ 1914-03-23     0     3   2.47   177.   72.4
    ##  7     7 1958-11-16 1976-12-01 othe~ 1913-09-26     0    11   2.56   169.   64.2
    ##  8     8 1965-05-16 1976-12-01 othe~ 1914-12-11     0     5   2.99   166.   73.8
    ##  9     9 1959-02-16 1962-01-10 othe~ 1892-01-10     0     2   2.31   166.   49.1
    ## 10    10 1964-07-16 1974-05-16 othe~ 1904-05-16     0     7   3.12   181.   78.3
    ## # ... with 322 more rows, and 5 more variables: fat <dbl>, fibre <dbl>,
    ## #   hieng <fct>, tse <dbl>, htgrp <fct>

``` r
diet_htgrp %>%
  group_by(htgrp) %>%
  summarise(count = n()) %>%
  mutate(percent = count * 100 / sum(count),
         cum = cumsum(percent))
```

    ## # A tibble: 4 x 4
    ##   htgrp count percent   cum
    ##   <fct> <int>   <dbl> <dbl>
    ## 1 0        92    27.7  27.7
    ## 2 1       102    30.7  58.4
    ## 3 2        83    25    83.4
    ## 4 3        55    16.6 100

## ΙΙ. Rate Ratios.

Αρχικά, κωδικοποιούμε κατάλληλα τις μεταβλητες **hieng** και **chd**.
Στη συνέχεια, φτιάχνουμε τον πίνακα με την μεταβλητή
**analysis\_time** που είναι τα person years των cases σε κάθε κατηγορία
της **hieng**. Στο τέλος, με την `mutate()` δημιουργούμε τα incidence
rates κάθε κατηγορίας, μαζί με τα 95% διαστήματα εμπιστοσύνης τους.

Όσον αφορά τo πακέτο **lubridate**, χρησιμοποιήθηκε η συνάρτηση
`time_length()` με την οποία μπορεί κανείς να κάνει υπολογίσει χρόνο σε
χρόνια, μήνες, βδομάδες, μέρες με το argument *(… ,unit = ’’)*.
Χρησιμοποιήθηκε επίσης, η `difftime()` με την οποία υπολογίζουμε
διαφορές μεταξύ ημερομηνιών.

``` r
diet_tse %>%
  mutate(analysis_time = lubridate::time_length(difftime(dox, doe), unit = 'years')) %>%
  group_by(hieng) %>%
  summarise(total = sum(chd == 'chd'), pyrs = sum(analysis_time)) %>%
  mutate(Inc_Rate = 1000 * total/pyrs,
         lower_l = Inc_Rate / exp(1.96 * sqrt(1/total)),
         upper_l = Inc_Rate * exp(1.96 * sqrt(1/ total))) -> by_hieng

by_hieng
```

    ## # A tibble: 2 x 6
    ##   hieng total  pyrs Inc_Rate lower_l upper_l
    ##   <fct> <int> <dbl>    <dbl>   <dbl>   <dbl>
    ## 1 high     18 2544.     7.07    4.46    11.2
    ## 2 low      28 2059.    13.6     9.39    19.7

Για τα rate ratio και το
<img src="https://render.githubusercontent.com/render/math?math=X^2">
test της εκτίμησης θα χρησιμοποιήσουμε το πακέτο **epiR**, που
περιλαμβάνει την συνάρτηση `epi.2by2()`. Η `epi.2by2()` δέχετε
**table** ως argument γι’ αυτό και κάνουμε την μετατροπη με την
`as.table()`. Από το data.frame **by\_hieng** θέλουμε την δεύτερη και
τρίτη στήλη (cases, pyrs), το οποίο λαμβάνουμε με ένα απλό slicing
του data.frame. Δηλώνοντας ως argument *(…, method = ‘cohort.time’)*,
λαμβάνουμε υπόψιν τα person years των cases.

Τέλος, το output της `epi.2by2()` δίνει πολλά αποτελέσματα που δεν
χρειαζόμαστε. Χρησιμοποιώντας τον accessor *$* μπορούμε να
εξάγουμε αυτά που μας ενδιαφέρουν και να δημιουργήσουμε μία
λίστα για να τα συνοψίσουμε σε ένα object.

``` r
RR <- as.table(as.matrix(by_hieng[, 2:3])) %>%
  epi.2by2('cohort.time')

RR
```

    ##              Outcome +    Time at risk        Inc rate *
    ## Exposed +           18            2544             0.707
    ## Exposed -           28            2059             1.360
    ## Total               46            4604             0.999
    ## 
    ## Point estimates and 95% CIs:
    ## -------------------------------------------------------------------
    ## Inc rate ratio                               0.52 (0.27, 0.97)
    ## Attrib rate *                                -0.65 (-1.25, -0.05)
    ## Attrib rate in population *                  -0.36 (-0.94, 0.22)
    ## Attrib fraction in exposed (%)               -92.17 (-268.89, -2.62)
    ## Attrib fraction in population (%)            -36.07 (-47.43, -23.50)
    ## -------------------------------------------------------------------
    ##  Wald confidence limits
    ##  CI: confidence interval
    ##  * Outcomes per 100 units of population time at risk

``` r
res_hieng <- list(table = RR$tab, RR = RR$res$RR.crude.wald, ChiSquare_Test = RR$res$chisq.crude)

res_hieng
```

    ## $table
    ##              Outcome +    Time at risk        Inc rate *
    ## Exposed +           18            2544             0.707
    ## Exposed -           28            2059             1.360
    ## Total               46            4604             0.999
    ## 
    ## $RR
    ##         est    lower    upper
    ## 1 0.5237295 0.290522 0.944137
    ## 
    ## $ChiSquare_Test
    ##   test.statistic df    p.value
    ## 1        4.79282  1 0.02857861

## III. Exposure with more than two levels

Κατά τα γνωστά, δημιουργούμε την μεταβλητη **energy\_cat** και
φτιάχνουμε το παρακάτω πινακάκι για τα περιγραφικά της.

  - Πολυ χρήσιμη η `cut()` με την οποία μπορούμε να δηλώσουμε breaks για
    την νέα μεταβλητή, καθώς και τα αντίστοιχα labels της.

<!-- end list -->

``` r
energycat_sum <- diet %>%
  mutate(energy_cat = cut(energy, breaks=c(1.5, 2.5, 3.0, 4.5), 
                        labels=c("low","middle","high"))) %>%
  group_by(energy_cat) %>%
  summarise(freq = n()) %>%
  mutate(percent = freq * 100 / sum(freq),
         cum = cumsum(percent))

energycat_sum
```

    ## # A tibble: 3 x 4
    ##   energy_cat  freq percent   cum
    ##   <fct>      <int>   <dbl> <dbl>
    ## 1 low           75    22.3  22.3
    ## 2 middle       150    44.5  66.8
    ## 3 high         112    33.2 100

Για την εντολή **strate eng3, per(1000)** του STATA, η οποία μας δίνει
τα rates για κάθε κατηγορία της μεταβλητής μαζί με τα 95% CI τους
κάνουμε τα εξής:

``` r
by_engcat <- diet_tse %>%
  mutate(energy_cat = cut(energy, breaks=c(1.5, 2.5, 3.0, 4.5), 
                        labels=c("low","middle","high"))) %>%
  group_by(energy_cat) %>%
  summarise(cases = sum(chd == 'chd'), 
            pyrs = sum(tse),
            rate = cases / pyrs,
            lower = rate / exp(1.96 * sqrt(1 / cases)),
            upper = rate * exp(1.96 * sqrt(1 / cases)))

by_engcat
```

    ## # A tibble: 3 x 6
    ##   energy_cat cases  pyrs    rate   lower   upper
    ##   <fct>      <int> <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 low           16  947. 0.0169  0.0104  0.0276 
    ## 2 middle        22 2017. 0.0109  0.00718 0.0166 
    ## 3 high           8 1640. 0.00488 0.00244 0.00976

Για την δημιουργία του γραφήματος θα χρησιμοποιήσουμε το πακέτο
**ggplot2**. Παρακάτω παρουσιάζεται η λογική για το χτίσιμο του
γραφήματος.

``` r
theme_set(theme_light()) # Θέτουμε θέμα για τα γραφήματα μας, σε περίπτωση που δεν αρέσει σε κάποιον το basic

energycat_plot <- by_engcat %>%
  ggplot(aes(energy_cat, rate)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                size = 1, width = 0.2, color = 'red') +
  labs(x = 'Energy categories (eng3)',
       y = 'Rate (per 1000)',
       title = ' Rates of CHD among three levels of exposure')

energycat_plot
```

![](lab2_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

Για το ylog plot αρκεί στο `aes()` να βάλουμε αντί για *rate*,
*log(rate)*

``` r
by_engcat %>%
  ggplot(aes(energy_cat, log(rate))) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = log(lower), ymax = log(upper)),
                size = 1, width = 0.2, color = 'red') +
  labs(x = 'Energy categories (eng3)',
       y = 'log(Rate) (per 1000)',
       title = ' Rates of CHD among three levels of exposure')
```

![](lab2_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

Για την εντολή **stmh eng3, c(1,0)** του STATA, η οποία μας δίνει το RR
για την πρώτη και την δεύτερη κατηγορία, μαζί με το 95% CI, κάνουμε τα
εξής:

``` r
RR_engcat <- as.table(as.matrix(by_engcat[1:2, 2:3])) %>%
  epi.2by2('cohort.time')

RR_engcat
```

    ##              Outcome +    Time at risk        Inc rate *
    ## Exposed +           16             947              1.69
    ## Exposed -           22            2017              1.09
    ## Total               38            2964              1.28
    ## 
    ## Point estimates and 95% CIs:
    ## -------------------------------------------------------------------
    ## Inc rate ratio                               1.55 (0.76, 3.09)
    ## Attrib rate *                                0.60 (-0.35, 1.54)
    ## Attrib rate in population *                  0.19 (-0.42, 0.80)
    ## Attrib fraction in exposed (%)               35.48 (-31.44, 67.63)
    ## Attrib fraction in population (%)            14.94 (6.17, 24.67)
    ## -------------------------------------------------------------------
    ##  Wald confidence limits
    ##  CI: confidence interval
    ##  * Outcomes per 100 units of population time at risk

Συνοψίζοντας τα αποτελέσματα της `epi.2by2()` σε ένα object, έχουμε:

``` r
res_engcat <- list(table = RR_engcat$tab, 
                   RR = RR_engcat$res$RR.crude.wald, 
                   ChiSquare_Test = RR_engcat$res$chisq.crude)

res_engcat
```

    ## $table
    ##              Outcome +    Time at risk        Inc rate *
    ## Exposed +           16             947              1.69
    ## Exposed -           22            2017              1.09
    ## Total               38            2964              1.28
    ## 
    ## $RR
    ##        est     lower    upper
    ## 1 1.540669 0.8128729 2.920088
    ## 
    ## $ChiSquare_Test
    ##   test.statistic df   p.value
    ## 1       1.780102  1 0.1821368

Όμοια, για την εντολή **stmh eng3, c(2,0)**:

``` r
RR_engcat1vs3 <- as.table(as.matrix(by_engcat[c(3,1), 2:3])) %>%
  epi.2by2('cohort.time')

RR_engcat1vs3
```

    ##              Outcome +    Time at risk        Inc rate *
    ## Exposed +            8            1640             0.488
    ## Exposed -           16             947             1.690
    ## Total               24            2586             0.928
    ## 
    ## Point estimates and 95% CIs:
    ## -------------------------------------------------------------------
    ## Inc rate ratio                               0.29 (0.11, 0.71)
    ## Attrib rate *                                -1.20 (-2.10, -0.31)
    ## Attrib rate in population *                  -0.76 (-1.67, 0.15)
    ## Attrib fraction in exposed (%)               -246.44 (-835.02, -39.89)
    ## Attrib fraction in population (%)            -82.15 (-98.80, -62.49)
    ## -------------------------------------------------------------------
    ##  Wald confidence limits
    ##  CI: confidence interval
    ##  * Outcomes per 100 units of population time at risk

Συνοψίζοντας τα αποτελέσματα της `epi.2by2()` σε ένα object, έχουμε:

``` r
res_engcat1vs3 <- list(table = RR_engcat1vs3$tab,
                       RR = RR_engcat1vs3$res$RR.crude.wald, 
                       ChiSquare_Test = RR_engcat1vs3$res$chisq.crude)

res_engcat1vs3
```

    ## $table
    ##              Outcome +    Time at risk        Inc rate *
    ## Exposed +            8            1640             0.488
    ## Exposed -           16             947             1.690
    ## Total               24            2586             0.928
    ## 
    ## $RR
    ##         est     lower     upper
    ## 1 0.2921015 0.1254798 0.6799763
    ## 
    ## $ChiSquare_Test
    ##   test.statistic df     p.value
    ## 1       9.234604  1 0.002374838

## IV. Controlling for Confounding.

Η λογική για τα παρακάτω δεν έχει πολλές διαφορές με τα προηγούμενα. Η
διαφορά εδώ είναι η έξτρα μεταβλητη (*job*) για την οποία θελουμε και
να κάνουμε control. Επομένως, θα γκρουπάρουμε ως προς *hieng* και ως
προς *job*. Τα αποτελέσματα που παίρνουμε μετα το `mutate()` είναι με
βάση τους συνδυασμούς των δύο κατηγοριών. Άρα, ξανα γκρουπάρουμε ως προς
*job* για να γίνει ομαδοποίηση και στη συνέχεια με `summarise()`
παίρνουμε το επιθυμητό αποτέλεσμα.

``` r
diet_tse$job <- factor(diet_tse$job, levels = c('0', '1', '2'))

diet_tse %>%
  group_by(job, hieng) %>%
  summarise(cases = sum(chd == 'chd'),
            pyrs = sum(tse),
            inc_rate = cases * 1000 / pyrs) %>%
  mutate(rr = c(1, inc_rate[1]/inc_rate[2]),
         lower = rr / exp(1.96 * sqrt(1 / cases[1] + 1 / cases[2])),
         upper = rr * exp(1.96 * sqrt(1 / cases[1] + 1 / cases[2]))) %>%
  group_by(job) %>%
  summarise(RR = rr[hieng == 'low'],
            ll = lower[hieng == 'low'],
            ul = upper[hieng == 'low']) -> by_job

by_job
```

    ## # A tibble: 3 x 4
    ##   job      RR    ll    ul
    ##   <fct> <dbl> <dbl> <dbl>
    ## 1 0     0.410 0.124  1.36
    ## 2 1     0.655 0.227  1.89
    ## 3 2     0.518 0.212  1.27

Όμοια με πριν, γκρουπάροντας δύο φορές και φιλτράροντας για ΝΑ και 0
τιμές στα cases, παίρνουμε τον παρακάτω πίνακα για ως προς *job* και
*htgrp*

``` r
by_job_htgrp <- diet_htgrp %>%
  group_by(job, htgrp, hieng) %>%
  summarize(cases = sum(chd == "chd"), py = sum(tse)) %>%
  filter(is.na(htgrp) == FALSE & cases != 0) %>%
  mutate(Rate = cases * 1000 / py) %>%
  group_by(job, htgrp) %>%
  summarise(
    RR = Rate[hieng == 'high'] / Rate[hieng == 'low'],
    ll = (Rate[hieng == 'high'] / Rate[hieng == 'low']) /
      exp(1.96 * sqrt((1 / cases[hieng == 'high']) +
                        (1 / cases[hieng == 'low']))),
    ul = (Rate[hieng == 'high'] / Rate[hieng == 'low']) *
      exp(1.96 * sqrt((1 / cases[hieng == 'high']) +
                        (1 / cases[hieng == 'low'])))
  )

by_job_htgrp
```

    ## # A tibble: 8 x 5
    ## # Groups:   job [3]
    ##     job htgrp    RR     ll    ul
    ##   <dbl> <fct> <dbl>  <dbl> <dbl>
    ## 1     0 0     0.687 0.0967  4.88
    ## 2     0 1     0.411 0.0428  3.95
    ## 3     0 2     0.224 0.0233  2.15
    ## 4     1 0     0.594 0.149   2.37
    ## 5     1 1     1.09  0.182   6.52
    ## 6     2 0     0.414 0.0759  2.26
    ## 7     2 1     0.911 0.184   4.51
    ## 8     2 2     0.496 0.111   2.22

Στην R, δεν υπάρχει έτοιμη εντολή (ή τουλάχιστον δεν την έχουμε βρει
εμείς) για το overall estimate του *hieng* κοντρολάροντας για *job*
και *htgrp*. Ο πιο γρήγορος τρόπος για αυτή την πληροφορία είναι μέσω
μοντέλου, διαφορετικά μπορεί να γίνει χειροκίνητα με τους τύπους.

Θα φτιάξουμε poisson μοντέλο για τα cases με offset τα person-years,
κάνοντας adjust για *job*, *hieng* και *htgrp*

``` r
diet_htgrp %>%
  group_by(job, htgrp, hieng) %>%
  summarise(cases = sum(chd == 'chd'),
            py = sum(tse)) %>%
  filter(!is.na(htgrp) & cases != 0) -> model_data

# Overall estimate controlling for job htgrp
overall_RR <- model_data %>%
  glm(cases ~ job + I(hieng == 'high') + htgrp + offset(log(py)),
      family = 'poisson',
      data = .) %>%
  summary() %>%
  .$coef %>%
  .[3,]

# εκθετικό στο συντελεστή για να πάρουμε το RR
overall_RR[1] <- exp(overall_RR[1])

overall_RR
```

    ##    Estimate  Std. Error     z value    Pr(>|z|) 
    ##  0.56668470  0.30503872 -1.86190208  0.06261689

προφανώς αν θέλουμε το Χ^2, κάνουμε (z-test)^2 και παίρνουμε 3.427124

``` r
# 95% CI for the estimate (RR)
model_data %>%
  glm(cases ~ job + I(hieng == 'high') + htgrp + offset(log(py)),
      family = 'poisson',
      data = .) %>%
  confint() %>%
  exp() %>%
  .[3,]
```

    ##     2.5 %    97.5 % 
    ## 0.3065947 1.0226985

-----

### Extras:

—\> Interactive plot με το πακέτο **plotly** (πολύ χρήσιμο για
παρουσιάσεις markdown)

``` r
p <- ggplotly(energycat_plot)
```

Kατεβάστε τα αρχεία από το GitHub και ανοίξτε το **.html** για να δείτε
το interactive γράφημα.

<!--html_preserve-->

<iframe src="C:/Users/argykar/Statistical Methods in Epidemiology/Lab 2/index.html" width="100%" height="600" scrolling="no" seamless="seamless" frameBorder="0">

</iframe>

<!--/html_preserve-->

-----

## References

<div id="refs" class="references">

<div id="ref-R4DS">

Hadley Wickham, Garret Grolemund. 2017. *R for Data Science*. O’ Reily.
<https://r4ds.had.co.nz/>.

</div>

</div>
