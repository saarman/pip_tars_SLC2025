Weekly **Cx. pipiens** and **Cx. tarsalis** abundance: SLC 2025 field
season
================
Norah Saarman
2026-03-13

- [Load data](#load-data)
- [GLMM by season and urbanization, each mosquito species
  separately](#glmm-by-season-and-urbanization-each-mosquito-species-separately)
  - [Random effects (1 \|
    site_name/disease_week)](#random-effects-1--site_namedisease_week)
  - [Random effects (1 \|
    site_name/collection_date)](#random-effects-1--site_namecollection_date)
  - [Marginal effects](#marginal-effects)

**Research Topic:** testing whether habitat and seasonal partitioning
between Culex pipiens s.l. and Culex tarsalis shapes West Nile Virus
(WNV) dynamics across urban–rural gradients.

**Core hypothesis:** early/mid-season amplification dominated by pipiens
in urban areas, later spillover involving tarsalis moving into
urban/peri-urban areas.

**Approach:** Preliminary results visualized via mapping, with species
identity and abundance as primary response variables. Model mosquito
abundance and proportions using GLMMs count ~ season*urbanization +
trap_type + (1\|site/date), family = poisson(link = “log”):  
- Response variable = mosquito abundance  
- Predictors = season*urbanization  
- The trap type could be important, so we will add that as a fixed
effect (covariate)… is this correct? We do think that the response
variable of count of mosquitoes depends on trap type, since tarsalis
seems to be more attracted to CO2 than pipiens, and we want to quantify
that effect. Note that poisson model does not give a fixed offset (due
to the log link)… The structure of this model means that it will
estimate an effect that scales with the total number of mosquitos
caught, which is exactly what we want.  
- The data are grouped into sites and are also linked through time, so
we’ll add those as random effects. I think the sites should be coded as
factors, **but I’m not sure what format to use for the date. I think it
should be disease week so that week 18 is treated closer to 19 than 20,
etc., but I’m not totally confident in this.**  
- The family = poisson (link = “log”)… why again?

Load libraries

``` r
library(tidyverse) # for data wrangling
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.2     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(glmmTMB)   # for model fitting
library(DHARMa)    # for residual plots
```

    ## This is DHARMa 0.4.7. For overview type '?DHARMa'. For recent changes, type news(package = 'DHARMa')

``` r
library(emmeans)   # for estimating marginal effects
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
library(multcomp)  # for statistical comparisons on fitted models
```

    ## Loading required package: mvtnorm
    ## Loading required package: survival
    ## Loading required package: TH.data
    ## Loading required package: MASS
    ## 
    ## Attaching package: 'MASS'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select
    ## 
    ## 
    ## Attaching package: 'TH.data'
    ## 
    ## The following object is masked from 'package:MASS':
    ## 
    ##     geyser

## Load data

``` r
## tarsalis datasets from SLCMAD:
tarsalis <- read.csv("../data/tarsalis_2025.csv")

head(tarsalis, 5)
```

    ##   site_code     site_name longitude latitude trap_type collection_date
    ## 1       224 1700 E Church  -111.842 40.72952      GRVD      2025-07-10
    ## 2       224 1700 E Church  -111.842 40.72952      GRVD      2025-07-17
    ## 3       224 1700 E Church  -111.842 40.72952      GRVD      2025-08-21
    ## 4       224 1700 E Church  -111.842 40.72952      GRVD      2025-08-21
    ## 5       224 1700 E Church  -111.842 40.72952      GRVD      2025-09-05
    ##   disease_week        species count season urban_cat
    ## 1           28 Culex tarsalis    NA    mid     urban
    ## 2           29 Culex tarsalis     1    mid     urban
    ## 3           34 Culex tarsalis    NA   late     urban
    ## 4           34 Culex tarsalis    NA   late     urban
    ## 5           36 Culex tarsalis    NA   late     urban

``` r
## pipiens datasets from SLCMAD:
pipiens <- read.csv("../data/pipiens_2025.csv")
head(pipiens, 5)
```

    ##   site_code     site_name longitude latitude trap_type collection_date
    ## 1       224 1700 E Church  -111.842 40.72952      GRVD      2025-05-30
    ## 2       224 1700 E Church  -111.842 40.72952      GRVD      2025-06-05
    ## 3       224 1700 E Church  -111.842 40.72952      GRVD      2025-06-12
    ## 4       224 1700 E Church  -111.842 40.72952      GRVD      2025-06-20
    ## 5       224 1700 E Church  -111.842 40.72952      GRVD      2025-06-26
    ##   disease_week       species count season urban_cat
    ## 1           22 Culex pipiens     3  early    urban 
    ## 2           23 Culex pipiens     9    mid    urban 
    ## 3           24 Culex pipiens    10    mid    urban 
    ## 4           25 Culex pipiens     5    mid    urban 
    ## 5           26 Culex pipiens    30    mid    urban

# GLMM by season and urbanization, each mosquito species separately

Start by binning into discrete seasons: early, mid, late

``` r
table(tarsalis$season, tarsalis$disease_week)
```

    ##        
    ##          15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32
    ##   early  28  33  29  27  32  33  75  79   0   0   0   0   0   0   0   0   0   0
    ##   late    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  89
    ##   mid     0   0   0   0   0   0   0   0  80  82  48  88  82  88  88  76 119   0
    ##        
    ##          33  34  35  36  37  38  39  40
    ##   early   0   0   0   0   0   0   0   0
    ##   late   89  94  48  73  83  90  82  39
    ##   mid     0   0   0   0   0   0   0   0

First, look at the response distribution

``` r
# sub in tarsalis dataset here instead of salamander 
hist(tarsalis$count) 
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/response-1.png)<!-- -->

``` r
hist(log(tarsalis$count))
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/response-2.png)<!-- -->

## Random effects (1 \| site_name/disease_week)

Fit the model with poisson(link=“log”), using random effects (1 \|
site_name/disease_week)

``` r
# tarsalis count data, starting with poisson model:
fit_pois <- glmmTMB(count ~ season*urban_cat + trap_type + (1 | site_name/disease_week),
                    family = poisson(link = "log"),
                    data = tarsalis)
summary(fit_pois)
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## count ~ season * urban_cat + trap_type + (1 | site_name/disease_week)
    ## Data: tarsalis
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##  300226.7  300292.2 -150101.4  300202.7      1721 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                 Name        Variance Std.Dev.
    ##  disease_week:site_name (Intercept) 2.1537   1.4675  
    ##  site_name              (Intercept) 0.2329   0.4826  
    ## Number of obs: 1733, groups:  disease_week:site_name, 1088; site_name, 56
    ## 
    ## Conditional model:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                2.60997    0.21341  12.230  < 2e-16 ***
    ## seasonlate                 2.87337    0.21572  13.320  < 2e-16 ***
    ## seasonmid                  3.16284    0.21480  14.724  < 2e-16 ***
    ## urban_catrural            -0.06776    0.27131  -0.250 0.802791    
    ## urban_caturban            -1.06669    0.36373  -2.933 0.003361 ** 
    ## trap_typeGRVD             -2.75024    0.14330 -19.192  < 2e-16 ***
    ## seasonlate:urban_catrural  0.62952    0.27461   2.292 0.021880 *  
    ## seasonmid:urban_catrural   0.75489    0.27301   2.765 0.005692 ** 
    ## seasonlate:urban_caturban -1.41319    0.37043  -3.815 0.000136 ***
    ## seasonmid:urban_caturban  -1.66376    0.36876  -4.512 6.43e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Look at the residuals:

``` r
simulateResiduals(fit_pois, plot = T)
```

    ## DHARMa:testOutliers with type = binomial may have inflated Type I error rates for integer-valued distributions. To get a more exact result, it is recommended to re-run testOutliers with type = 'bootstrap'. See ?testOutliers for details

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/resids-poisson-log-link-week-1.png)<!-- -->

    ## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
    ##  
    ## Scaled residual values: 0.3757723 0.6631096 0.8546529 0.472019 0.5208264 0.5681077 0.4845933 0.7686813 0.3847694 0.1291647 0.1217306 0.3483495 0.7032616 0.06760113 0.6509805 0.98 0.992 0.992 0.09160973 0.204 ...

## Random effects (1 \| site_name/collection_date)

Fit the model with poisson(link=“log”), using random effects (1 \|
site_name/collection_date)

``` r
# tarsalis count data, starting with poisson model:
fit_pois_date <- glmmTMB(count ~ season*urban_cat + trap_type + (1 | site_name/collection_date),
                    family = poisson(link = "log"),
                    data = tarsalis)
summary(fit_pois_date)
```

    ##  Family: poisson  ( log )
    ## Formula:          
    ## count ~ season * urban_cat + trap_type + (1 | site_name/collection_date)
    ## Data: tarsalis
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##   22644.3   22709.8  -11310.2   22620.3      1721 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups                    Name        Variance Std.Dev.
    ##  collection_date:site_name (Intercept) 2.8366   1.6842  
    ##  site_name                 (Intercept) 0.2325   0.4822  
    ## Number of obs: 1733, groups:  collection_date:site_name, 1668; site_name, 56
    ## 
    ## Conditional model:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                2.77945    0.21334  13.029  < 2e-16 ***
    ## seasonlate                 2.45969    0.20531  11.980  < 2e-16 ***
    ## seasonmid                  2.66287    0.20259  13.144  < 2e-16 ***
    ## urban_catrural            -0.01155    0.27059  -0.043  0.96595    
    ## urban_caturban            -1.19776    0.37425  -3.200  0.00137 ** 
    ## trap_typeGRVD             -2.75428    0.14530 -18.956  < 2e-16 ***
    ## seasonlate:urban_catrural  0.76447    0.26047   2.935  0.00334 ** 
    ## seasonmid:urban_catrural   0.67223    0.25683   2.617  0.00886 ** 
    ## seasonlate:urban_caturban -1.04970    0.37584  -2.793  0.00522 ** 
    ## seasonmid:urban_caturban  -1.21580    0.37190  -3.269  0.00108 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Look at the residuals
simulateResiduals(fit_pois_date, plot = T)
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/poisson-log-link-date-1.png)<!-- -->

    ## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
    ##  
    ## Scaled residual values: 0.451173 0.570916 0.8657983 0.4242838 0.6669795 0.5438227 0.3057994 0.7070638 0.2808567 0.137651 0.06285764 0.2788065 0.6499876 0.1160842 0.576 0.92 0.984 0.992 0.076 0.284 ...

I think it still doesn’t look great… KS test shows significant
deviation, dispersion test shows significant deviation, and there is a
funnel shape. Maybe the negative binomial distribution can help with
this as it helped with the salamanders? But this is all beyond me… give
it a shot!

``` r
fit_nbinom <- glmmTMB(count ~ season*urban_cat + trap_type + (1 | site_name/collection_date),
                    family = nbinom1(link = "log"),
                    data = tarsalis)

simulateResiduals(fit_nbinom, plot = T)
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/nbinom-log-link-date-1.png)<!-- -->

    ## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
    ##  
    ## Scaled residual values: 0.4321881 0.6137061 0.8431486 0.4433061 0.6047028 0.6146553 0.4889945 0.7055691 0.357675 0.158536 0.1110842 0.3167209 0.6085761 0.0709376 0.5811151 0.944 0.968 0.968 0.128 0.28 ...

I think it looks the same!!! So, stick with the poisson log-link? But I
do think that random effects (1 \| site_name/collection_date) looked
better than disease_week. Not sure if that is good justification to use
collection_date though.

## Marginal effects

Sticking with model with poisson(link=“log”), using random effects (1 \|
site_name/collection_date) `fit_pois_date`

glmmTMB(count ~ season\*urban_cat + trap_type + (1 \|
site_name/collection_date), family = nbinom1(link = “log”), data =
tarsalis)

``` r
# First, just look at the effect of season, doing a 1-way pairwise comparison
em_season = emmeans(fit_pois_date, pairwise ~ season, type = "response")
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
em_season
```

    ## $emmeans
    ##  season  rate    SE  df asymp.LCL asymp.UCL
    ##  early   2.72 0.411 Inf      2.02      3.65
    ##  late   28.90 3.400 Inf     22.95     36.40
    ##  mid    32.49 3.810 Inf     25.81     40.89
    ## 
    ## Results are averaged over the levels of: urban_cat, trap_type 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## 
    ## $contrasts
    ##  contrast      ratio     SE  df null z.ratio p.value
    ##  early / late 0.0940 0.0128 Inf    1 -17.356  <.0001
    ##  early / mid  0.0836 0.0113 Inf    1 -18.418  <.0001
    ##  late / mid   0.8895 0.0896 Inf    1  -1.163  0.4756
    ## 
    ## Results are averaged over the levels of: urban_cat, trap_type 
    ## P value adjustment: tukey method for comparing a family of 3 estimates 
    ## Tests are performed on the log scale

``` r
cl = cld(em_season, Letters = letters)

ggplot(data = cl, aes(x = season, y = rate)) +
    geom_point(size=2.5, color="black") +
    geom_errorbar(aes(x=season, ymin = asymp.LCL,
                      ymax = asymp.UCL),
                  width = 0.2, size=1, color="black") +
    geom_text(aes(label = gsub(" ", "", .group)),
              position = position_nudge(x = 0.3)) +
    ggeasy::easy_remove_axes("both", "title")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Next, let's see how the effects vary with different mined levels in a 2-way emmeans comparison
em_spp2 = emmeans(fit_pois_date, pairwise ~ season | urban_cat, type = "response")
em_spp2
```

    ## $emmeans
    ## urban_cat = peri:
    ##  season   rate     SE  df asymp.LCL asymp.UCL
    ##  early    4.06  0.916 Inf     2.613      6.32
    ##  late    47.56  9.250 Inf    32.490     69.62
    ##  mid     58.28 11.200 Inf    40.024     84.85
    ## 
    ## urban_cat = rural:
    ##  season   rate     SE  df asymp.LCL asymp.UCL
    ##  early    4.02  0.730 Inf     2.814      5.74
    ##  late   100.98 16.100 Inf    73.877    138.02
    ##  mid    112.83 17.700 Inf    82.940    153.49
    ## 
    ## urban_cat = urban:
    ##  season   rate     SE  df asymp.LCL asymp.UCL
    ##  early    1.23  0.382 Inf     0.667      2.26
    ##  late     5.03  1.000 Inf     3.398      7.43
    ##  mid      5.22  1.030 Inf     3.537      7.69
    ## 
    ## Results are averaged over the levels of: trap_type 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## 
    ## $contrasts
    ## urban_cat = peri:
    ##  contrast      ratio      SE  df null z.ratio p.value
    ##  early / late 0.0855 0.01750 Inf    1 -11.980  <.0001
    ##  early / mid  0.0697 0.01410 Inf    1 -13.144  <.0001
    ##  late / mid   0.8161 0.13700 Inf    1  -1.210  0.4473
    ## 
    ## urban_cat = rural:
    ##  contrast      ratio      SE  df null z.ratio p.value
    ##  early / late 0.0398 0.00638 Inf    1 -20.110  <.0001
    ##  early / mid  0.0356 0.00562 Inf    1 -21.122  <.0001
    ##  late / mid   0.8950 0.11800 Inf    1  -0.842  0.6771
    ## 
    ## urban_cat = urban:
    ##  contrast      ratio      SE  df null z.ratio p.value
    ##  early / late 0.2441 0.07690 Inf    1  -4.478  <.0001
    ##  early / mid  0.2353 0.07340 Inf    1  -4.638  <.0001
    ##  late / mid   0.9636 0.20600 Inf    1  -0.174  0.9835
    ## 
    ## Results are averaged over the levels of: trap_type 
    ## P value adjustment: tukey method for comparing a family of 3 estimates 
    ## Tests are performed on the log scale

``` r
cl2 <- cld(em_spp2$emmeans, Letters = letters)

ggplot(data = cl2, aes(x = season, y = rate)) +
    geom_point(size=2.5, color="black") +
    geom_errorbar(aes(x=season, ymin = asymp.LCL,
                      ymax = asymp.UCL),
                  width = 0.2, size=1, color="black") +
    facet_wrap(~ urban_cat, #nrow = 2,
               labeller = label_both) +
    geom_text(aes(label = gsub(" ", "", .group)),
              position = position_nudge(x = 0.4)) +
    ggeasy::easy_remove_axes("y", "title")
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->
\# Salamander example

Just to understand a bit better what it looks like when residuals don’t
“look good”:

``` r
# use the salamanders dataset from glmmTMB:
data("Salamanders")
head(Salamanders)
```

    ##   site mined      cover sample        DOP       Wtemp       DOY spp count
    ## 1 VF-1   yes -1.4423172      1 -0.5956834 -1.22937861 -1.497003  GP     0
    ## 2 VF-2   yes  0.2984104      1 -0.5956834  0.08476529 -1.497003  GP     0
    ## 3 VF-3   yes  0.3978806      1 -1.1913668  1.01417627 -1.294467  GP     0
    ## 4  R-1    no -0.4476157      1  0.0000000 -3.02335795 -2.712216  GP     2
    ## 5  R-2    no  0.5968209      1  0.5956834 -0.14434533 -0.686860  GP     2
    ## 6  R-3    no  1.3428470      1  0.5956834 -0.01466007 -0.686860  GP     1

``` r
# we'll model salamander count as a function of species (spp) and mining presence (mined).
# Because the data are grouped into sites, we'll add that as a random effect.
# first, look at the response distribution

hist(Salamanders$count)
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/salamander-1.png)<!-- -->

``` r
hist(log(Salamanders$count))
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/salamander-2.png)<!-- -->

``` r
# count data, so start with poisson model:
fit_pois <- glmmTMB(count ~ spp*mined + (1 | site),
                    family = poisson(link = "log"),
                    data = Salamanders)

summary(fit_pois)
```

    ##  Family: poisson  ( log )
    ## Formula:          count ~ spp * mined + (1 | site)
    ## Data: Salamanders
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##    1940.2    2007.3    -955.1    1910.2       629 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups Name        Variance Std.Dev.
    ##  site   (Intercept) 0.3316   0.5759  
    ## Number of obs: 644, groups:  site, 23
    ## 
    ## Conditional model:
    ##                  Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)       -3.3771     0.7340  -4.601 4.20e-06 ***
    ## sppPR              0.9163     0.8367   1.095 0.273435    
    ## sppDM              2.2513     0.7434   3.028 0.002458 ** 
    ## sppEC-A            0.6931     0.8660   0.800 0.423494    
    ## sppEC-L            1.7047     0.7687   2.218 0.026576 *  
    ## sppDES-L           2.5257     0.7348   3.437 0.000588 ***
    ## sppDF              2.5257     0.7348   3.437 0.000588 ***
    ## minedno            4.1109     0.7587   5.418 6.02e-08 ***
    ## sppPR:minedno     -2.4887     0.8688  -2.864 0.004178 ** 
    ## sppDM:minedno     -2.1526     0.7554  -2.850 0.004377 ** 
    ## sppEC-A:minedno   -1.5279     0.8838  -1.729 0.083853 .  
    ## sppEC-L:minedno   -1.1212     0.7782  -1.441 0.149670    
    ## sppDES-L:minedno  -1.9527     0.7448  -2.622 0.008748 ** 
    ## sppDF:minedno     -2.6674     0.7485  -3.563 0.000366 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# look at the residuals:
simulateResiduals(fit_pois, plot = T)
```

    ## DHARMa:testOutliers with type = binomial may have inflated Type I error rates for integer-valued distributions. To get a more exact result, it is recommended to re-run testOutliers with type = 'bootstrap'. See ?testOutliers for details

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/salamander-3.png)<!-- -->

    ## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
    ##  
    ## Scaled residual values: 0.8971645 0.5376291 0.6865468 0.4771795 0.5695308 0.3930258 0.3732537 0.5270645 0.7092329 0.3035174 0.938414 0.03329872 0.3608276 0.02056911 0.06117746 0.8165105 0.3426504 0.1332761 0.2699177 0.8026385 ...

``` r
# this doesn't look great. It seems like there is a trend in the residuals with
# greater error at higher values. Maybe the negative binomial distribution can help with this
```

``` r
#continue salamander dataset example
# this doesn't look great. It seems like there is a trend in the residuals with
# greater error at higher values. Maybe the negative binomial distribution can help with this

fit_nbinom <- glmmTMB(count ~ spp*mined + (1 | site),
                      family = nbinom1(link = "log"),
                      data = Salamanders)
simulateResiduals(fit_nbinom, plot = T)
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/salamanders-nbinom-1.png)<!-- -->

    ## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
    ##  
    ## Scaled residual values: 0.4589582 0.7493933 0.6690189 0.5948831 0.51507 0.4854772 0.52614 0.6080341 0.8140457 0.3224168 0.3786706 0.05314925 0.4407776 0.0385783 0.07113632 0.942124 0.4041715 0.3563878 0.382026 0.234601 ...

``` r
# this looks much better. Now we can look at the marginal effects:
```

``` r
# First, just look at the effect of species, doing a 1-way pairwise comparison
em_spp = emmeans(fit_nbinom, pairwise ~ spp, type = "response")
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
em_spp
```

    ## $emmeans
    ##  spp   response     SE  df asymp.LCL asymp.UCL
    ##  GP       0.289 0.1500 Inf     0.105     0.797
    ##  PR       0.225 0.0776 Inf     0.114     0.442
    ##  DM       0.965 0.2020 Inf     0.640     1.455
    ##  EC-A     0.314 0.1040 Inf     0.164     0.600
    ##  EC-L     1.017 0.2160 Inf     0.671     1.543
    ##  DES-L    1.311 0.2530 Inf     0.898     1.914
    ##  DF       1.015 0.1970 Inf     0.694     1.484
    ## 
    ## Results are averaged over the levels of: mined 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## 
    ## $contrasts
    ##  contrast         ratio     SE  df null z.ratio p.value
    ##  GP / PR          1.287 0.7720 Inf    1   0.421  0.9996
    ##  GP / DM          0.300 0.1600 Inf    1  -2.254  0.2667
    ##  GP / (EC-A)      0.921 0.5450 Inf    1  -0.139  1.0000
    ##  GP / (EC-L)      0.284 0.1520 Inf    1  -2.351  0.2200
    ##  GP / (DES-L)     0.221 0.1170 Inf    1  -2.859  0.0643
    ##  GP / DF          0.285 0.1500 Inf    1  -2.377  0.2082
    ##  PR / DM          0.233 0.0861 Inf    1  -3.941  0.0016
    ##  PR / (EC-A)      0.715 0.3200 Inf    1  -0.749  0.9895
    ##  PR / (EC-L)      0.221 0.0818 Inf    1  -4.079  0.0009
    ##  PR / (DES-L)     0.171 0.0619 Inf    1  -4.884  <.0001
    ##  PR / DF          0.221 0.0797 Inf    1  -4.190  0.0006
    ##  DM / (EC-A)      3.073 1.0900 Inf    1   3.154  0.0269
    ##  DM / (EC-L)      0.948 0.2390 Inf    1  -0.211  1.0000
    ##  DM / (DES-L)     0.736 0.1750 Inf    1  -1.290  0.8566
    ##  DM / DF          0.951 0.2250 Inf    1  -0.214  1.0000
    ##  (EC-A) / (EC-L)  0.309 0.1100 Inf    1  -3.296  0.0170
    ##  (EC-A) / (DES-L) 0.240 0.0831 Inf    1  -4.119  0.0008
    ##  (EC-A) / DF      0.309 0.1070 Inf    1  -3.391  0.0124
    ##  (EC-L) / (DES-L) 0.776 0.1850 Inf    1  -1.061  0.9394
    ##  (EC-L) / DF      1.003 0.2380 Inf    1   0.011  1.0000
    ##  (DES-L) / DF     1.291 0.2870 Inf    1   1.149  0.9127
    ## 
    ## Results are averaged over the levels of: mined 
    ## P value adjustment: tukey method for comparing a family of 7 estimates 
    ## Tests are performed on the log scale

``` r
cl = cld(em_spp, Letters = letters)

ggplot(data = cl, aes(x = spp, y = response)) +
  geom_point(size=2.5, color="black") +
  geom_errorbar(aes(x=spp, ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.2, size=1, color="black") +
  geom_text(aes(label = gsub(" ", "", .group)),
            position = position_nudge(x = 0.3)) +
  ggeasy::easy_remove_axes("both", "title")
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/salamanders-effect-1.png)<!-- -->

``` r
# Next, let's see how the effects vary with different mined levels in a 2-way emmeans comparison
em_spp2 = emmeans(fit_nbinom, pairwise ~ spp | mined, type = "response")
em_spp2
```

    ## $emmeans
    ## mined = yes:
    ##  spp   response     SE  df asymp.LCL asymp.UCL
    ##  GP      0.0369 0.0374 Inf   0.00506     0.269
    ##  PR      0.1099 0.0662 Inf   0.03379     0.358
    ##  DM      0.3843 0.1400 Inf   0.18835     0.784
    ##  EC-A    0.1100 0.0661 Inf   0.03391     0.357
    ##  EC-L    0.3241 0.1220 Inf   0.15474     0.679
    ##  DES-L   0.4918 0.1640 Inf   0.25553     0.947
    ##  DF      0.5565 0.1770 Inf   0.29881     1.036
    ## 
    ## mined = no:
    ##  spp   response     SE  df asymp.LCL asymp.UCL
    ##  GP      2.2677 0.4580 Inf   1.52691     3.368
    ##  PR      0.4587 0.1520 Inf   0.24003     0.877
    ##  DM      2.4227 0.4840 Inf   1.63813     3.583
    ##  EC-A    0.8958 0.2370 Inf   0.53302     1.506
    ##  EC-L    3.1942 0.6060 Inf   2.20256     4.632
    ##  DES-L   3.4930 0.6520 Inf   2.42296     5.036
    ##  DF      1.8510 0.3950 Inf   1.21830     2.812
    ## 
    ## Confidence level used: 0.95 
    ## Intervals are back-transformed from the log scale 
    ## 
    ## $contrasts
    ## mined = yes:
    ##  contrast          ratio     SE  df null z.ratio p.value
    ##  GP / PR          0.3353 0.3870 Inf    1  -0.948  0.9647
    ##  GP / DM          0.0959 0.1010 Inf    1  -2.232  0.2782
    ##  GP / (EC-A)      0.3350 0.3860 Inf    1  -0.949  0.9645
    ##  GP / (EC-L)      0.1138 0.1200 Inf    1  -2.065  0.3737
    ##  GP / (DES-L)     0.0750 0.0780 Inf    1  -2.491  0.1624
    ##  GP / DF          0.0663 0.0685 Inf    1  -2.624  0.1186
    ##  PR / DM          0.2861 0.1890 Inf    1  -1.889  0.4874
    ##  PR / (EC-A)      0.9992 0.8140 Inf    1  -0.001  1.0000
    ##  PR / (EC-L)      0.3392 0.2260 Inf    1  -1.624  0.6666
    ##  PR / (DES-L)     0.2235 0.1440 Inf    1  -2.320  0.2342
    ##  PR / DF          0.1976 0.1260 Inf    1  -2.550  0.1418
    ##  DM / (EC-A)      3.4922 2.3100 Inf    1   1.892  0.4857
    ##  DM / (EC-L)      1.1857 0.5530 Inf    1   0.365  0.9998
    ##  DM / (DES-L)     0.7813 0.3410 Inf    1  -0.566  0.9977
    ##  DM / DF          0.6906 0.2920 Inf    1  -0.876  0.9761
    ##  (EC-A) / (EC-L)  0.3395 0.2260 Inf    1  -1.625  0.6662
    ##  (EC-A) / (DES-L) 0.2237 0.1440 Inf    1  -2.323  0.2330
    ##  (EC-A) / DF      0.1977 0.1260 Inf    1  -2.551  0.1414
    ##  (EC-L) / (DES-L) 0.6590 0.2910 Inf    1  -0.943  0.9656
    ##  (EC-L) / DF      0.5824 0.2490 Inf    1  -1.262  0.8692
    ##  (DES-L) / DF     0.8838 0.3500 Inf    1  -0.312  0.9999
    ## 
    ## mined = no:
    ##  contrast          ratio     SE  df null z.ratio p.value
    ##  GP / PR          4.9434 1.6300 Inf    1   4.844  <.0001
    ##  GP / DM          0.9360 0.1880 Inf    1  -0.329  0.9999
    ##  GP / (EC-A)      2.5313 0.6720 Inf    1   3.500  0.0085
    ##  GP / (EC-L)      0.7099 0.1380 Inf    1  -1.757  0.5772
    ##  GP / (DES-L)     0.6492 0.1230 Inf    1  -2.278  0.2545
    ##  GP / DF          1.2251 0.2630 Inf    1   0.944  0.9653
    ##  PR / DM          0.1893 0.0623 Inf    1  -5.057  <.0001
    ##  PR / (EC-A)      0.5121 0.1900 Inf    1  -1.802  0.5470
    ##  PR / (EC-L)      0.1436 0.0468 Inf    1  -5.961  <.0001
    ##  PR / (DES-L)     0.1313 0.0424 Inf    1  -6.290  <.0001
    ##  PR / DF          0.2478 0.0836 Inf    1  -4.136  0.0007
    ##  DM / (EC-A)      2.7044 0.7140 Inf    1   3.766  0.0031
    ##  DM / (EC-L)      0.7585 0.1460 Inf    1  -1.433  0.7842
    ##  DM / (DES-L)     0.6936 0.1300 Inf    1  -1.948  0.4485
    ##  DM / DF          1.3089 0.2790 Inf    1   1.261  0.8696
    ##  (EC-A) / (EC-L)  0.2805 0.0727 Inf    1  -4.905  <.0001
    ##  (EC-A) / (DES-L) 0.2565 0.0656 Inf    1  -5.316  <.0001
    ##  (EC-A) / DF      0.4840 0.1330 Inf    1  -2.646  0.1125
    ##  (EC-L) / (DES-L) 0.9145 0.1650 Inf    1  -0.496  0.9989
    ##  (EC-L) / DF      1.7257 0.3550 Inf    1   2.651  0.1111
    ##  (DES-L) / DF     1.8871 0.3820 Inf    1   3.136  0.0285
    ## 
    ## P value adjustment: tukey method for comparing a family of 7 estimates 
    ## Tests are performed on the log scale

``` r
cl2 = cld(em_spp2, Letters = letters)

ggplot(data = cl2, aes(x = spp, y = response)) +
  geom_point(size=2.5, color="black") +
  geom_errorbar(aes(x=spp, ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.2, size=1, color="black") +
  facet_wrap(~ mined, #nrow = 2,
             labeller = label_both) +
  geom_text(aes(label = gsub(" ", "", .group)),
            position = position_nudge(x = 0.4)) +
  ggeasy::easy_remove_axes("y", "title")
```

![](02_GLMM_pip_tars_SLC2025_files/figure-gfm/salamanders-effect-2.png)<!-- -->
