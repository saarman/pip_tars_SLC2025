GAM_Cx.pip/tar_abundance_SLC2025
================
Katie Graybeal
2026-06-03

- [Setup](#setup)
- [Prepare Data](#prepare-data)
  - [Combined data](#combined-data)
- [GAM with each species separately, count by urbanization + trap_type
  (random effect = site
  name)](#gam-with-each-species-separately-count-by-urbanization--trap_type-random-effect--site-name)
  - [Fit model](#fit-model)

# Setup

**Research Topic:** testing whether habitat and seasonal partitioning
between Culex pipiens s.l. and Culex tarsalis shapes West Nile Virus
(WNV) dynamics across urban–rural gradients.

**Core hypothesis:** early/mid-season amplification dominated by pipiens
in urban areas, later spillover involving tarsalis moving into
urban/peri-urban areas.

**Approach:** Preliminary results visualized via mapping, with species
identity and abundance as primary response variables. Model mosquito
abundance and proportions using GLMM with GAM smoothing:

count ~ season\*urbanization + trap_type + (1\|site/date), family =
poisson(link = “log”):

- Response variable = mosquito abundance  
- Predictors = season\*urbanization  
- The trap type could be important, so we will add that as a fixed
  effect (covariate)… is this correct? We do think that the response
  variable of count of mosquitoes depends on trap type, since tarsalis
  seems to be more attracted to CO2 than pipiens, and we want to
  quantify that effect. Note that poisson model does not give a fixed
  offset (due to the log link)… The structure of this model means that
  it will estimate an effect that scales with the total number of
  mosquitos caught, which is exactly what we want.  
- The data are grouped into sites and are also linked through time, so
  we’ll add those as random effects. I think the sites should be coded
  as factors, **but I’m not sure what format to use for the date. I
  think it should be disease week so that week 18 is treated closer to
  19 than 20, etc., but I’m not totally confident in this.**
- The family = poisson (link = “log”)… why again?

**For simple model:** count ~ disease_week\*urbanization +
(1\|site/date), family = poisson(link = “log”)

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
library(mgcViz)    # for residual plots
```

    ## Loading required package: mgcv
    ## Loading required package: nlme
    ## 
    ## Attaching package: 'nlme'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     collapse
    ## 
    ## This is mgcv 1.9-4. For overview type '?mgcv'.
    ## Loading required package: qgam
    ## Registered S3 method overwritten by 'mgcViz':
    ##   method from   
    ##   +.gg   ggplot2
    ## 
    ## Attaching package: 'mgcViz'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     qqline, qqnorm, qqplot

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

``` r
library(dplyr)     # for mutating dataframe to change labels in dataset
library(mgcv)      # fits GAM
library(broom)
library(ggplot2)
```

# Prepare Data

## Combined data

``` r
## tarsalis datasets from SLCMAD:
tarsalis <- read.csv("../data/tarsalis_2025.csv")
## pipiens datasets from SLCMAD:
pipiens <- read.csv("../data/pipiens_2025.csv")
## combine

combined <- bind_rows(tarsalis, pipiens)

## Set factor levels
combined <- combined %>%
  mutate(
    species = factor(
      species,
      levels = c("Culex pipiens", "Culex tarsalis")
    ),
    urban_cat = trimws(tolower(urban_cat)),
    urbanization = factor(
      urban_cat,
      levels = c("rural", "peri", "urban")
    ),
    season = factor(
      season,
      levels = c("early", "mid", "late")
    )
  )

#check
table(combined$species)
```

    ## 
    ##  Culex pipiens Culex tarsalis 
    ##           1394           1774

``` r
table(combined$season, combined$species)
```

    ##        
    ##         Culex pipiens Culex tarsalis
    ##   early           193            336
    ##   mid             653            751
    ##   late            548            687

``` r
combined <- combined %>%
  mutate(
    species = factor(species),
    urbanization = factor(urbanization),
    trap_type = factor(trap_type),
    site_name = factor(site_name),
    disease_week = as.numeric(disease_week)
  )
```

# GAM with each species separately, count by urbanization + trap_type (random effect = site name)

## Fit model

``` r
# pull out pipiens
pipiens <- combined[combined$species == "Culex pipiens",]
head(pipiens$site_name)
```

    ## [1] 1700 E Church 1700 E Church 1700 E Church 1700 E Church 1700 E Church
    ## [6] 1700 E Church
    ## 59 Levels: 1700 E Church 300 E Church 700 S 200 W ... Wingpointe

``` r
# Fit GAM model with site_name random effect
pip_gam <- gam(
  count ~ urbanization + trap_type + 
    s(disease_week, 
      by = urbanization, 
      bs = "fs",  # bs = "fs" for independent smooths
      k=10, 
      m=3) +  # m = 3 restricts wigglyness
    s(site_name, bs = "re"),
  family = nb(),
  data = pipiens,
  method = "REML"
)

# pull out tarsalis
tarsalis <- combined[combined$species == "Culex tarsalis",]
head(tarsalis$site_name)
```

    ## [1] 1700 E Church 1700 E Church 1700 E Church 1700 E Church 1700 E Church
    ## [6] 1700 E Church
    ## 59 Levels: 1700 E Church 300 E Church 700 S 200 W ... Wingpointe

``` r
# Fit GAM model with site_name only
tar_gam <- gam(
  count ~ urbanization + trap_type + 
    s(disease_week, 
      by = urbanization, 
      bs = "fs",  # bs = "fs" for independent smooths
      k=20, 
      m=3) +
    s(site_name, bs = "re"),
  family = nb(),
  data = tarsalis,
  method = "REML"
)
```

\#Splitting up tarsalis by trap type

``` r
table(tarsalis$urbanization,tarsalis$trap_type)
```

    ##        
    ##         CO2 GRVD
    ##   rural 862    0
    ##   peri  529    0
    ##   urban 273  110

``` r
tarsalis_CO2 <- combined %>%
  filter(species == "Culex tarsalis",
         trap_type == "CO2")

tarsalis_GRVD <- combined %>%
  filter(species == "Culex tarsalis",
         trap_type == "GRVD")

head(tarsalis_CO2)
```

    ##   site_code  site_name longitude latitude trap_type collection_date
    ## 1        27 Ambassador -112.0292 40.84601       CO2      2025-04-11
    ## 2        27 Ambassador -112.0292 40.84601       CO2      2025-04-16
    ## 3        27 Ambassador -112.0292 40.84601       CO2      2025-04-23
    ## 4        27 Ambassador -112.0292 40.84601       CO2      2025-05-02
    ## 5        27 Ambassador -112.0292 40.84601       CO2      2025-05-08
    ## 6        27 Ambassador -112.0292 40.84601       CO2      2025-05-15
    ##   disease_week        species count season urban_cat urbanization
    ## 1           15 Culex tarsalis     7  early     rural        rural
    ## 2           16 Culex tarsalis     2  early     rural        rural
    ## 3           17 Culex tarsalis     1  early     rural        rural
    ## 4           18 Culex tarsalis     6  early     rural        rural
    ## 5           19 Culex tarsalis    28  early     rural        rural
    ## 6           20 Culex tarsalis     1  early     rural        rural

``` r
head(tarsalis_GRVD)
```

    ##   site_code     site_name longitude latitude trap_type collection_date
    ## 1       224 1700 E Church  -111.842 40.72952      GRVD      2025-07-10
    ## 2       224 1700 E Church  -111.842 40.72952      GRVD      2025-07-17
    ## 3       224 1700 E Church  -111.842 40.72952      GRVD      2025-08-21
    ## 4       224 1700 E Church  -111.842 40.72952      GRVD      2025-08-21
    ## 5       224 1700 E Church  -111.842 40.72952      GRVD      2025-09-05
    ## 6       224 1700 E Church  -111.842 40.72952      GRVD      2025-09-18
    ##   disease_week        species count season urban_cat urbanization
    ## 1           28 Culex tarsalis    NA    mid     urban        urban
    ## 2           29 Culex tarsalis     1    mid     urban        urban
    ## 3           34 Culex tarsalis    NA   late     urban        urban
    ## 4           34 Culex tarsalis    NA   late     urban        urban
    ## 5           36 Culex tarsalis    NA   late     urban        urban
    ## 6           38 Culex tarsalis     2   late     urban        urban

``` r
#tarsalis_C02 gam
gam_tarsalis_CO2 <- gam(
  count ~ urbanization +
    s(disease_week, k = 10, m = 2) +
    s(site_name, bs = "re"),
  family = nb(),
  data = tarsalis_CO2,
  method = "REML"
)

table(tarsalis_CO2$urbanization)
```

    ## 
    ## rural  peri urban 
    ##   862   529   273

``` r
table(tarsalis_CO2$site_name)
```

    ## 
    ##                  1700 E Church                   300 E Church 
    ##                              0                              0 
    ##                    700 S 200 W                     Ambassador 
    ##                              0                             42 
    ##       Amelia Earhart (F.S. #9)                            ATV 
    ##                             40                             45 
    ##                        Audubon                         Avocet 
    ##                             38                             38 
    ##               Bird Reclamation                      Blackhawk 
    ##                             44                             41 
    ##                      Blue Lake              Canyon Rim Church 
    ##                             36                              0 
    ##                Deseret Nursery                 Downington Ave 
    ##                             43                             20 
    ##                       Drechsel                Elephant Island 
    ##                             45                             36 
    ##                      Fair Park                 Fire Station 1 
    ##                             41                              0 
    ##                Fire Station 13                 Fire Station 2 
    ##                             13                             14 
    ##                 Fire Station 3                 Fire Station 4 
    ##                              0                             17 
    ##                 Fire Station 5                 Fire Station 6 
    ##                             18                             19 
    ##                 Fire Station 8                Frontage Corral 
    ##                             19                             39 
    ##           Glendale Golf Course                      Goat Farm 
    ##                             37                             42 
    ##                         Goggin               Graystone Church 
    ##                             44                              0 
    ##                       Harrison                       Hinckley 
    ##                             44                             44 
    ##                      Hogle Zoo Indian Hills Elementary School 
    ##                             13                              0 
    ##                     Inland Sea                           Jade 
    ##                             42                             41 
    ##                     KETO Pumps                      Lakefront 
    ##                             39                             40 
    ##   Lee Kay Dog Training Grounds                      New State 
    ##                             39                             40 
    ##             Nibley Golf Course                         Noggin 
    ##                             20                             40 
    ##                     Northpoint                      O'Reillys 
    ##                             41                             37 
    ##                Power Pole Road                         Prison 
    ##                             39                             44 
    ##                     RAC Soccer               Rose Park Church 
    ##                             44                              0 
    ##          Rose Park Golf Course                           Rudy 
    ##                             42                             45 
    ##                         Runway                  Senior Center 
    ##                             29                              0 
    ##                SLCMAD Facility                 South Saratoga 
    ##                             42                             40 
    ##                 Stratford Ward                   Tracy Aviary 
    ##                              0                              0 
    ##                   Train Tracks                   White Church 
    ##                             44                              0 
    ##                     Wingpointe 
    ##                             44

``` r
str(tarsalis_CO2)
```

    ## 'data.frame':    1664 obs. of  12 variables:
    ##  $ site_code      : int  27 27 27 27 27 27 27 27 27 27 ...
    ##  $ site_name      : Factor w/ 59 levels "1700 E Church",..: 4 4 4 4 4 4 4 4 4 4 ...
    ##  $ longitude      : num  -112 -112 -112 -112 -112 ...
    ##  $ latitude       : num  40.8 40.8 40.8 40.8 40.8 ...
    ##  $ trap_type      : Factor w/ 2 levels "CO2","GRVD": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ collection_date: chr  "2025-04-11" "2025-04-16" "2025-04-23" "2025-05-02" ...
    ##  $ disease_week   : num  15 16 17 18 19 20 21 21 22 22 ...
    ##  $ species        : Factor w/ 2 levels "Culex pipiens",..: 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ count          : int  7 2 1 6 28 1 19 248 660 776 ...
    ##  $ season         : Factor w/ 3 levels "early","mid",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ urban_cat      : chr  "rural" "rural" "rural" "rural" ...
    ##  $ urbanization   : Factor w/ 3 levels "rural","peri",..: 1 1 1 1 1 1 1 1 1 1 ...

\#splitting up pipiens by trap type

``` r
table(pipiens$urbanization,pipiens$trap_type)
```

    ##        
    ##         CO2 GRVD
    ##   rural 357    0
    ##   peri  401    0
    ##   urban 244  392

``` r
pipiens_CO2 <- combined %>%
  filter(species == "Culex pipiens",
         trap_type == "CO2")

#pipiens_CO2 gam
gam_pipiens_CO2 <- gam(
  count ~ urbanization +
    s(disease_week, k = 10, m = 2) +
    s(site_name, bs = "re"),
  family = nb(),
  data = pipiens_CO2,
  method = "REML"
)


#pipiens_GRVD <- combined %>%
  #filter(species == "Culex pipiens",
        # trap_type == "GRVD")

head(pipiens_CO2)
```

    ##   site_code  site_name longitude latitude trap_type collection_date
    ## 1        27 Ambassador -112.0292 40.84601       CO2      2025-04-11
    ## 2        27 Ambassador -112.0292 40.84601       CO2      2025-04-23
    ## 3        27 Ambassador -112.0292 40.84601       CO2      2025-05-15
    ## 4        27 Ambassador -112.0292 40.84601       CO2      2025-05-28
    ## 5        27 Ambassador -112.0292 40.84601       CO2      2025-06-10
    ## 6        27 Ambassador -112.0292 40.84601       CO2      2025-06-18
    ##   disease_week       species count season urban_cat urbanization
    ## 1           15 Culex pipiens     1  early     rural        rural
    ## 2           17 Culex pipiens     1  early     rural        rural
    ## 3           20 Culex pipiens     1  early     rural        rural
    ## 4           22 Culex pipiens    76  early     rural        rural
    ## 5           24 Culex pipiens    79    mid     rural        rural
    ## 6           25 Culex pipiens     5    mid     rural        rural

``` r
#head(pipiens_GRVD)

table(pipiens_CO2$urbanization)
```

    ## 
    ## rural  peri urban 
    ##   357   401   244

``` r
table(pipiens_CO2$site_name)
```

    ## 
    ##                  1700 E Church                   300 E Church 
    ##                              0                              0 
    ##                    700 S 200 W                     Ambassador 
    ##                              0                             14 
    ##       Amelia Earhart (F.S. #9)                            ATV 
    ##                             32                             36 
    ##                        Audubon                         Avocet 
    ##                             12                             11 
    ##               Bird Reclamation                      Blackhawk 
    ##                             26                             17 
    ##                      Blue Lake              Canyon Rim Church 
    ##                             16                              0 
    ##                Deseret Nursery                 Downington Ave 
    ##                             37                             20 
    ##                       Drechsel                Elephant Island 
    ##                             32                             15 
    ##                      Fair Park                 Fire Station 1 
    ##                             28                              0 
    ##                Fire Station 13                 Fire Station 2 
    ##                             12                             15 
    ##                 Fire Station 3                 Fire Station 4 
    ##                              0                             15 
    ##                 Fire Station 5                 Fire Station 6 
    ##                             18                             17 
    ##                 Fire Station 8                Frontage Corral 
    ##                             20                             22 
    ##           Glendale Golf Course                      Goat Farm 
    ##                             36                             24 
    ##                         Goggin               Graystone Church 
    ##                             22                              0 
    ##                       Harrison                       Hinckley 
    ##                             13                             29 
    ##                      Hogle Zoo Indian Hills Elementary School 
    ##                             13                              0 
    ##                     Inland Sea                           Jade 
    ##                             14                             13 
    ##                     KETO Pumps                      Lakefront 
    ##                             36                             13 
    ##   Lee Kay Dog Training Grounds                      New State 
    ##                             28                             12 
    ##             Nibley Golf Course                         Noggin 
    ##                             18                             16 
    ##                     Northpoint                      O'Reillys 
    ##                             17                             30 
    ##                Power Pole Road                         Prison 
    ##                             20                             20 
    ##                     RAC Soccer               Rose Park Church 
    ##                             42                              0 
    ##          Rose Park Golf Course                           Rudy 
    ##                             32                             17 
    ##                         Runway                  Senior Center 
    ##                             10                              0 
    ##                SLCMAD Facility                 South Saratoga 
    ##                             25                             13 
    ##                 Stratford Ward                   Tracy Aviary 
    ##                              0                              0 
    ##                   Train Tracks                   White Church 
    ##                             36                              0 
    ##                     Wingpointe 
    ##                             38

``` r
str(pipiens_CO2)
```

    ## 'data.frame':    1002 obs. of  12 variables:
    ##  $ site_code      : int  27 27 27 27 27 27 27 27 27 27 ...
    ##  $ site_name      : Factor w/ 59 levels "1700 E Church",..: 4 4 4 4 4 4 4 4 4 4 ...
    ##  $ longitude      : num  -112 -112 -112 -112 -112 ...
    ##  $ latitude       : num  40.8 40.8 40.8 40.8 40.8 ...
    ##  $ trap_type      : Factor w/ 2 levels "CO2","GRVD": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ collection_date: chr  "2025-04-11" "2025-04-23" "2025-05-15" "2025-05-28" ...
    ##  $ disease_week   : num  15 17 20 22 24 25 26 26 27 28 ...
    ##  $ species        : Factor w/ 2 levels "Culex pipiens",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ count          : int  1 1 1 76 79 5 231 7 33 69 ...
    ##  $ season         : Factor w/ 3 levels "early","mid",..: 1 1 1 1 2 2 2 2 2 2 ...
    ##  $ urban_cat      : chr  "rural" "rural" "rural" "rural" ...
    ##  $ urbanization   : Factor w/ 3 levels "rural","peri",..: 1 1 1 1 1 1 1 1 1 1 ...

\#Compare models on the same dataset. If your scientific question is
“Does trap type matter?”, fit both models to the full dataset:

``` r
m1 <- gam(
  count ~ species * urbanization +
    s(disease_week, by = species, bs = "fs", k = 10, m = 2) +
    s(site_name, bs = "re"),
  family = nb(),
  data = combined,
  method = "REML"
)

m2 <- gam(
  count ~ species * urbanization + trap_type +
    s(disease_week, by = species, bs = "fs", k = 10, m = 2) +
    s(site_name, bs = "re"),
  family = nb(),
  data = combined,
  method = "REML"
)

AIC(m1, m2)
```

    ##          df      AIC
    ## m1 75.70741 34615.28
    ## m2 76.13194 34570.08

``` r
#Likelihood ratio test (nested models) If both models are fit to the same data and one is nested within the other.
anova(m1, m2, test = "Chisq")
```

    ## Analysis of Deviance Table
    ## 
    ## Model 1: count ~ species * urbanization + s(disease_week, by = species, 
    ##     bs = "fs", k = 10, m = 2) + s(site_name, bs = "re")
    ## Model 2: count ~ species * urbanization + trap_type + s(disease_week, 
    ##     by = species, bs = "fs", k = 10, m = 2) + s(site_name, bs = "re")
    ##   Resid. Df Resid. Dev      Df Deviance Pr(>Chi)    
    ## 1    3036.1      34464                              
    ## 2    3035.2      34418 0.86367    46.05 8.02e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

\#Compare Cx.tarsalis abundance by trap type Urban

``` r
#Analyze Urban sites only:If your question is simply:

"Within Urban sites, does trap type affect abundance?fit the model on Urban observations only"
```

    ## [1] "Within Urban sites, does trap type affect abundance?fit the model on Urban observations only"

``` r
urban_dat <- subset(tarsalis, urbanization == "urban")
nrow(urban_dat)
```

    ## [1] 383

``` r
summary(urban_dat)
```

    ##    site_code                     site_name     longitude         latitude    
    ##  Min.   : 31.0   Rose Park Golf Course: 42   Min.   :-111.9   Min.   :40.71  
    ##  1st Qu.: 33.0   Fair Park            : 41   1st Qu.:-111.9   1st Qu.:40.73  
    ##  Median :213.0   Glendale Golf Course : 37   Median :-111.9   Median :40.75  
    ##  Mean   :159.7   Fire Station 8       : 27   Mean   :-111.9   Mean   :40.75  
    ##  3rd Qu.:219.0   Downington Ave       : 25   3rd Qu.:-111.9   3rd Qu.:40.78  
    ##  Max.   :231.0   Nibley Golf Course   : 25   Max.   :-111.8   Max.   :40.80  
    ##                  (Other)              :186                                   
    ##  trap_type  collection_date     disease_week             species   
    ##  CO2 :273   Length:383         Min.   :15.00   Culex pipiens :  0  
    ##  GRVD:110   Class :character   1st Qu.:26.00   Culex tarsalis:383  
    ##             Mode  :character   Median :31.00                       
    ##                                Mean   :30.38                       
    ##                                3rd Qu.:35.00                       
    ##                                Max.   :40.00                       
    ##                                                                    
    ##      count          season     urban_cat         urbanization
    ##  Min.   :  1.00   early: 44   Length:383         rural:  0   
    ##  1st Qu.:  2.00   mid  :163   Class :character   peri :  0   
    ##  Median :  9.00   late :176   Mode  :character   urban:383   
    ##  Mean   : 46.55                                              
    ##  3rd Qu.: 42.25                                              
    ##  Max.   :797.00                                              
    ##  NA's   :41

``` r
urban_gam <- mgcv::gam(
  count ~ trap_type +
    s(disease_week, k = 5) +
    s(site_name, bs = "re"),
  family = mgcv::nb(),
  data = urban_dat,
  method = "REML"
)

summary(urban_gam)
```

    ## 
    ## Family: Negative Binomial(1.046) 
    ## Link function: log 
    ## 
    ## Formula:
    ## count ~ trap_type + s(disease_week, k = 5) + s(site_name, bs = "re")
    ## 
    ## Parametric coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     3.3231     0.1824   18.22   <2e-16 ***
    ## trap_typeGRVD  -2.7347     0.2053  -13.32   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                    edf Ref.df Chi.sq p-value    
    ## s(disease_week)  3.647  3.926  99.63  <2e-16 ***
    ## s(site_name)    14.558 21.000 164.91  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.413   Deviance explained = 63.7%
    ## -REML = 1402.2  Scale est. = 1         n = 342

``` r
#Then test the trap effect:

#drop1(urban_gam, test = "Chisq")
```

\#Compare Cx.tarsalis abundance by trap type Urban

``` r
#Analyze Urban sites only:If your question is simply:
"Within Urban sites, does trap type affect abundance?fit the model on Urban observations only"
```

    ## [1] "Within Urban sites, does trap type affect abundance?fit the model on Urban observations only"

``` r
urban_dat <- subset(tarsalis, urbanization == "urban")
nrow(urban_dat)
```

    ## [1] 383

``` r
summary(urban_dat)
```

    ##    site_code                     site_name     longitude         latitude    
    ##  Min.   : 31.0   Rose Park Golf Course: 42   Min.   :-111.9   Min.   :40.71  
    ##  1st Qu.: 33.0   Fair Park            : 41   1st Qu.:-111.9   1st Qu.:40.73  
    ##  Median :213.0   Glendale Golf Course : 37   Median :-111.9   Median :40.75  
    ##  Mean   :159.7   Fire Station 8       : 27   Mean   :-111.9   Mean   :40.75  
    ##  3rd Qu.:219.0   Downington Ave       : 25   3rd Qu.:-111.9   3rd Qu.:40.78  
    ##  Max.   :231.0   Nibley Golf Course   : 25   Max.   :-111.8   Max.   :40.80  
    ##                  (Other)              :186                                   
    ##  trap_type  collection_date     disease_week             species   
    ##  CO2 :273   Length:383         Min.   :15.00   Culex pipiens :  0  
    ##  GRVD:110   Class :character   1st Qu.:26.00   Culex tarsalis:383  
    ##             Mode  :character   Median :31.00                       
    ##                                Mean   :30.38                       
    ##                                3rd Qu.:35.00                       
    ##                                Max.   :40.00                       
    ##                                                                    
    ##      count          season     urban_cat         urbanization
    ##  Min.   :  1.00   early: 44   Length:383         rural:  0   
    ##  1st Qu.:  2.00   mid  :163   Class :character   peri :  0   
    ##  Median :  9.00   late :176   Mode  :character   urban:383   
    ##  Mean   : 46.55                                              
    ##  3rd Qu.: 42.25                                              
    ##  Max.   :797.00                                              
    ##  NA's   :41

``` r
urban_gam <- mgcv::gam(
  count ~ trap_type +
    s(disease_week, k = 5) +
    s(site_name, bs = "re"),
  family = mgcv::nb(),
  data = urban_dat,
  method = "REML"
)
summary(urban_gam)
```

    ## 
    ## Family: Negative Binomial(1.046) 
    ## Link function: log 
    ## 
    ## Formula:
    ## count ~ trap_type + s(disease_week, k = 5) + s(site_name, bs = "re")
    ## 
    ## Parametric coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     3.3231     0.1824   18.22   <2e-16 ***
    ## trap_typeGRVD  -2.7347     0.2053  -13.32   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Approximate significance of smooth terms:
    ##                    edf Ref.df Chi.sq p-value    
    ## s(disease_week)  3.647  3.926  99.63  <2e-16 ***
    ## s(site_name)    14.558 21.000 164.91  <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## R-sq.(adj) =  0.413   Deviance explained = 63.7%
    ## -REML = 1402.2  Scale est. = 1         n = 342

``` r
#Then test the trap effect:
#drop1(urban_gam, test = "Chisq")
```
