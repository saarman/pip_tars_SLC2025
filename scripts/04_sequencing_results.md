Oxford nanopore results for CT and UT
================
Norah Saarman
2026-05-13

- [Setup](#setup)
- [Map for UT replicating the one for
  CT](#map-for-ut-replicating-the-one-for-ct)

# Setup

``` r
library(tidyverse)
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
# -------------------------
# Connecticut
# -------------------------
ct <- tribble(
  ~City, ~Code, ~N, ~ITS2, ~kdr, ~pip_M, ~pip_its2, ~rest_M, ~rest_its2,
  "Bridgeport",   "BP39", 128, 122, 118, 109, 119, 15,  2,
  "Glastonbury",  "GL77", 118, 117, 115,  23, 111, 44,  4,
  "Hartford",     "HA78", 225, 221, 221, 177, 209, 32, 11,
  "Killingworth", "K27",   30,  29,  28,   8,   1, 22, 27,
  "New Haven",    "NH54", 174, 168, 170, 110, 164, 30,  2,
  "Salem",        "SA105", 60,  57,  57,  25,  23, 23, 32,
  "West Haven",   "WH104",195, 189, 185, 123, 161, 34, 25,
  "Westport",     "WP69", 134, 129, 130,  99,  84, 31, 43
)


ct_plot <- ct %>%
  mutate(
    pipiens    = pip_its2,
    restuans   = rest_its2,
    unassigned = ITS2 - pip_its2 - rest_its2
  ) %>%
  dplyr::select(Code, ITS2, pipiens, restuans, unassigned) %>%
  pivot_longer(cols = c(pipiens, restuans, unassigned),
               names_to = "species", values_to = "count") %>%
  mutate(prop = count / ITS2,
         dataset = "Connecticut")

# ADD CONNECTICUT ORDERING
ct_order <- tribble(
  ~Code, ~urbanization,
  "WP69", "Rural",
  "GL77", "Rural",
  "K27",  "Rural",
  "SA105","Rural",
  "BP39", "Urban",
  "HA78", "Urban",
  "NH54", "Urban",
  "WH104","Urban"
)

ct_plot <- ct_plot %>%
  left_join(ct_order, by = "Code") %>%
  mutate(
    urbanization = factor(urbanization, levels = c("Rural", "Urban")),
    Code = factor(Code, levels = ct_order$Code)
  )

# -------------------------
# Utah
# -------------------------
ut <- tribble(
  ~City, ~Code, ~N, ~cqm1, ~ace2, ~pip_ace2, ~hyb_ace2, ~other_ace2, ~pip_sl_M, ~tars_M,
  "Audubon",              "AUDU",   1,   1,   1,   1, 0, 0,   462, 17063,
  "Avocet",               "AVOC",   1,   1,   1,   1, 0, 0,  1645, 40675,
  "Bird Rec",             "BIRD",   5,   4,   4,   3, 0, 1,  1439, 49430,
  "Blue Lake",            "BLUE",   5,   4,   4,   3, 0, 0,   368, 20997,
  "Conduit",              "COND",  13,  13,  13,  12, 1, 0,   815, 47352,
  "Downington",           "DOWA",  96,  95,  90,  89, 1, 0,   289, 14699,
  "Elephant Island",      "ELIS",   4,   4,   4,   3, 0, 1,  1458, 51701,
  "Fair Park",            "FAIR",  17,  17,  17,  17, 0, 0,   627, 48786,
  "Frontage Corral",      "FRCO",   7,   6,   6,   5, 0, 1,   325, 14443,
  "F.S.#13 Parley Way",   "FS13",  48,  46,  43,  43, 0, 0,   368, 12756,
  "F.S.#2 West High",     "FS2W",  34,  33,  33,  33, 0, 0,   139,   725,
  "F.S.#4 Avenues",       "FS4A",  20,  20,  19,  19, 0, 0,   768, 11848,
  "F.S.#5 Harvey Milk",   "FS5H",  52,  45,  46,  44, 2, 0,  2559, 13305,
  "F.S.#6 Indiana",       "FS6I",  57,  55,  55,  55, 0, 0,  2780,  2644,
  "F.S.#8 Ballpark",      "FS8B",  35,  31,  28,  27, 1, 0,  2914,  4035,
  "F.S.#9 Amelia Earhart","FS9A",  43,  42,  41,  41, 0, 0,  1573, 21291,
  "Glendale",             "GLEN",  80,  76,  74,  74, 0, 0,  5910, 13971,
  "Graystone Church",     "GRAYS",  1,   1,   1,   1, 0, 0,   258,  9847,
  "Hogle Zoo - Primates", "HOGL",  31,  30,  24,  24, 0, 0,  8842, 31066,
  "Jade",                 "JADE",   8,   7,   8,   7, 0, 1, 10940, 25024,
  "KETO",                 "KETO",  64,  64,  64,  64, 0, 0,   249,     5,
  "Lakefront",            "LAKE",   1,   1,   1,   1, 0, 0,  1986,  1280,
  "Lee Kay",              "LEEK",  35,  32,  32,  31, 0, 1,   481,    64,
  "New State",            "NEWS",  12,  12,  12,  12, 0, 0,   368,   188,
  "Nibley Golf Course",   "NIBL",  36,  35,  33,  33, 0, 0,   361,   746,
  "O'Reillys",            "OREI",   3,   3,   3,   3, 0, 0,   514,   281,
  "Power Pole",           "POPO",   5,   4,   5,   4, 0, 1,   861,  1062,
  "SLCMAD",               "SLCM",  14,  14,  14,  14, 0, 0,  1212,   563,
  "South Saratoga",       "SOSA",   7,   7,   6,   6, 0, 0,   323,     0,
  "Train Tracks",         "TRTR", 110, 103, 105, 101, 0, 4,   377,   219,
  "Wingpointe",           "WING", 119, 116, 115, 114, 0, 1,   642,   928
)

ut_plot <- ut %>%
  mutate(
    pipiens    = pip_ace2,
    hybrid     = hyb_ace2,
    other      = other_ace2,
    unassigned = ace2 - pip_ace2 - hyb_ace2 - other_ace2
  ) %>%
  dplyr::select(Code, ace2, pipiens, hybrid, other, unassigned) %>%
  pivot_longer(cols = c(pipiens, hybrid, other, unassigned),
               names_to = "species", values_to = "count") %>%
  mutate(prop = count / ace2,
         dataset = "Utah")

# ADD UTAH ORDERING
ut_order <- tribble(
  ~Code, ~urbanization,
  "AUDU","Rural",
  "AVOC","Rural",
  "BIRD","Rural",
  "BLUE","Rural",
  "ELIS","Rural",
  "JADE","Rural",
  "LAKE","Rural",
  "NEWS","Rural",
  "POPO","Rural",
  "SOSA","Rural",
  "FAIR","Peri",
  "FRCO","Peri",
  "FS9A","Peri",
  "GLEN","Peri",
  "KETO","Peri",
  "LEEK","Peri",
  "OREI","Peri",
  "SLCM","Peri",
  "TRTR","Peri",
  "WING","Peri",
  "COND","Urban",
  "DOWA","Urban",
  "FS13","Urban",
  "FS2W","Urban",
  "FS4A","Urban",
  "FS5H","Urban",
  "FS6I","Urban",
  "FS8B","Urban",
  "GRAYS","Urban",
  "HOGL","Urban",
  "NIBL","Urban"
)

ut_plot <- ut_plot %>%
  left_join(ut_order, by = "Code") %>%
  mutate(
    urbanization = factor(urbanization, levels = c("Rural", "Peri", "Urban")),
    Code = factor(Code, levels = ut_order$Code)
  )


# -------------------------
# combine
# -------------------------
plot_df <- bind_rows(
  ct_plot %>% rename(genotyped_n = ITS2),
  ut_plot %>% rename(genotyped_n = ace2)
)

plot_df$species <- factor(
  plot_df$species,
  levels = c("pipiens", "restuans", "hybrid", "other", "unassigned")
)

ggplot(plot_df, aes(x = Code, y = prop, fill = species)) +
  geom_col(width = 0.8) +
  facet_wrap(~dataset, scales = "free_x", ncol = 1) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(
    values = c(
      pipiens = "#00E5FF",
      restuans = "#FDAE61",
      hybrid = "#7B3294",
      other = "#FF00AA",
      unassigned = "grey80"
    ),
    na.value = "grey80"
  ) +
  
  labs(
    x = "Site",
    y = "Proportion of genotyped individuals",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
)
```

![](../figures/unnamed-chunk-1-1.png)<!-- -->

``` r
ct_long <- ct %>%
  left_join(ct_order, by = "Code") %>%
  mutate(
    pipiens = pip_its2,
    restuans = rest_its2
  ) %>%
  pivot_longer(cols = c(pipiens, restuans),
               names_to = "species", values_to = "count") %>%
  uncount(count) %>%
  mutate(species_bin = ifelse(species == "pipiens", 1, 0))


# simple glm
glm_fit <- glm(species_bin ~ urbanization, data = ct_long, family = binomial)

summary(glm_fit)
```

    ## 
    ## Call:
    ## glm(formula = species_bin ~ urbanization, family = binomial, 
    ##     data = ct_long)
    ## 
    ## Coefficients:
    ##                   Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)         0.7256     0.1183   6.133 8.64e-10 ***
    ## urbanizationUrban   2.0671     0.2013  10.267  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 837.04  on 1017  degrees of freedom
    ## Residual deviance: 716.24  on 1016  degrees of freedom
    ## AIC: 720.24
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
# glm with site random effect
library(lme4)
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
glm_fit <- glmer(
  species_bin ~ urbanization + (1 | Code),
  data = ct_long,
  family = binomial
)

summary(glm_fit)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: binomial  ( logit )
    ## Formula: species_bin ~ urbanization + (1 | Code)
    ##    Data: ct_long
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     595.0     609.8    -294.5     589.0      1015 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -8.4347  0.1186  0.2127  0.3895  3.6679 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Code   (Intercept) 2.975    1.725   
    ## Number of obs: 1018, groups:  Code, 8
    ## 
    ## Fixed effects:
    ##                   Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)         0.1894     0.8927   0.212   0.8320  
    ## urbanizationUrban   3.1798     1.2675   2.509   0.0121 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr)
    ## urbnztnUrbn -0.706

``` r
ct_habitat_counts <- ct_plot %>%
  filter(species %in% c("pipiens", "restuans")) %>%
  group_by(urbanization, species) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(
    species = factor(species, levels = c("pipiens", "restuans")),
    urbanization = factor(urbanization, levels = c("Rural", "Urban"))
  )

ggplot(ct_habitat_counts, aes(x = urbanization, y = count, fill = species)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  scale_fill_manual(
    values = c(
      pipiens = "#00E5FF",
      restuans = "#FDAE61"
    ),
    labels = c(
      pipiens = "Culex pipiens",
      restuans = "Culex restuans"
    )
  ) +
  labs(
    x = "Habitat",
    y = "Abundance (count)",
    fill = NULL
  ) +
  theme_classic() +
  theme(legend.position = "right")
```

![](../figures/plot-abundance-1.png)<!-- --> Predicted response is:  
- P(species_bin = 1) = predicted relative abundance of pipiens  
- 1 - P(species_bin = 1) = predicted relative abundance of restuans

To get the overall effect of urbanization, you want to exclude the
random site effect, so use re.form = NA.

``` r
library(dplyr)
library(tidyr)
library(ggplot2)

ct_long$species_bin
```

    ##    [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##   [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##   [75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [112] 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [149] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [186] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [223] 1 1 1 1 1 1 1 1 1 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [260] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [297] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [334] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [371] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [408] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [445] 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ##  [482] 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [519] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [556] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [593] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [630] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [667] 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ##  [704] 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [741] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [778] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [815] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [852] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ##  [889] 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [926] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ##  [963] 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
    ## [1000] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

``` r
ct_long$species
```

    ##    [1] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##    [7] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [13] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [19] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [25] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [31] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [37] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [43] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [49] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [55] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [61] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [67] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [73] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [79] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [85] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [91] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##   [97] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [103] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [109] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [115] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "restuans"
    ##  [121] "restuans" "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [127] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [133] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [139] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [145] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [151] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [157] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [163] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [169] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [175] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [181] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [187] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [193] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [199] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [205] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [211] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [217] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [223] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [229] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "restuans" "restuans"
    ##  [235] "restuans" "restuans" "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [241] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [247] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [253] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [259] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [265] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [271] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [277] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [283] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [289] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [295] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [301] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [307] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [313] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [319] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [325] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [331] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [337] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [343] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [349] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [355] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [361] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [367] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [373] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [379] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [385] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [391] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [397] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [403] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [409] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [415] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [421] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [427] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [433] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [439] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [445] "pipiens"  "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [451] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [457] "pipiens"  "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [463] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [469] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [475] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [481] "restuans" "restuans" "restuans" "restuans" "pipiens"  "pipiens" 
    ##  [487] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [493] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [499] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [505] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [511] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [517] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [523] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [529] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [535] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [541] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [547] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [553] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [559] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [565] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [571] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [577] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [583] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [589] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [595] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [601] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [607] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [613] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [619] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [625] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [631] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [637] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [643] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [649] "restuans" "restuans" "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [655] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [661] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [667] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [673] "pipiens"  "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [679] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [685] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [691] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [697] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [703] "restuans" "restuans" "restuans" "pipiens"  "pipiens"  "pipiens" 
    ##  [709] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [715] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [721] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [727] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [733] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [739] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [745] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [751] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [757] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [763] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [769] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [775] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [781] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [787] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [793] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [799] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [805] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [811] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [817] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [823] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [829] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [835] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [841] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [847] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [853] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [859] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [865] "pipiens"  "pipiens"  "restuans" "restuans" "restuans" "restuans"
    ##  [871] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [877] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [883] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [889] "restuans" "restuans" "restuans" "pipiens"  "pipiens"  "pipiens" 
    ##  [895] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [901] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [907] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [913] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [919] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [925] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [931] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [937] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [943] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [949] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [955] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [961] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [967] "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens"  "pipiens" 
    ##  [973] "pipiens"  "pipiens"  "pipiens"  "restuans" "restuans" "restuans"
    ##  [979] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [985] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [991] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ##  [997] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ## [1003] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ## [1009] "restuans" "restuans" "restuans" "restuans" "restuans" "restuans"
    ## [1015] "restuans" "restuans" "restuans" "restuans"

``` r
# prediction data
newdat <- data.frame(
  urbanization = factor(levels(ct_long$urbanization),
                        levels = levels(ct_long$urbanization))
)

# predict overall effect (exclude random site effect)
pred_link <- predict(
  glm_fit,
  newdata = newdat,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)
```

    ## Warning in predict.merMod(glm_fit, newdata = newdat, type = "link", se.fit =
    ## TRUE, : se.fit computation uses an approximation to estimate the sampling
    ## distribution of the parameters

``` r
# add fitted values and CIs
newdat <- newdat %>%
  mutate(
    fit_link = pred_link$fit,
    se_link  = pred_link$se.fit,
    pred_pipiens   = plogis(fit_link),
    lower_pipiens  = plogis(fit_link - 1.96 * se_link),
    upper_pipiens  = plogis(fit_link + 1.96 * se_link),
    pred_restuans  = 1 - pred_pipiens,
    lower_restuans = 1 - upper_pipiens,
    upper_restuans = 1 - lower_pipiens
  )

# reshape ONLY the prediction columns
plotdat <- newdat %>%
  dplyr::select(
    urbanization,
    pred_pipiens, lower_pipiens, upper_pipiens,
    pred_restuans, lower_restuans, upper_restuans
  ) %>%
  pivot_longer(
    cols = -urbanization,
    names_to = c(".value", "species"),
    names_pattern = "(pred|lower|upper)_(pipiens|restuans)"
  ) %>%
  mutate(
    species = factor(
      species,
      levels = c("pipiens", "restuans"),
    )
  )

cols <- c(
  "pipiens"  = "#1bc8ea",
  "restuans" = "#FDAE61"
)

ggplot(plotdat, aes(x = urbanization, y = pred, color = species, group = species)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08) +
  scale_color_manual(values = cols) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Urbanization",
    y = "Predicted relative abundance",
    color = NULL
  ) +
  theme_classic() +
  theme(legend.position = "right")
```

    ## Warning: No shared levels found between `names(values)` of the manual scale and the
    ## data's colour values.

![](../figures/plot-pred-1.png)<!-- -->

``` r
library(dplyr)
library(tidyr)
library(ggplot2)

# prediction data
newdat <- data.frame(
  urbanization = factor(levels(ct_long$urbanization),
                        levels = levels(ct_long$urbanization))
)

# predict overall effect (exclude random effect)
pred_link <- predict(
  glm_fit,
  newdata = newdat,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)
```

    ## Warning in predict.merMod(glm_fit, newdata = newdat, type = "link", se.fit =
    ## TRUE, : se.fit computation uses an approximation to estimate the sampling
    ## distribution of the parameters

``` r
# add fitted values and CIs
newdat <- newdat %>%
  mutate(
    fit_link = pred_link$fit,
    se_link  = pred_link$se.fit,
    pred_pipiens   = plogis(fit_link),
    lower_pipiens  = plogis(fit_link - 1.96 * se_link),
    upper_pipiens  = plogis(fit_link + 1.96 * se_link),
    pred_restuans  = 1 - pred_pipiens,
    lower_restuans = 1 - upper_pipiens,
    upper_restuans = 1 - lower_pipiens
  )

# reshape
plotdat <- newdat %>%
  pivot_longer(
    cols = -urbanization,
    names_to = c(".value", "species"),
    names_pattern = "(pred|lower|upper)_(pipiens|restuans)"
  ) %>%
  mutate(
    species = factor(species,
      levels = c("pipiens", "restuans"),
      labels = c("Culex pipiens", "Culex restuans")
    )
  )

# colors
cols <- c(
  "Culex pipiens"  = "#1bc8ea",
  "Culex restuans" = "#FDAE61"
)

# quick check
print(unique(plotdat$species))
```

    ## factor()
    ## Levels: Culex pipiens Culex restuans

``` r
# plot
ggplot(plotdat, aes(x = urbanization, y = pred, color = species, group = species)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08) +
  scale_color_manual(values = cols) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Urbanization",
    y = "Predicted relative abundance",
    color = NULL
  ) +
  theme_classic() +
  theme(legend.position = "right")
```

    ## Warning: No shared levels found between `names(values)` of the manual scale and the
    ## data's colour values.

![](../figures/unnamed-chunk-3-1.png)<!-- -->

``` r
unique(plotdat$species)
```

    ## factor()
    ## Levels: Culex pipiens Culex restuans

``` r
library(dplyr)
library(tidyr)
library(ggplot2)

# prediction data: one row per urbanization level
newdat <- expand.grid(
  urbanization = levels(ct_long$urbanization)
)

newdat$urbanization <- factor(newdat$urbanization, levels = levels(ct_long$urbanization))

# predict on link scale, excluding random effect
pred_link <- predict(glm_fit, newdata = newdat, type = "link", se.fit = TRUE, re.form = NA)
```

    ## Warning in predict.merMod(glm_fit, newdata = newdat, type = "link", se.fit =
    ## TRUE, : se.fit computation uses an approximation to estimate the sampling
    ## distribution of the parameters

``` r
# convert to probabilities
newdat <- newdat %>%
  mutate(
    fit_link = pred_link$fit,
    se_link  = pred_link$se.fit,
    pred_pipiens = plogis(fit_link),
    lower_pipiens = plogis(fit_link - 1.96 * se_link),
    upper_pipiens = plogis(fit_link + 1.96 * se_link),
    pred_restuans = 1 - pred_pipiens,
    lower_restuans = 1 - upper_pipiens,
    upper_restuans = 1 - lower_pipiens
  )

# reshape for plotting
plotdat <- newdat %>%
  pivot_longer(
    cols = c(pred_pipiens, lower_pipiens, upper_pipiens,
             pred_restuans, lower_restuans, upper_restuans),
    names_to = c(".value", "species"),
    names_pattern = "(pred|lower|upper)_(pipiens|restuans)"
  ) %>%
  mutate(
    species = recode(species,
      pipiens = "Culex pipiens",
      restuans = "Culex restuans"
    )
  )

# colors
cols <- c(
  "pipiens"  = "#1bc8ea",
  "restuans" = "#FDAE61"
)

# plot
ggplot(plotdat, aes(x = urbanization, y = pred, color = species, group = species)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.08) +
  scale_color_manual(values = cols) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Urbanization",
    y = "Predicted relative abundance",
    color = NULL
  ) +
  theme_classic() +
  theme(legend.position = "right")
```

    ## Warning: No shared levels found between `names(values)` of the manual scale and the
    ## data's colour values.

![](../figures/unnamed-chunk-4-1.png)<!-- -->

``` r
unique(ct_long$species)
```

    ## [1] "pipiens"  "restuans"

``` r
library(dplyr)
library(tidyr)
library(ggplot2)

# prediction grid
newdat <- data.frame(
  urbanization = factor(
    c("Rural", "Urban"),
    levels = c("Rural", "Urban")
  )
)

# predict from mixed model
pred_link <- predict(
  glm_fit,
  newdata = newdat,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)
```

    ## Warning in predict.merMod(glm_fit, newdata = newdat, type = "link", se.fit =
    ## TRUE, : se.fit computation uses an approximation to estimate the sampling
    ## distribution of the parameters

``` r
# convert to probabilities + confidence intervals
newdat <- newdat %>%
  mutate(
    fit_link = pred_link$fit,
    se_link  = pred_link$se.fit,

    pred_pipiens  = plogis(fit_link),
    lower_pipiens = plogis(fit_link - 1.96 * se_link),
    upper_pipiens = plogis(fit_link + 1.96 * se_link),

    pred_restuans  = 1 - pred_pipiens,
    lower_restuans = 1 - upper_pipiens,
    upper_restuans = 1 - lower_pipiens
  )

# reshape for plotting
plotdat <- newdat %>%
  dplyr::select(
    urbanization,
    pred_pipiens, lower_pipiens, upper_pipiens,
    pred_restuans, lower_restuans, upper_restuans
  ) %>%
  pivot_longer(
    cols = -urbanization,
    names_to = c(".value", "species"),
    names_pattern = "(pred|lower|upper)_(pipiens|restuans)"
  ) %>%
  mutate(
    species = recode(
      species,
      pipiens = "Culex pipiens",
      restuans = "Culex restuans"
    )
  )

# colors
cols <- c(
  "Culex pipiens"  = "#1bc8ea",
  "Culex restuans" = "#FDAE61"
)

# plot
ggplot(plotdat,
       aes(x = urbanization,
           y = pred,
           color = species,
           group = species)) +

  geom_point(size = 4) +
  geom_line(linewidth = 1) +

  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.08
  ) +

  scale_color_manual(values = cols) +

  scale_y_continuous(
    limits = c(0, 1)
  ) +

  labs(
    x = "Urbanization",
    y = "Predicted relative abundance",
    color = NULL
  ) +

  theme_classic() +
  theme(
    legend.position = "right"
  )
```

![](../figures/plot-predicted-1.png)<!-- -->

``` r
library(broom.mixed)
library(dplyr)
library(ggplot2)

coef_df <- tidy(glm_fit, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(
    odds_ratio = exp(estimate),
    lower = exp(estimate - 1.96 * std.error),
    upper = exp(estimate + 1.96 * std.error),
    term_clean = recode(term,
      "urbanizationUrban" = "Urban vs Rural"
    )
  )

ggplot(coef_df, aes(x = term_clean, y = odds_ratio)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(trans = "log10") +
  labs(
    x = NULL,
    y = "Odds ratio (pipiens vs restuans)"
  ) +
  theme_classic()
```

![](../figures/plot-odds-1.png)<!-- -->

# Map for UT replicating the one for CT

Figure 2. Map of Utah 2025 collections. Red: urban, yellow: periurban,
green: rural.

``` r
library(tidyverse)
library(maps)
```

    ## 
    ## Attaching package: 'maps'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     map

``` r
library(ggrepel)

# -------------------------
# site data
# -------------------------
ut_sites <- tribble(
  ~Code,  ~urbanization, ~lat,        ~lon,
  "AUDU", "Rural",       40.8261295, -112.11365,
  "AVOC", "Rural",       40.7948361, -112.12243,
  "BIRD", "Rural",       40.8021326, -112.05023,
  "BLUE", "Rural",       40.8299255, -111.98620,
  "ELIS", "Rural",       40.8520373, -112.04947,
  "JADE", "Rural",       40.8237778, -112.06829,
  "LAKE", "Rural",       40.8699330, -112.02340,
  "NEWS", "Rural",       40.8711671, -111.98135,
  "POPO", "Rural",       40.8222004, -112.01486,
  "SOSA", "Rural",       40.7843680, -112.10903,
  "FAIR", "Peri",        40.7736869, -111.92002,
  "FRCO", "Peri",        40.7689356, -112.08336,
  "FS9A", "Peri",        40.7786825, -112.02853,
  "GLEN", "Peri",        40.7276405, -111.93180,
  "KETO", "Peri",        40.7480485, -111.96788,
  "LEEK", "Peri",        40.7278171, -112.04577,
  "OREI", "Peri",        40.7561468, -112.03857,
  "SLCM", "Peri",        40.8160060, -111.95438,
  "TRTR", "Peri",        40.7542365, -112.00355,
  "WING", "Peri",        40.7715194, -111.97745,
  "COND", "Urban",       40.7528000, -111.89830,
  "DOWA", "Urban",       40.7308000, -111.85720,
  "FS13", "Urban",       40.7236000, -111.82260,
  "FS2W", "Urban",       40.7762000, -111.89890,
  "FS4A", "Urban",       40.7815000, -111.86470,
  "FS5H", "Urban",       40.7500000, -111.86130,
  "FS6I", "Urban",       40.7523000, -111.91830,
  "FS8B", "Urban",       40.7411000, -111.89160,
  "GRAYS","Urban",       40.7108000, -111.85920,
  "HOGL", "Urban",       40.7501000, -111.81610,
  "NIBL", "Urban",       40.7103000, -111.87270
) %>%
  mutate(
    urbanization = factor(urbanization, levels = c("Urban", "Peri", "Rural"))
  )

# -------------------------
# Utah county map
# -------------------------
utah_map <- map_data("county") %>%
  filter(region == "utah", subregion %in% c(
    "salt lake", "utah", "davis", "weber", "tooele", "summit"
  ))

# -------------------------
# plot
# -------------------------
ggplot() +
  geom_polygon(
    data = utah_map,
    aes(x = long, y = lat, group = group),
    fill = "grey80",
    color = "white",
    linewidth = 0.3
  ) +
  geom_point(
    data = ut_sites,
    aes(x = lon, y = lat, fill = urbanization),
    shape = 21,
    size = 4.5,
    color = "grey25",
    stroke = 0.7
  ) +
  geom_text_repel(
    data = ut_sites,
    aes(x = lon, y = lat, label = Code),
    size = 7,
    family = "",
    box.padding = 0.25,
    point.padding = 0.2,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  scale_fill_manual(
    values = c(
      "Urban" = "red",
      "Peri"  = "#F4B400",
      "Rural" = "green3"
    )
  ) +
  coord_quickmap(
    xlim = c(-112.16, -111.78),
    ylim = c(40.69, 40.89),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )
```

![](../figures/map-1.png)<!-- -->

``` r
library(tidyverse)
library(ggmap)
```

    ## ℹ Google's Terms of Service: <https://mapsplatform.google.com>
    ##   Stadia Maps' Terms of Service: <https://stadiamaps.com/terms-of-service>
    ##   OpenStreetMap's Tile Usage Policy: <https://operations.osmfoundation.org/policies/tiles>
    ## ℹ Please cite ggmap if you use it! Use `citation("ggmap")` for details.

``` r
library(ggrepel)


# 1. Register your Google Maps API key
register_google(key = "AIzaSyCdMfZPEp_xxkhdQrN2gfcGdBNpl1ga7ng")

# 2. Build a bounding box
lon_min <- -112.8
lon_max <- -111.2
lat_min <- 40.6
lat_max <- 41

# 3. Download the satellite map using the bounding box
bg <- get_map(
  location = c(
    left   = lon_min,
    bottom = lat_min,
    right  = lon_max,
    top    = lat_max
  ),
  maptype = "terrain",
  zoom = 11
)
```

    ## ! Bounding box given to Google - spatial extent only approximate.
    ## ℹ <https://maps.googleapis.com/maps/api/staticmap?center=40.8,-112&zoom=11&size=640x640&scale=2&maptype=terrain&language=en-EN&key=xxx>

``` r
# -------------------------
# Plot
# -------------------------
ggmap(bg) +
  geom_point(
    data = ut_sites,
    aes(x = lon, y = lat, fill = urbanization),
    shape = 21,
    size = 2.5,
    color = "black",
    stroke = 0.7
  ) +
  geom_text_repel(
    data = ut_sites,
    aes(x = lon, y = lat, label = Code),
    size = 1.5,
    box.padding = 0.03,
    point.padding = 0.02,
    segment.color = NA,
    max.overlaps = Inf
  ) +
  scale_fill_manual(
    values = c(
      "Urban" = "red",
      "Peri"  = "#F4B400",
      "Rural" = "green3"
    )
  ) +
  coord_quickmap(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max)) +
  theme_void() +
  theme(legend.position = "none")
```

    ## Coordinate system already present. Adding new coordinate system, which will
    ## replace the existing one.

![](../figures/google-map-1.png)<!-- -->
