

# Research Topic: testing whether habitat and seasonal partitioning between Culex pipiens s.l. and Culex tarsalis shapes West Nile Virus (WNV) dynamics across urban–rural gradients. Core hypothesis: early/mid-season amplification dominated by pipiens in urban areas, later spillover involving tarsalis moving into urban/peri-urban areas.

# Approach: Preliminary results visualized via mapping, with species identity and abundance as primary response variables. Model mosquito abundance and proportions using GLMMs count ~ season*urbanization + trap_type + (1|site/date), family = poisson(link = "log")

## Packages needed:

library(tidyverse) # for data wrangling
library(glmmTMB)   # for model fitting
library(DHARMa)    # for residual plots
library(emmeans)   # for estimating marginal effects
library(multcomp)  # for statistical comparisons on fitted models


## tarsalis datasets from SLCMAD:
tarsalis <- read.csv("/cloud/project/tarsalis_2025.csv")
head(tarsalis, 5)

## pipiens datasets from SLCMAD:
pipiens <- read.csv("/cloud/project/pipiens_2025.csv")
head(pipiens, 5)

# we'll model salamander count as a function of species (spp) and mining presence (mined).
# Because the data are grouped into sites, we'll add that as a random effect.
# first, look at the response distribution

# sub in tarsalis dataset here instead of salamander 
hist(tarsalis$count)
hist(log(tarsalis$count))

# tarsalis count data, starting with poisson model:
fit_pois <- glmmTMB(count ~ season*urban_cat + trap_type + (1 | site_name),
                    family = poisson(link = "log"),
                    data = tarsalis)

summary(fit_pois)
# look at the residuals:
simulateResiduals(fit_pois, plot = T)

#################################################
#continue salamander dataset example
# this doesn't look great. It seems like there is a trend in the residuals with
# greater error at higher values. Maybe the negative binomial distribution can help with this

fit_nbinom <- glmmTMB(count ~ spp*mined + (1 | site),
                      family = nbinom1(link = "log"),
                      data = Salamanders)
simulateResiduals(fit_nbinom, plot = T)

# this looks much better. Now we can look at the marginal effects:

# First, just look at the effect of species, doing a 1-way pairwise comparison
em_spp = emmeans(fit_nbinom, pairwise ~ spp, type = "response")
em_spp

cl = cld(em_spp, Letters = letters)

ggplot(data = cl, aes(x = spp, y = response)) +
  geom_point(size=2.5, color="black") +
  geom_errorbar(aes(x=spp, ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.2, size=1, color="black") +
  geom_text(aes(label = gsub(" ", "", .group)),
            position = position_nudge(x = 0.3)) +
  ggeasy::easy_remove_axes("both", "title")


# Next, let's see how the effects vary with different mined levels in a 2-way emmeans comparison
em_spp2 = emmeans(fit_nbinom, pairwise ~ spp | mined, type = "response")
em_spp2

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
