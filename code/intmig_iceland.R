
##----------------- BDEF ---------------------
##
##  Reproduce examples of BDEF with STAN code
##
## Chapter 15: Internal migration in Iceland
## 
##  Author: Benjamin Schlüter
##  Date: July 2022
##
##
##  Notes:
##
## 
##----------------------------------------------------------------

rm(list = ls())


## Load packages -------------------------------------------------

packages <- c("tidyverse", "ggplot2", "viridis", 
              "scales", "rstan", "tidybayes")
invisible( lapply(packages, library, character.only = TRUE))


## Functions ----------------------------------------------------

round3 <- function(x) {
        if (x%%3 == 0) {
                out <- x
        } else {
                sign <- sample(c(-1,1), 1, prob = c(0.5, 0.5))
                dev <- sample(c(1,2), 1, prob = c(2/3, 1/3))
                out <- x + sign*dev
        }
        return(out)
}


## Load and tidy data -------------------------------------------

mig <- bdefdata::iceland_migration
pop <- bdefdata::iceland_population

# Convert to data frame for EDA
df.mig <- as.data.frame.table(mig) %>% 
        rename("mig" = Freq) %>% 
        # apply random rounding base 3
        mutate(mig_c = map(mig, round3))

df.pop <- as.data.frame.table(pop) %>% 
        rename("pop" = Freq,
               "region_orig" = region) 
# Need to compute exposure



df <- df.mig %>% 
        left_join(df.pop,
                  by = c("age", "sex", "time", "region_orig"))
