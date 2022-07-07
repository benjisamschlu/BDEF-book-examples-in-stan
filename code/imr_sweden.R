
##----------------- BDEF ---------------------
##
##  Reproduce examples of BDEF with STAN code
##
## Chapter 11: Infant mortality rates Sweden
## 
##  Author: Benjamin Schlüter
##  Date: June 2022
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


## Load and tidy data --------------------------------------------

dth <- bdefdata::sweden_deaths 
bth <- bdefdata::sweden_births

# Convert to data frame for EDA
df.bth <- bth %>% as_tibble() %>% 
        mutate(county = rownames(bth)) %>% 
        pivot_longer(!county, values_to = "bth", names_to = "year")

df.dth <- dth %>% as_tibble() %>% 
        mutate(county = rownames(dth)) %>% 
        pivot_longer(!county, values_to = "dth", names_to = "year")

# Join df
order.cty <- rownames(dth)
df <- df.bth %>% 
        left_join(df.dth,
                  by = c("county", "year")) %>% 
        mutate(year = as.numeric(year),
               imr = 1000*(dth/bth),
               # set levels order as data order
               county = factor(county,
                               levels = order.cty,
                               labels = order.cty)) %>%
        # used for order of county by nber of births
        group_by(county) %>% 
        mutate(tot_dth = sum(dth)) %>% 
        ungroup() 
ordered.cty <- rev(unique(df$county))



## EDA ------------------------------------------------------------

sel.cty <- c("Gotland", "Varmland", "Halland", "Vasterbotten", "Kalmar",
             "Dalarna", "Orebro", "Sodermanland", "Uppsala", "Stockholm")
# Death counts
df %>% 
        mutate(county = factor(county,
                               level = ordered.cty,
                               labels = ordered.cty)) %>% 
        filter(county %in% sel.cty) %>% 
        ggplot(aes(x = year, y = dth, group = 1)) +
        facet_wrap(~county, ncol = 5,
                   scale = "free_y") +
        geom_line() +
        theme_bw() +
        scale_x_continuous(breaks=seq(1995, 2015, 4))        

# Birth counts
df %>% 
        mutate(county = factor(county,
                               level = ordered.cty,
                               labels = ordered.cty),
               bth = bth/1000) %>% 
        filter(county %in% sel.cty) %>% 
        ggplot(aes(x = year, y = bth, group = 1)) +
        facet_wrap(~county, ncol = 5,
                   scale = "free_y") +
        geom_line() +
        theme_bw() +
        scale_x_continuous(breaks=seq(1995, 2015, 4))        

# Infant mortality rates
df %>% 
        mutate(county = factor(county,
                               level = ordered.cty,
                               labels = ordered.cty)) %>% 
        filter(county %in% sel.cty) %>% 
        ggplot(aes(x = year, y = imr, group = 1)) +
        facet_wrap(~county, ncol = 5,
                   scale = "free_y") +
        geom_line() +
        theme_bw() +
        scale_x_continuous(breaks=seq(1995, 2015, 4))    


## STAN data ------------------------------------------------------


stan_data <- list(T = length(unique(df$year)),
                  C = length(unique(df$county)),
                  y = dth,
                  w = bth,
                  nprojyears = 10)

options(mc.cores = parallel::detectCores()-1)

# Pars HMC
niter <- 3000
nwarmup <- niter/2
nchains <- 4
nsim <- (niter-nwarmup)*nchains

fit = stan("./code/stan/imr_sweden.stan",
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data, 
           control = list(adapt_delta  = 0.95, 
                          max_treedepth= 12)
           )
# Save/load data
# saveRDS(fit,
#         "./estimates/imr_sweden.rda")
# fit <- readRDS("./estimates/imr_sweden.rda")

## Model's diagnostics -------------------------------------------



## Model's outputs -----------------------------------------------

# Define factors in df.sim so that
# they are recognized by recover_types
df <- df %>% 
        mutate(year = factor(year))

# Obtain posterior median and 95% CI
# by strata
# Would be nice to have 95% AND 50% CI !!
df.fit <- fit %>%
        # recover the factors from data
        # Careful: county order in dth must be the same in df !!
        recover_types(df) %>%
        spread_draws(gamma[county, year]) %>% 
        median_qi()  

# Similar to Fig 10.11 from BDEF book (no 50% CI)
ggplot() +
        facet_wrap(~county, ncol = 5) +
        geom_line(data = df %>% 
                          mutate(year = as.numeric(as.character(year)),
                                 county = factor(county,
                                                 level = ordered.cty,
                                                 labels = ordered.cty)) %>% 
                          filter(county %in% sel.cty),
                  aes(x = year, y = imr, group = 1)) +
        geom_ribbon(data = df.fit %>% 
                          mutate(year = as.numeric(as.character(year)),
                                 imr = gamma*1000,
                                 county = factor(county,
                                                 level = ordered.cty,
                                                 labels = ordered.cty)) %>% 
                          filter(county %in% sel.cty),
                  aes(x = year, ymin = .lower*1000, ymax = .upper*1000), fill = "skyblue3",
                  col = "skyblue3", alpha = 0.25) +
        geom_line(data = df.fit %>% 
                          mutate(year = as.numeric(as.character(year)),
                                 imr = gamma*1000,
                                 county = factor(county,
                                                 level = ordered.cty,
                                                 labels = ordered.cty)) %>% 
                          filter(county %in% sel.cty),
                  aes(x = year, y = imr), col = "white") +
        theme_bw() +
        scale_x_continuous(breaks=seq(1995, 2015, 5)) +
        labs(y = "IMR (per 1000)")

# extract forecast from stan object
df.forecast <- fit %>%
        recover_types(df) %>% 
        spread_draws(gamma_proj[county, nprojyears]) %>% 
        median_qi() %>% 
        mutate(year = nprojyears + 2015)

# Similar to Fig 11.14 from BDEF book (no 50% CI)
ggplot() +
        facet_wrap(~county, ncol = 5) +
        geom_line(data = df %>% 
                          mutate(year = as.numeric(as.character(year)),
                                 county = factor(county,
                                                 level = ordered.cty,
                                                 labels = ordered.cty)) %>% 
                          filter(county %in% sel.cty),
                  aes(x = year, y = imr, group = 1)) +
        geom_ribbon(data = df.fit %>% 
                            mutate(year = as.numeric(as.character(year)),
                                   county = factor(county,
                                                   level = ordered.cty,
                                                   labels = ordered.cty)) %>% 
                            filter(county %in% sel.cty),
                    aes(x = year, ymin = .lower*1000, ymax = .upper*1000), fill = "skyblue3",
                    col = "skyblue3", alpha = 0.25) +
        geom_line(data = df.fit %>% 
                          mutate(year = as.numeric(as.character(year)),
                                 imr = gamma*1000,
                                 county = factor(county,
                                                 level = ordered.cty,
                                                 labels = ordered.cty)) %>% 
                          filter(county %in% sel.cty),
                  aes(x = year, y = imr), col = "white") +
        theme_bw() +
        # Add forecast to previous plot
        geom_ribbon(data = df.forecast %>% 
                            mutate(county = factor(county,
                                                   level = ordered.cty,
                                                   labels = ordered.cty)) %>% 
                            filter(county %in% sel.cty),
                    aes(x = year, ymin = .lower*1000, ymax = .upper*1000), fill = "red4",
                    col = "red4", alpha = 0.25) +
        geom_line(data = df.forecast %>% 
                          mutate(imr = gamma_proj*1000,
                                 county = factor(county,
                                                 level = ordered.cty,
                                                 labels = ordered.cty)) %>% 
                          filter(county %in% sel.cty),
                  aes(x = year, y = imr), col = "white", linetype = "dashed") +
        theme_bw() +
        scale_x_continuous(breaks=seq(1995, 2025, 10)) +
        labs(y = "IMR (per 1000)")





