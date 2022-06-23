
##----------------- ZHANG BRYANT MODEL ---------------------
##
##  Estimate Zhang and Bryant Bayesian small area Demography 
##  model on the first year of the simulated data
## 
##  Author: Benjamin Schlüter
##  Date: April 2022
##
##
##  Notes
##
## 
##-------------------------------------------------------

rm(list = ls())

## Load packages ----------------------------------------

packages <- c("tidyverse", "ggplot2", "HMDHFDplus", "viridis", "scales",
              "forecast", "cowplot", "RColorBrewer", "rstan", "tidybayes")
invisible( lapply(packages, library, character.only = TRUE))


## HMD identifiers --------------------------------------

userID <- "benjamin-samuel.schluter@uclouvain.be"
pwd <- "Marieke0917"


## Load data --------------------------------------------

# Simulated data
df <- readRDS("../../../PhD/Research Projects/subnat models comparison/data/simulation/df_simulated.rda")

# Age composition SVK 2019
age.comp <- readHMDweb("SVK", "Exposures_5x1",
                      userID, pwd) %>%
        filter(Year == 2019) %>% 
        dplyr::select(Year, Age, Male) %>% 
        rename("Pop" = Male) %>% 
        mutate(Age = ifelse(Age >100, 100, Age)) %>% 
        # Sum age > 100
        group_by(Year, Age) %>% 
        summarise(Pop = sum(Pop)) %>% 
        ungroup() %>%  
        group_by(Year) %>% 
        mutate(sumExp = sum(Pop),
               prop = Pop/sumExp) %>% 
        pull(prop)


## Simulate death counts for a given N ------------------

# Simulated nmx
nmx <- df %>% 
        filter(unit == "SVK-1",
               year == 1) %>% 
        pull(nmx)

# Generate death counts
set.seed(02052022)
N <- 1000000
dth <- rpois(length(nmx), nmx*(N*age.comp))

# Store in df
df.dth <- tibble(age = c(0, 1, seq(5, 100, 5)),
                 dth = dth,
                 pop = N*age.comp) %>% 
        mutate(nmx = ifelse(dth == 0, NA, dth/pop)) 

ggplot(df.dth, aes(x = age, y = nmx)) +
        geom_point() +
        theme_bw() +
        scale_y_log10()


## Estimate simplified Bryant & Zhang model for age schedule -------

stan_data <- list(A = dim(df.dth)[1],
                  n = df.dth$pop,
                  d = df.dth$dth,
                  age0 = ifelse(df.dth$age == 0, 1, 0))


options(mc.cores = parallel::detectCores()-1)

# Pars HMC
niter <- 4000
nwarmup <- niter/2
nchains <- 4
nsim <- (niter-nwarmup)*nchains

fit = stan("./code/stan/bryant_zhang_model_t_dim.stan",
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data,
           control=list(adapt_delta  = 0.90, 
                        max_treedepth= 15))


## Models' diagnostics --------------------------------------

traceplot(fit, pars = c("beta_a[1:10]"),
          ncol = 2)

## Models' outputs --------------------------------------


# With tidyvayes ?

nmx.hat <- exp(extract(fit)$eta)


df.fit <- apply(nmx.hat, 2, quantile, probs = c(0.025, 0.5, 0.975)) %>% 
        t() %>% 
        as.data.frame() %>% 
        mutate(age = c(0, 1, seq(5, 85, 5)),
               nmx.true = nmx[1:19]) %>% 
        rename("nmx.l" = "2.5%",
               "nmx.m" = "50%",
               "nmx.u" = "97.5%")


ggplot() +
        geom_point(data = df.dth, aes(x = age.gp, y = nmx),
                   size = 2) +
        geom_ribbon(data = df.fit, aes(x = age, ymin = nmx.l, ymax = nmx.u),
                    alpha = 0.3,
                    col = "skyblue3",
                    fill = "skyblue3") +
        geom_line(data = df.fit, aes(x = age, y = nmx.m),
                  col = "white") +
        geom_line(data = df.fit, aes(x = age, y = nmx.true),
                  linetype = "dashed") +
        theme_bw() +
        scale_y_log10()



## Models' performance ------------------------------------

# RMSE
df.fit <- df.fit %>% 
        mutate(sq.diff = (nmx.m - nmx.true)^2)
rmse <- sqrt( (1/dim(df.fit)[1])*sum(df.fit$sq.diff) )

# Coverage
df.fit <- df.fit %>% 
        mutate(hthanl = nmx.true >= nmx.l,
               lthanu = nmx.true < nmx.u,
               inCI = hthanl * lthanu)
cov <- ( (1/dim(df.fit)[1])*sum(df.fit$inCI) ) * 100 






## ---------------------------------------------------------
## Model does not offer great flexibility at young ages
## hence test with 5y age groups, 2 sexes and 10 years
## for France and reducing its exposure
## ---------------------------------------------------------


period <- 2010:2019

## Data for FRA --------------------------------------------

# Age composition by 5y age groups
age.comp <- readHMDweb("FRATNP", "Exposures_5x1",
                       userID, pwd) %>%
        filter(Year %in% period) %>% 
        dplyr::select(Year, Age, Male, Female) %>% 
        pivot_longer(Male:Female, names_to = "sex", values_to = "exp") %>% 
        group_by(Year, sex) %>% 
        mutate(sumExp = sum(exp),
               prop = exp/sumExp) %>% 
        arrange(sex, Year, Age) %>% 
        pull(prop)

# Obtain 5mx
nmx_f <- readHMDweb("FRATNP", "fltper_5x1",
                    userID, pwd) %>%
        filter(Year %in% period) %>% 
        dplyr::select(Year, Age, mx) %>% 
        pull(mx)
nmx_m <- readHMDweb("FRATNP", "mltper_5x1",
                  userID, pwd) %>%
        filter(Year %in% period) %>% 
        dplyr::select(Year, Age, mx) %>% 
        pull(mx)

nmx <- c(nmx_f, nmx_m)


## Simulate death counts for a given N ------------------

set.seed(02052022)
# Set pop size to Maori in Zhang Bryant Nissen (2019)
N <- 328000
dth <- rpois(length(nmx), nmx*(N*age.comp))

# Store in df
age.gp <- c(0, 1, seq(5, 110, 5))
sex = c("f", "m")
df.sim <- expand.grid(age = age.gp, year = period, sex = sex) %>% 
        mutate(d = dth,
               e = N*age.comp,
               mx_sim = ifelse(d == 0, NA, d/e),
               mx_true = nmx) %>% 
        # !! order matters for stan modeling !!
        # --> storing in arrays
        arrange(year, age, sex)
# Check simulated outputs
df.sim %>% 
        pivot_longer(mx_sim:mx_true, names_to = "type", values_to = "nmx") %>% 
ggplot(aes(x = age, y = nmx, col = sex)) +
        geom_point() +
        theme_bw() +
        scale_y_log10() +
        facet_wrap(~ type)
# Reducing exposure increased data variability
# as expected and create some zero counts (NAs)

## Estimate full Bryant & Zhang model ------------------------

stan_data <- list(A = length(age.gp),
                  T = length(period),  
                  S = 2,
                  E = array(df.sim$e, 
                            dim = c(2, length(age.gp), length(period)) ),
                  D = array(df.sim$d, 
                            dim = c(2, length(age.gp), length(period)) ),
                  age0 = ifelse(age.gp == 0, 1, 0),
                  sex = c(0, 1)
                  )


options(mc.cores = parallel::detectCores()-1)


fit = stan("./code/stan/bryant_zhang_age_schedules_ast_dims.stan",
           iter = 6000,
           chains = 4,
           data = stan_data)


## Models' outputs --------------------------------------

# Define factors in df.sim so that
# they are recognized by recover_types
df.sim <- df.sim %>% 
        mutate(year = factor(year),
               age = factor(age))

# Obtain posterior median and 95% CI
# by strata
df.fit <- fit %>%
        # recover the factors from data
        recover_types(df.sim) %>%
        spread_draws(eta[sex, age, year]) %>% 
        median_qi()
# Add true mx
df.fit %>% 
        left_join(df.sim %>% 
                          select(age, year, sex, mx_true, mx_sim),
                  by = c("age", "year", "sex")) %>% 
        mutate(ln_mx_true = log(mx_true),
               ln_mx_sim = log(mx_sim)) %>% 
        ggplot(aes(x = age, y = eta, col = year, group = year)) +
        geom_line() +
        # geom_point(aes(y = ln_mx_true)) +
        geom_point(aes(y = ln_mx_sim)) +
        theme_bw() +
        facet_wrap(~sex)
# Drop at the end because counts goes to zero. 
# For small population, go until age 100.

