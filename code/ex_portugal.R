
##----------------- BDEF ---------------------
##
##  Reproduce examples of BDEF with STAN code
##
## Chapter 12: Life expectancy Portugal
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


## Load and tidy data --------------------------------------------

dth <- bdefdata::portugal_deaths 
exp <- bdefdata::portugal_exposure

# Convert to data frame for EDA
df.exp <- as.data.frame.table(exp) %>% 
        rename("exp" = Freq)

df.dth <- as.data.frame.table(dth) %>% 
        rename("dth" = Freq)

# Join df
df <- df.exp %>% 
        left_join(df.dth,
                  by = c("age", "sex", "year")) %>% 
        mutate(age = sub("\\-.*", "", as.character(age)),
               age = gsub('[+]', '', age),
               age = as.numeric(age),
               mx = dth/exp)



## EDA ------------------------------------------------------------

# Direct estimates
df %>% 
        filter(year %in% c("1950", "1970", "1990", "2010")) %>% 
        ggplot(aes(x = age, y = mx, group = sex, linetype = sex)) +
        facet_wrap(~year, ncol = 4) +
        geom_line() +
        theme_bw() +
        theme(legend.position = "top") +
        scale_y_log10() +
        scale_x_continuous(breaks = seq(0, 100, 10))


## STAN data ------------------------------------------------------

# Only use 1950-1995 to estimate model (heldback data)
df.stan <- df %>% 
        filter(year %in% as.character(1950:1995))

# Order D & E vector according to STAN code
# Female = 1st position
D <- df.stan %>%
        arrange(sex, year, age) %>% 
        pull(dth)
E <- df.stan %>% 
        arrange(sex, year, age) %>% 
        pull(exp)

stan_data <- list(n = dim(df.stan)[1],
                  A = length(unique(df.stan$age)),
                  T = length(unique(df.stan$year)),
                  S = length(unique(df.stan$sex)),
                  Tproj = 20,
                  D = D,
                  E = E,
                  age0 = c(1, rep(0, 21)),
                  sex = c(0, 1)
                  )

# saveRDS(stan_data,
#         "./server/in/stan_data.rda")


## Fit in STAN -------------------------------------------------

options(mc.cores = parallel::detectCores()-1)

# Pars HMC
niter <- 2000
nwarmup <- niter/2
nchains <- 4
nsim <- (niter-nwarmup)*nchains

model <- "./code/stan/ex_portugal.stan"

fit = stan(model,
           pars = c("eta_proj"),
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data, 
           control = list(adapt_delta  = 0.95, 
                          max_treedepth= 12)
)

# print(proc.time())

fit <- array(as.matrix(fit,'eta_proj'), c(nsim, 2, 22, 20),
             dimnames=list(paste0('sample',1:nsim),
                           1:2,
                           c(0, 1, seq(5, 100, 5)),
                           1996:2015))


## Forecasted mx -----------------------------------------------

# Load STAN output
# fit <- readRDS("./server/out/eta_proj.rda")


eta.proj <- apply(exp(fit), c(2, 3, 4), quantile, 
                  probs = c(0.1, 0.5, 0.9))
df.eta.proj <- as.data.frame.table(eta.proj) %>% 
        rename("sex" = Var2,
               "age" = Var3,
               "year" = Var4) %>% 
        pivot_wider(names_from = Var1, values_from = Freq)
df.eta.proj %>% 
        mutate(year = as.character(year) %>% as.numeric()) %>% 
        filter(year %in% sample(1996:2015, 6)) %>% 
        ggplot(aes(x = age, y = `50%`, group = sex, linetype = sex)) +
        facet_wrap(~ year) +
        geom_line() +
        theme_bw() +
        scale_y_log10()


## Forecasted ex -----------------------------------------------

ages <- c(0, 1, seq(5, 100, 5))
n_age <- length(ages)
dth.proj <- array(NA, dim = c(4000, 2, 22, 20))

# Exposure for forecast
E_proj <- df %>% 
        filter(year %in% as.character(1996:2015)) %>% 
        arrange(year, age, sex) %>% 
        pull(exp) %>% 
        array(dim = c(2, 22, 20))

# Generate deaths counts from Poisson dist. 
set.seed(08052022)
for (y in 1:20) {
        for (s in 1:2) {
                for (it in 1:4000) {
                        dth.proj[it, s, , y] <- rpois(n_age, exp(fit[it, s, , y]) * E_proj[s, , y])
                                      }
        }
}
# Compute finite-population mx
mx.proj <- array(NA, dim = c(4000, 2, 22, 20),
                 dimnames = list(1:4000,
                              1:2, 
                              ages,
                              1996:2015)) 
for (it in 1:4000) {
        mx.proj[it, , , ] <- dth.proj[it, , , ]/E_proj
}

# Compute life expectancy
source("./code/function/derive_ex.R")
ex.proj <- apply(mx.proj, c(1,2,4), get.ex, ages = ages)

ex.proj.ci <- apply(ex.proj, c(1,3,4), quantile, 
                    probs = c(0.1, 0.5, 0.9),
                    # one e100 is NA (0/0)
                    na.rm = TRUE)

# Store in a data frame
df.ex.proj <- as.data.frame.table(ex.proj.ci) %>% 
        rename("age" = Var2,
               "sex" = Var3,
               "year" = Var4) %>% 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        mutate(age = paste0("e", as.character(age)),
               sex = ifelse(sex == 1, "Female", "Male"))

# Compute ex from raw data
df.nested <- df %>% 
        select(-c(dth, exp)) %>% 
        group_by(sex, year) %>% 
        nest()

df.nested <- df.nested %>% 
        mutate(ex = map(data, ~{get.ex(.x$age, .x$mx)}),
               e0 = map_dbl(ex, ~{.x[1]}),
               e65 = map_dbl(ex, ~{.x[15]})) %>% 
        pivot_longer(e0:e65, names_to = "age", values_to = "le")


# Similar to Fig 12.10 from BDEF book 
ggplot() +
        facet_wrap(age ~ sex, scales = "free_y") +
        geom_line(data = df.nested,
                  aes(x = year, y = le, group = 1)) +
        geom_ribbon(data = df.ex.proj %>% filter(age %in% c("e0", "e65")),
                    aes(x = year, ymin = `10%`, ymax = `90%`, group = 1),
                    col = "blue", fill = "blue", alpha = 0.3) +
        geom_line(data = df.ex.proj %>% filter(age %in% c("e0", "e65")),
                  aes(x = year, y = `50%`, group = 1),
                  col = "white") +
        theme_bw()
        
