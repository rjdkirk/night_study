# Final model and analysis pipeline.

library(tidyverse)
library(GGally)
library(ggpubr)
library(reshape2)
library(brms) 
library(marginaleffects) 
library(ggdist) 
library(cowplot) 


# Set up ------------------------------------------------------------------
setwd()

# Load data
data_sleep <-
  read.csv('data_sleep.csv')

# Previous data cleaning:
# Realigned night numbers where participants had skipped/missed nights (four occurrences)

# Preliminary tests -------------------------------------------------------

# Look at the distributions of our variables.
ggplot(melt(data_sleep[,c(5:7,28,9,10,15)], id = 'cndtn'), aes(value)) +
  geom_histogram() +
  facet_wrap(~variable, scales = 'free') +
  theme_minimal()


# Scatter plot matrix using numerical versions of each variable.
ggpairs(data_sleep, columns = c(6,7,28,9,10,15), aes(colour = cndtn, alpha = 0.7),
                    upper = list(continuous = wrap(ggally_cor, method = "kendall", size = 3)),
                    lower = list(continuous = wrap("points", alpha = 0.7, size = 0.3)),
                    progress = FALSE) +
  theme(panel.background = element_blank())


# SQ mdl ------------------------------------------------------------------

mdl_sq <-
  brm(
    sleep_qual ~
      cndtn * (alertness + tension + valence) +
      (cndtn | ID), # random slope of condition by participant for within comparison
    data = data_sleep,
    family = cumulative(link = "probit"), # family for ordinal variable
    prior = set_prior("normal(0, 1)", class = "b"), 
    control = list(adapt_delta = 0.99), # increased to avoid convergence warnings
    file = "mdl_sq"
    )

# Check posterior predictions
pp_check(mdl_sq, ndraws = 100)
ggsave('sq_pp.jpeg')

# Extract some effects of interest using marginaleffects
avg_predictions(mdl_sq, by = "cndtn", type = "link") 

# Plot on the link scale
plot_predictions(mdl_sq, by = "cndtn", type = "prediction") + 
  labs(x = 'Condition', y = 'SQ') +
  theme_classic()
ggsave('sq.jpeg')

# Look at pairwise comparisons
avg_comparisons(mdl_sq, variables = list(cndtn = "pairwise"), type = "link") |>
  print() |>
  posterior_draws(shape = "DxP") |>
  hypothesis(c("b1 > 0",
               "b2 > 0",
               "b3 > 0",
               "b4 < 0",
               "b5 < 0",
               "b6 > 0"))


# SOL mdl -----------------------------------------------------------------

# We've chosen a gamma distribution for this model. However, there is one 0 value in our data set, which causes
# an issue because the gamma distribution doesn't work on zero values. We'll add a +1 constant to the variable
# to get around this. 

data_sleep$sleep_onset_add1 <- data_sleep$sleep_onset + 1

mdl_sol_gamma <-
  brm(
    sleep_onset_add1 ~
      cndtn * (alertness + tension + valence) +
      (cndtn | ID),
    data = data_sleep,
    family = "gamma", 
    prior = set_prior("normal(0, 1)", class = "b"),
    control = list(adapt_delta = 0.99), 
    file = "mdl_sol_gamma"
    )

# Check posterior predictions
pp_check(mdl_sol_gamma, ndraws = 100)
ggsave('sol_pp.jpeg')

# Extract some effects of interest using marginaleffects
avg_predictions(mdl_sol_gamma, by = "cndtn", type = "link") 

# Plot on the link scale
plot_predictions(mdl_sol_gamma, by = "cndtn", type = "prediction") + 
  labs(x = 'Condition', y = 'SOL') +
  theme_classic() 
ggsave('sol.jpeg')

# Look at pairwise comparisons
avg_comparisons(mdl_sol_gamma, variables = list(cndtn = "pairwise"), type = "link") |>
  print() |>
  posterior_draws(shape = "DxP") |>
  hypothesis(c("b1 > 0",
               "b2 > 0",
               "b3 > 0",
               "b4 < 0",
               "b5 < 0",
               "b6 > 0"))

# SE mdl ------------------------------------------------------------------

# We will use a gamma distribution for this model as well. Again, we need to make some adjustments to the 
# variable in order to proceed; in order for it to work we need to invert the variable so that it is right
# rather than left skewed. This isn't too bad as there are hard limits to the scale, 0-100, so we can 
# reasonably flip the values and the outcomes will still make some intuitive sense. Likewise, we can add a 
# truncation to model to avoid a spuriously long tail.

# There is also a single data point with SE > 100, which is not possible, so must be an anomaly and we'll
# remove it.

data_sleep$sleep_effncy_inv <- 100 - data_sleep$sleep_effncy

mdl_se_gamma_trunc <-
  brm(
    sleep_effncy_inv | trunc(ub = 100) ~
      cndtn * (alertness + tension + valence) +
      (cndtn | ID),
    data = data_sleep,
    family = "gamma", 
    prior = set_prior("normal(0, 1)", class = "b"),
    control = list(adapt_delta = 0.99), 
    file = "mdl_se"
    )

# Check posterior predictions
pp_check(mdl_se_gamma_trunc, ndraws = 100)
ggsave('se_pp.jpeg')

# Extract some effects of interest using marginaleffects
avg_predictions(mdl_se_gamma_trunc, by = "cndtn", type = "link") 

# Plot on the link scale
plot_predictions(mdl_se_gamma_trunc, by = "cndtn", type = "prediction") + 
  labs(x = 'Condition', y = 'SE') +
  theme_classic() 
ggsave('se.jpeg')

# Look at pairwise comparisons
avg_comparisons(mdl_se_n, variables = list(cndtn = "pairwise"), type = "link") |>
  print() |>
  posterior_draws(shape = "DxP") |>
  hypothesis(c("b1 > 0",
               "b2 > 0",
               "b3 > 0",
               "b4 < 0",
               "b5 < 0",
               "b6 > 0"))


# Preference comparisons --------------------------------------------------

# Let's look at the balance of preference. The main df pref column is repeated for each participant across all 
# their nights, so we need to condense this down to the unique response per participant. 
pref_playlist_df <- unique(data_sleep[,c('ID','pref_playlist')])

ggplot(pref_playlist_df, aes(pref_playlist)) +
  geom_histogram(stat = 'count', alpha = 0.7) +
  theme_minimal()

ggsave('pref_playlist.jpeg')

# Now we'll repeat the models but with music preference as our grouping variable instead of condition. 

# SQ pref mdl -------------------------------------------------------------

mdl_sq_pref <-
  brm(
    sleep_qual ~
      cndtn_pref * (alertness + tension + valence) +
      (cndtn_pref | ID),
    data = data_sleep,
    family = cumulative(link = "probit"), 
    prior = set_prior("normal(0, 1)", class = "b"),
    control = list(adapt_delta = 0.99), 
    file = "mdl_sq_pref"
    )

# Check posterior predictions
pp_check(mdl_sq_pref, ndraws = 100)
ggsave('sq_pref_pp.jpeg')

# Extract some effects of interest using marginaleffects
avg_predictions(mdl_sq_pref, by = "cndtn_pref", type = "link") 

# Plot on the link scale
plot_predictions(mdl_sq_pref, by = "cndtn_pref", type = "prediction") + 
  labs(x = 'Condition', y = 'SQ') +
  theme_classic() 
ggsave('sq_pref.jpeg')

# Look at pairwise comparisons
avg_comparisons(mdl_sq_pref, variables = list(cndtn_pref = "pairwise"), type = "link") |>
  print() |>
  posterior_draws(shape = "DxP") |>
  hypothesis(c("b1 > 0",
               "b2 > 0",
               "b3 > 0",
               "b4 < 0",
               "b5 < 0",
               "b6 > 0"))


# SOL pref mdl ------------------------------------------------------------

mdl_solpref_gamma <-
  brm(
    sleep_onset_add1 ~
      cndtn_pref * (alertness + tension + valence) +
      (cndtn | ID),
    data = data_sleep,
    family = "gamma", 
    prior = set_prior("normal(0, 1)", class = "b"),
    control = list(adapt_delta = 0.99), 
    file = "mdl_sol_pref"
    )

# Check posterior predictions
pp_check(mdl_solpref_gamma, ndraws = 100)
ggsave('sol_pref_pp.jpeg')

# Extract some effects of interest using marginaleffects
avg_predictions(mdl_solpref_gamma, by = "cndtn_pref", type = "link") 

# Plot on the link scale
plot_predictions(mdl_solpref_gamma, by = "cndtn_pref", type = "prediction") + 
  labs(x = 'Condition', y = 'SOL') +
  theme_classic() 
ggsave('sol_pref.jpeg')

# Look at pairwise comparisons
avg_comparisons(mdl_solpref_gamma, variables = list(cndtn_pref = "pairwise"), type = "link") |>
  print() |>
  posterior_draws(shape = "DxP") |>
  hypothesis(c("b1 > 0",
               "b2 > 0",
               "b3 > 0",
               "b4 < 0",
               "b5 < 0",
               "b6 > 0"))


# SE pref mdl -------------------------------------------------------------

mdl_sepref_gamma_trunc <-
  brm(
    sleep_effncy_inv | trunc(ub = 100) ~
      cndtn_pref * (alertness + tension + valence) +
      (cndtn | ID),
    data = data_sleep,
    family = "gamma", 
    prior = set_prior("normal(0, 1)", class = "b"),
    control = list(adapt_delta = 0.99), 
    file = "mdl_se_pref"
    )

# Check posterior predictions
pp_check(mdl_sepref_gamma_trunc, ndraws = 100)
ggsave('se_pref_pp.jpeg')

# Extract some effects of interest using marginaleffects
avg_predictions(mdl_sepref_gamma_trunc, by = "cndtn_pref", type = "link") 

# Plot on the link scale
plot_predictions(mdl_sepref_gamma_trunc, by = "cndtn_pref", type = "prediction") + 
  labs(x = 'Condition', y = 'SE') +
  theme_classic() 
ggsave('se_pref.jpeg')

# Look at pairwise comparisons
avg_comparisons(mdl_sepref_gamma_trunc, variables = list(cndtn_pref = "pairwise"), type = "link") |>
  print() |>
  posterior_draws(shape = "DxP") |>
  hypothesis(c("b1 > 0",
               "b2 > 0",
               "b3 > 0",
               "b4 < 0",
               "b5 < 0",
               "b6 > 0"))


