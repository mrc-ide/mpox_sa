library(tidyverse)
source("support.R")

# Data
data <- list(
  mpox_cases = list(
    sa_observed_ahiv = 15,
    sa_observed = 22,
    cd4_u100 = 85,      # Global Mitja et al
    cd4_u200 = 179,     # Global Mitja et al
    us_hiv_pos = 5350,  # Eustacio et al
    us_total = 9492  ), # Eustacio et al
  mpox_deaths = list(
    sa_observed = 3,
    cd4_u100 = 23,      # Global Mitja et al
    cd4_u200 = 27)      # Global Mitja et al
)

# Probabilities
probs <- list(
  hiv_us = c(q2.5 = 10910.5, mean = 12372.9, q97.5 = 13835.3) / 1e5,  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9668254/
  hiv_sa = c(q2.5 = 0.255, mean = 0.294, q97.5 = 0.335),              # Thembisa v4.7
  cd4_u200_sa = 0.066
)
probs$cd4_u100_sa <- probs$cd4_u200_sa / 2 # Assume half CD4<200/mm3 are <100/mm3

# Begin simulation
n <- 1e4
samples <- data.frame(sample = seq_len(n))
set.seed(123)

# Compute CFR for CD4 counts
samples <- samples %>%
  mutate(
    cfr_cd4_u100 = rbeta(n, data$mpox_deaths$cd4_u100,
                         data$mpox_cases$cd4_u100 - data$mpox_deaths$cd4_u100),
    cfr_cd4_u200 = rbeta(n, data$mpox_deaths$cd4_u200,
                         data$mpox_cases$cd4_u200 - data$mpox_deaths$cd4_u200)
  )

mean_ci(samples$cfr_cd4_u100)

# Calculate beta distribution parameters
p_hiv_sa_beta_params <- fit_beta(probs$hiv_sa["mean"], probs$hiv_sa["q2.5"],
                                 probs$hiv_sa["q97.5"])
p_hiv_us_beta_params <- fit_beta(probs$hiv_us["mean"], probs$hiv_us["q2.5"],
                                 probs$hiv_us["q97.5"])

# Generate samples
samples <- samples %>%
  mutate(
    p_hiv_sa = rbeta(n, p_hiv_sa_beta_params$alpha,
                     p_hiv_sa_beta_params$beta),
    p_hiv_us = rbeta(n, p_hiv_us_beta_params$alpha,
                     p_hiv_us_beta_params$beta),
    odds_hiv_us = odds(p_hiv_us),
    odds_hiv_sa = odds(p_hiv_sa),
    p_hiv_mpox_us = rbeta(n, data$mpox_cases$us_hiv_pos,
                          data$mpox_cases$us_total - data$mpox_cases$us_hiv_pos),
    odds_hiv_mpox_us = odds(p_hiv_mpox_us)
  )

# Calculate relative odds
samples <- samples %>%
  mutate(
    relative_odds_hiv = odds_hiv_mpox_us / odds_hiv_us,
    odds_hiv_mpox_sa = odds_hiv_sa * relative_odds_hiv,
    p_hiv_mpox_sa = odds_hiv_mpox_sa / (1 + odds_hiv_mpox_sa)
  )

# Calculate cases for different scenarios
samples <- samples %>%
  mutate(
    cases_sa_ahiv_u100_raw = data$mpox_deaths$sa_observed / cfr_cd4_u100,
    cases_sa_ahiv_u200_raw = data$mpox_deaths$sa_observed / cfr_cd4_u200,
    cases_sa_ahiv_u100 = pmax(data$mpox_cases$sa_observed_ahiv, cases_sa_ahiv_u100_raw),
    cases_sa_ahiv_u200 = pmax(data$mpox_cases$sa_observed_ahiv, cases_sa_ahiv_u200_raw),    
    cases_sa_hiv_given_u100 = cases_sa_ahiv_u100 / probs$cd4_u100_sa,
    cases_sa_hiv_given_u200 = cases_sa_ahiv_u200 / probs$cd4_u200_sa,
    cases_sa_pessimistic_u100 = cases_sa_hiv_given_u100 / p_hiv_sa,
    cases_sa_pessimistic_u200 = cases_sa_hiv_given_u200 / p_hiv_sa,
    cases_sa_neutral_u100 = cases_sa_hiv_given_u100 / p_hiv_mpox_sa,
    cases_sa_neutral_u200 = cases_sa_hiv_given_u200 / p_hiv_mpox_sa,
    ifr_sa_neutal_u100 = data$mpox_deaths$sa_observed / cases_sa_neutral_u100 * 100,
    ifr_sa_neutal_u200 = data$mpox_deaths$sa_observed / cases_sa_neutral_u200 * 100,
    ifr_sa_pessimistic_u100 = data$mpox_deaths$sa_observed / cases_sa_pessimistic_u100 * 100,
    ifr_sa_pessimistic_u200 = data$mpox_deaths$sa_observed / cases_sa_pessimistic_u200 * 100,
    pc_detected_sa_neutal_u100 = data$mpox_cases$sa_observed / cases_sa_neutral_u100 * 100,
    pc_detected_sa_neutal_u200 = data$mpox_cases$sa_observed / cases_sa_neutral_u200 * 100,
    pc_detected_sa_pessimistic_u100 = data$mpox_cases$sa_observed / cases_sa_pessimistic_u100 * 100,
    pc_detected_sa_pessimistic_u200 = data$mpox_cases$sa_observed / cases_sa_pessimistic_u200 * 100,
  )

# Summarize results
summary <- samples %>%
  select(-sample) %>% 
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(value = sprintf("%.3f (%.3f - %.3f)",
                            mean(value),
                            quantile(value, 0.025),
                            quantile(value, 0.975)))

print(summary, n = Inf)

