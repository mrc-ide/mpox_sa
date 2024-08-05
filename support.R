# Calculate Binomial mean and exact 95% CI
calc_p <- function(x, n) {
  ci <- binom.test(x, n)$conf.int
  c(ci[1], x/n, ci[2])
}

# Odds calculation function
odds <- function(p) p / (1 - p)

# Fit beta distribution function via MOM
fit_beta <- function(mean, lower_ci, upper_ci, ci_level = 0.95) {
  z_score <- qnorm(1 - (1 - ci_level) / 2)
  se <- (upper_ci - lower_ci) / (2 * z_score)
  variance <- se ^ 2
  alpha <- mean * ((mean * (1 - mean) / variance) - 1)
  beta <- (1 - mean) * ((mean * (1 - mean) / variance) - 1)
  list(alpha = alpha, beta = beta)
}

# Empirical mean and 95% CI
mean_ci <- function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975)))

# Parametric Mean and 95% CI of beta dist
beta_mean_ci <- function(alpha, beta, format = TRUE) {
  x <- c(mean = alpha / beta, qbeta(c(0.025, 0.975), alpha, beta)) * 100
  if (format) sprintf("%.1f (%.1f - %.1f)", x[1], x[2], x[3])
}
