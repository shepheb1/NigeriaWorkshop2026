# Day 1 Simulation Exercise

rm(list = ls())
library(MASS)
set.seed(123)

n <- 100
rho <- 0.5
beta <- 0.5

x <- c(rep(0, n / 2), rep(1, n / 2))
mu <- c(5, 5)

Sigma <- matrix(
  c(1, rho,
    rho, 1),
  nrow = 2,
  byrow = TRUE
)

# Generate Y0 and Y1 
y <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
y

#vc Baseline outcome
y0 <- y[, 1] 

# Post-baseline outcome with treatment effect added
y1 <- y[, 2] + beta * x

# Change from baseline
change <- y1 - y0

# Create data frame
dat <- data.frame(
  id = 1:n,
  x = x,
  y0 = y0,
  y1 = y1,
  change = change
)
head(dat)
tail(dat)


# a. T-test comparing Y1 for X = 0 vs X = 1
fitA <- t.test(y1 ~ x, data = dat)
estA <- mean(dat$y1[dat$x == 1]) - mean(dat$y1[dat$x == 0])
pA <- fitA$p.value

# b. T-test comparing Y1 - Y0 for X = 0 vs X = 1
fitB <- t.test(change ~ x, data = dat)
estB <- mean(dat$change[dat$x == 1]) - mean(dat$change[dat$x == 0])
pB <- fitB$p.value

# c. Linear model regressing Y1 on X
fitC <- lm(y1 ~ x, data = dat)
estC <- coef(fitC)["x"]
summary(fitC)
pC <- summary(fitC)$coefficients["x", "Pr(>|t|)"]

# d. Linear model regressing Y1 on X and Y0 (ANCOVA model)
fitD <- lm(y1 ~ x + y0, data = dat)
summary(fitD)
estD <- coef(fitD)["x"]
pD <- summary(fitD)$coefficients["x", "Pr(>|t|)"]

# e. Linear model regressing Y1 - Y0 on X
fitE <- lm(change ~ x, data = dat)
summary(fitE)
estE <- coef(fitE)["x"]
pE <- summary(fitE)$coefficients["x", "Pr(>|t|)"]

# f. Linear model regressing Y1 - Y0 on X and Y0 (ANCOVA model)
fitF <- lm(change ~ x + y0, data = dat)
summary(fitF) 
estF <- coef(fitF)["x"]
pF <- summary(fitF)$coefficients["x", "Pr(>|t|)"]


#all results together
results <- data.frame(
  method = c(
    "A: t-test comparing Y1 by X",
    "B: t-test comparing change by X",
    "C: Linear model: Y1 ~ X",
    "D: Linear model: Y1 ~ X + Y0",
    "E: Linear model: change ~ X",
    "F: Linear model: change ~ X + Y0"
  ),
  
  estimate = c(estA, estB, estC, estD, estE, estF),
  
  p_value = c(pA, pB, pC, pD, pE, pF)
)
#results
results

#summaries
fitA
fitB
summary(fitC)
summary(fitD)
summary(fitE)
summary(fitF)


############################################################
# 2. Repeat 1000 times to estimate the power of each analysis strategy.
#Vary ρ ∈ {0,0.25,0.50,0.75,0.90}.
############################################################

nsims <- 1000
n <- 100
beta <- 0.5
rho_values <- c(0, 0.25, 0.50, 0.75, 0.90)

# Empty data frame to store final results
all_results <- data.frame()

# Loop over rho values
for (rho in rho_values) {
  # Storage for estimates and p-values for this rho
  estimates <- matrix(NA, nrow = nsims, ncol = 6)
  p_values  <- matrix(NA, nrow = nsims, ncol = 6)
  
  colnames(estimates) <- c("A", "B", "C", "D", "E", "F")
  colnames(p_values)  <- c("A", "B", "C", "D", "E", "F")
  
  # Repeat simulation 1000 times
  for (i in 1:nsims) {
    x <- c(rep(0, n / 2), rep(1, n / 2))
    mu <- c(5, 5)
    Sigma <- matrix(
      c(1, rho,
        rho, 1),
      nrow = 2,
      byrow = TRUE
    )
    y <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
    y0 <- y[, 1]
    y1 <- y[, 2] + beta * x
    change <- y1 - y0
    dat <- data.frame(
      id = 1:n,
      x = x,
      y0 = y0,
      y1 = y1,
      change = change
    )
    
    # a. T-test comparing Y1 by X
    fitA <- t.test(y1 ~ x, data = dat)
    estimates[i, "A"] <- mean(dat$y1[dat$x == 1]) - mean(dat$y1[dat$x == 0])
    p_values[i, "A"] <- fitA$p.value
    
    # b. T-test comparing change by X
    fitB <- t.test(change ~ x, data = dat)
    estimates[i, "B"] <- mean(dat$change[dat$x == 1]) - mean(dat$change[dat$x == 0])
    p_values[i, "B"] <- fitB$p.value
    
    # c. Linear model: Y1 ~ X
    fitC <- lm(y1 ~ x, data = dat)
    estimates[i, "C"] <- coef(fitC)["x"]
    p_values[i, "C"] <- summary(fitC)$coefficients["x", "Pr(>|t|)"]
    
    # D. Linear model: Y1 ~ X + Y0
    fitD <- lm(y1 ~ x + y0, data = dat)
    estimates[i, "D"] <- coef(fitD)["x"]
    p_values[i, "D"] <- summary(fitD)$coefficients["x", "Pr(>|t|)"]
    
    # E. Linear model: change ~ X
    fitE <- lm(change ~ x, data = dat)
    estimates[i, "E"] <- coef(fitE)["x"]
    p_values[i, "E"] <- summary(fitE)$coefficients["x", "Pr(>|t|)"]
    
    # F. Linear model: change ~ X + Y0
    fitF <- lm(change ~ x + y0, data = dat)
    estimates[i, "F"] <- coef(fitF)["x"]
    p_values[i, "F"] <- summary(fitF)$coefficients["x", "Pr(>|t|)"]
  }

  # Calculate power and mean estimate for this rho
  mean_estimate <- colMeans(estimates)
  power <- colMeans(p_values < 0.05)
  
  result_rho <- data.frame(
    rho = rho,
    method = c(
      "A: t-test comparing Y1 by X",
      "B: t-test comparing change by X",
      "C: Linear model: Y1 ~ X",
      "D: Linear model: Y1 ~ X + Y0",
      "E: Linear model: change ~ X",
      "F: Linear model: change ~ X + Y0"
    ),
    mean_estimate = mean_estimate,
    power = power
  )
  
  all_results <- rbind(all_results, result_rho)
}

# Final results
all_results


