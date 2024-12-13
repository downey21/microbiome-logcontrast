
# -*- coding: utf-8 -*-

rm(list = ls())

library(MASS)
library(dirmult)

# Following settings in
# A fast small‐sample kernel independence test for microbiome community‐level association analysis, Biometrics (2017)
# Compositional knockoff filter for high‐dimensional regression analysis of microbiome data, Biometrics (2021) 에서 위 논문 DM follow.
# (The first scheme was to generate microbiome counts from the Dirichlet-multinomial (DM) distribution following a previous design (Zhan et al., 2017).)

# The Dirichlet-multinomial (DM) distribution
n <- 20            # n: sample size
p <- 7             # p: number of features
noise_sigma <- 1    # noise_sigma: noise level for response
depth <- stats::rnbinom(n, mu = 10000, size = 25)

# covariate generation
set.seed(1234)
sex <- sample(c(0, 1), n, replace = TRUE)  # sex: 0 or 1
age <- rnorm(n, mean = 50, sd = 10)        # age: N(50, 10^2)

# treatment generation
treatment <- sample(c(0, 1), n, replace = TRUE) # treatment: 0 or 1

# Create a random binary tree for the p features
random_tree <- ape::rcoal(p)
random_tree$tip.label <- paste0("x", 1:p)

# Compute the distance matrix using cophenetic distances
dist_matrix <- stats::cophenetic(random_tree)

# variance-covariance setting
cov_matrix <- exp(-dist_matrix) / 2

# coefficients
# base_alpha <- rep(1, p)
# alpha_sex <- rnorm(p, mean = 0.2, sd = 0.1)
# alpha_age <- rnorm(p, mean = 0.01, sd = 0.005)
# alpha_treatment <- rnorm(p, mean = 0.2, sd = 0.1)

# coefficients
# base_alpha <- rep(1, p)
# alpha_sex <- 0
# alpha_age <- 0
# alpha_treatment <- 0

# coefficients
base_alpha <- exp(MASS::mvrnorm(1, mu = rep(0, p), Sigma = cov_matrix))
alpha_sex <- 0
alpha_age <- 0
alpha_treatment <- 0

# calculate mu for each i
alpha_matrix <- matrix(NA, n, p)
for (i in 1:n) {
    alpha_matrix[i, ] <- base_alpha * exp(alpha_treatment * treatment[i] + alpha_sex * sex[i] + alpha_age * age[i])
}

# generation of Dirichlet-multinomial samples
z <- matrix(NA, nrow = n, ncol = p)
for (i in 1:n) {
    dirichlet_sample <- dirmult::rdirichlet(1, alpha_matrix[i, ])
    z[i, ] <- rmultinom(1, size = depth[i], prob = dirichlet_sample)
}

# label
colnames(z) <- paste0("x", 1:p)

# add the psudeo-count 0.5
z <- ifelse(z == 0, 0.5, z)

# calculate proportion
z <- z / rowSums(z)

# check
head(z)
apply(z, 1, sum)

ape::plot.phylo(random_tree, type = "fan", main = paste("Random Binary Tree with", p, "Variables"))

# log transformation
z <- log(z)

# coefficients
beta_non_zero <-
    c(
        -3, 3, 2.5, -1, -1.5,
        3, 3, -2, -2, -2,
        1, -1, 3, -2, -1
    )

if (length(beta_non_zero) >= p) {
    beta <- beta_non_zero[1:p]
} else {
    beta <- c(beta_non_zero, rep(0, p - length(beta_non_zero)))
}

beta_sex <- 1
beta_age <- 1
beta_treatment <- 1
base_y <- rep(1, n)

y <- base_y + beta_treatment * treatment + beta_sex * sex + beta_age * age + as.vector(z %*% beta) + stats::rnorm(n, 0, sd = noise_sigma)
