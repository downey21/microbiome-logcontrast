
# -*- coding: utf-8 -*-

rm(list = ls())

set.seed(1234)

# The logistic normal (LN) distribution
n <- 20            # n: sample size
p <- 7             # p: number of features
noise_sigma <- 1    # noise_sigma: noise level for response

# covariate generation
sex <- sample(c(0, 1), n, replace = TRUE)  # sex: 0 or 1
age <- rnorm(n, mean = 50, sd = 10)        # age: N(50, 10^2)

# treatment generation
treatment <- sample(c(0, 1), n, replace = TRUE) # treatment: 0 or 1

# parameters for normal distribution
base_mu <- rep(0, p)

# Create a random binary tree for the p features
random_tree <- ape::rcoal(p)
random_tree$tip.label <- paste0("x", 1:p)

# Compute the distance matrix using cophenetic distances
dist_matrix <- stats::cophenetic(random_tree)

# variance-covariance setting
sigma <- exp(-dist_matrix) / 2

# sigma <- matrix(0, p, p)
# gamma <- 0.5
# for(i in 1:nrow(sigma)){
#     for(j in 1:nrow(sigma)){
#         sigma[i,j] <- gamma^(abs(i-j))   
#     }
# }

# coefficients
# alpha_sex <- rnorm(p, mean = 0.2, sd = 0.1)
# alpha_age <- rnorm(p, mean = 0.01, sd = 0.005)
# alpha_treatment <- rnorm(p, mean = 0.2, sd = 0.1)

# coefficients
alpha_sex <- 0
alpha_age <- 0
alpha_treatment <- 0

# calculate mu for each i
mu_matrix <- matrix(NA, n, p)
for (i in 1:n) {
    mu_matrix[i, ] <- base_mu + alpha_treatment * treatment[i] + alpha_sex * sex[i] + alpha_age * age[i]
}

# generation of logistic normal samples
z <- matrix(NA, n, p)
for (i in 1:n) {
    normal_sample <- MASS::mvrnorm(1, mu_matrix[i, ], sigma)
    exp_sample <- exp(normal_sample)
    z[i, ] <- exp_sample / sum(exp_sample)
}

# label
colnames(z) <- paste0("x", 1:p)

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
