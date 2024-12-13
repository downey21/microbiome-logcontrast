
# -*- coding: utf-8 -*-

rm(list = ls())

source("/Users/dahunseo/programming/microbiome-logcontrast/Simulation/functions.R")

set.seed(1234)

# The logistic normal (LN) distribution
n <- 100            # n: sample size
p <- 70             # p: number of features
noise_sigma <- 1    # noise_sigma: noise level for response

# covariate generation
sex <- sample(c(0, 1), n, replace = TRUE)  # sex: 0 or 1
age <- rnorm(n, mean = 50, sd = 10)        # age: N(50, 10^2)

# treatment generation
treatment <- sample(c(0, 1), n, replace = TRUE) # treatment: 0 or 1

# parameters for normal distribution
# base_mu <- rep(0, p)

# parameters for normal distribution
# Some taxa are often significantly more abundant than others, as commonly seen in real microbiome compositional data.
base_mu <- c(rep(p/2, 5), rep(1, p-5)) 

# Create a random binary tree for the p features
random_tree <- ape::rcoal(p)

# Compute the distance matrix using cophenetic distances
dist_matrix <- stats::cophenetic(random_tree)

# variance-covariance setting
sigma <- exp(-dist_matrix) / 2

# variance-covariance setting
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

# add the psudeo-count 0.5
z <- ifelse(z == 0, 0.5, z)

# check
head(z)
apply(z, 1, sum)

# log transformation
log_z <- log(z)

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

beta_sex <- 0
beta_age <- 0
beta_treatment <- 0
base_y <- rep(100, n)

y <- base_y + beta_treatment * treatment + beta_sex * sex + beta_age * age + as.vector(log_z %*% beta) + stats::rnorm(n, 0, sd = noise_sigma)




tree_info <- create_tree_structure(tips = 1:p, edges = random_tree$edge)
tree_info$levels

g <- igraph::graph_from_edgelist(random_tree$edge, directed = FALSE)

layout <- igraph::layout_as_tree(g, root = tree_info$levels[[length(tree_info$levels)]], mode = "out")

plot(g,
    layout = layout,
    vertex.label = igraph::V(g)$name,  
    vertex.size = 5,         
    vertex.color = "skyblue", 
    edge.color = "black"
)

expanded_z <- construct_features(z, tips = 1:p, edges = random_tree$edge, eta = 0, tau = 1)

source("/Users/dahunseo/programming/microbiome-logcontrast/Biometrika_2014/ConstLasso.R")

C <- matrix(1, p, ncol = 1)

res2 <- ConstrLassoCrossVal(y = y, x = log_z, C = C, nfolds = 10)

plot(res2$cvm)

res2$bet.sel
t(C) %*% res2$bet.sel

res2$cvm[res2$sel]
res2$Rsq.sel

{plot(
    beta, res2$bet.sel,
    xlab = "true beta", 
    ylab = "estimated beta",
    pch = 16,
    col = "blue",
    asp = 1
)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)}

{plot(
    y, log_z %*% res2$bet.sel + res2$int.sel,
    xlab = "true y", 
    ylab = "y hat",
    pch = 16,
    col = "blue",
    asp = 1
)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)}

res3 <- ConstrLassoCrossVal(y = y, x = expanded_z$data, C = expanded_z$C, nfolds = 10)

plot(res3$cvm)

res3$bet[, which.min(res3$cvm)]

t(expanded_z$C) %*% res3$bet[, which.min(res3$cvm)]

expanded_z$data %*% res3$bet[,min(which(res3$cvm <= min(res3$cvm+res3$cvsd)))] + res3$int[min(which(res3$cvm <= min(res3$cvm+res3$cvsd)))]
expanded_z$data %*% res3$bet.sel + res3$int.sel

res3$cvm[res3$sel]
res3$Rsq.sel

{plot(
    y, expanded_z$data %*% res3$bet.sel + res3$int.sel,
    xlab = "true y", 
    ylab = "y hat",
    pch = 16,
    col = "blue",
    asp = 1
)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)}
