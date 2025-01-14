
# -*- coding: utf-8 -*-

rm(list = ls())

source("/Users/dahunseo/programming/microbiome-logcontrast/Simulation/functions.R")

suppressPackageStartupMessages({
    library(ggtree)
    library(ape)
})

set.seed(1234)

# The logistic normal (LN) distribution
n <- 100            # n: sample size
p <- 70             # p: number of features
noise_sigma <- 0.3    # noise_sigma: noise level for response

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

expanded_z <- construct_features(z, tips = 1:p, edges = random_tree$edge, eta = 0, tau = 2)

tree_info <- create_tree_structure(tips = 1:p, edges = random_tree$edge)

non_zero_parent <- sample_non_zero_parent(tree_info, tau = 2, n_per_level = 4)

beta <- setNames(rep(0, length(colnames(expanded_z$data))), gsub("^.+_x", "", colnames(expanded_z$data)))

for (i in non_zero_parent) {
    update_indices <- as.character(tree_info$children[[i]])
    length <- length(update_indices)
    values <- generate_beta_vector(length)

    beta[update_indices] <- values
}

non_zero_node <- as.numeric(names(beta[beta != 0]))
names(beta) <- colnames(expanded_z$data)

y <- base_y + beta_treatment * treatment + beta_sex * sex + beta_age * age + as.vector(expanded_z$data %*% beta) + stats::rnorm(n, 0, sd = noise_sigma)

g <- ggtree::ggtree(random_tree, layout = "circular", branch.length="none", color = "black", size = 0.5, linetype = 1, alpha = 0.3)

height_length <- length(tree_info$levels)

g$data$level <- NA
for (h in seq_len(height_length)) {
    nodes <- tree_info$levels[[h]]
    g$data$level[g$data$node %in% nodes] <- h
}
g$data$level <- factor(g$data$level, levels = 1:height_length)

palette_function <- colorRampPalette(RColorBrewer::brewer.pal(n = min(height_length, 8), name = "Set1"))
colors <- palette_function(height_length)

g$data$label_text <- g$data$node
g$data$level_label <- factor(as.numeric(as.character(g$data$level)) - 1)
g$data$alpha_value <- 0.3
g$data$alpha_value[g$data$node %in% non_zero_node] <- 1

g +
    ggtree::geom_tippoint(aes(fill = level_label, alpha = alpha_value), size = 4, shape = 21, stroke = 0.5) +
    ggtree::geom_nodepoint(aes(fill = level_label, alpha = alpha_value), size = 4, shape = 21, stroke = 0.5) +
    ggtree::scale_fill_manual(
        values = colors,
        drop = FALSE
    ) +
    ggrepel::geom_text_repel(aes(label = label_text, alpha = alpha_value), size = 4, box.padding = 0.3, point.padding = 0.5, segment.linetype = "dotted", segment.alpha = 0.2) +
    ggtree::theme(legend.position = "right") +
    ggplot2::labs(fill = "Level") +
    ggplot2::guides(alpha = "none")

source("/Users/dahunseo/programming/microbiome-logcontrast/Biometrika_2014/ConstLasso.R")

C <- matrix(1, p, ncol = 1)

res2 <- ConstrLassoCrossVal(y = y, x = log_z, C = C, nfolds = 10)

plot(res2$cvm)

res2$bet.sel
t(C) %*% res2$bet.sel

res2$cvm[res2$sel]
res2$Rsq.sel

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
    beta, res3$bet.sel,
    xlab = "true beta", 
    ylab = "estimated beta",
    pch = 16,
    col = "blue",
    asp = 1
)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)}

{plot(
    y, expanded_z$data %*% res3$bet.sel + res3$int.sel,
    xlab = "true y", 
    ylab = "y hat",
    pch = 16,
    col = "blue",
    asp = 1
)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)}
