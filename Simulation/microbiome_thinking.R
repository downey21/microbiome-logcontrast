
# -*- encoding: utf-8 -*-

set.seed(48105)
n <- 500
p <- 3

Omega_inv <- diag(p)
for(i in 1:p) {
    for(j in 1:p) {
        Omega_inv[i,j] <- 0.7^abs(i-j)
    }
}    
out <- base::eigen(Omega_inv, symmetric = TRUE)
Omega_inv_sprt <- base::tcrossprod(out$vec*rep(out$val^(0.5), each=p), out$vec)
Omega <- base::tcrossprod(out$vec*rep(out$val^(-1), each=p), out$vec)
round(Omega, 2)

Y <- matrix(rnorm(n*p), nrow=n, ncol=p) %*% Omega_inv_sprt
apply(Y, 2, mean)

cov(Y)
Omega_inv

plot_hist_with_density <- function(data, omega_inv) {
    
    sd_vector <- sqrt(diag(omega_inv))
    
    par(mfrow = c(1, ncol(data)))
    
    for (i in 1:ncol(data)) {
        hist(data[, i], main = paste("Y", i), xlab = "", prob = TRUE)
        curve(dnorm(x, mean = 0, sd = sd_vector[i]), col = "red", lwd = 2, add = TRUE)
    }
    
    par(mfrow = c(1, 1))
}

plot_hist_with_density(data = Y, omega_inv = Omega_inv)
