
create_tree_structure <- function(tips, edges) {
    
    all_nodes <- unique(c(edges))
    n_nodes <- length(all_nodes)
    
    parents <- vector("list", n_nodes)
    children <- vector("list", n_nodes)
    
    for (i in seq_len(nrow(edges))) {
        parent <- edges[i, 1]
        child <- edges[i, 2]
        
        if (is.null(children[[parent]])) {
            children[[parent]] <- c(child)
        } else {
            children[[parent]] <- c(children[[parent]], child)
        }
        
        parents[[child]] <- parent
    }
    
    levels <- list()
    node_levels <- rep(NA, n_nodes)
    
    for (tip in tips) {
        node_levels[tip] <- 1
    }
    levels[[1]] <- tips
    
    remaining_nodes <- setdiff(all_nodes, tips)
    while (length(remaining_nodes) > 0) {
        for (node in remaining_nodes) {

            child_nodes <- children[[node]]
            if (!is.null(child_nodes) && all(!is.na(node_levels[child_nodes]))) {
            
                node_levels[node] <- max(node_levels[child_nodes]) + 1
                
                current_level <- node_levels[node]
                if (length(levels) < current_level) {
                    levels[[current_level]] <- c(node)
                } else {
                    levels[[current_level]] <- c(levels[[current_level]], node)
                }
            }
        }
        
        remaining_nodes <- remaining_nodes[is.na(node_levels[remaining_nodes])]
    }
    
    return(list(levels = levels, parents = parents, children = children))
}

construct_features <- function(x, tips, edges, eta = 0, tau = NULL) {
    
    tree_structure <- create_tree_structure(tips, edges)
    levels <- tree_structure$levels
    parents <- tree_structure$parents
    children <- tree_structure$children
    
    C <- NULL

    max_height <- length(levels) - 1
    if (is.null(tau) || tau > max_height - 1) {
        tau <- max_height - 1
    }
    
    if (eta < 0 || eta > tau) {
        stop("Invalid range for eta and tau")
    }
    
    n_samples <- nrow(x)
    x <- ifelse(x == 0, 0.5, x)
    z <- matrix(NA, nrow = n_samples, ncol = length(parents))
    z[, 1:ncol(x)] <- x
    
    features_output <- unlist(levels[(eta + 1):(tau + 1)])

    p <- length(features_output)
    transformed_flag <- rep(FALSE, length(parents))
    
    if (!(eta == 0 & tau == 0)) {

        for (h in seq(1, tau + 1)) {
        
            current_nodes <- levels[[h + 1]]

            for (parent in current_nodes) {
                
                temp_C <- matrix(0, nrow = p, ncol = 1)

                child_indices <- children[[parent]]
                
                if (!is.null(child_indices)) {
                    child_features <- z[, child_indices]
                    parent_abundance <- rowSums(child_features, na.rm = TRUE)
                    z[, parent] <- parent_abundance
                    
                    z[, child_indices] <- log(z[, child_indices] / parent_abundance)
                    transformed_flag[child_indices] <- TRUE

                    temp_C[which(features_output %in% child_indices), ] <- 1
                    C <- cbind(C, temp_C)
                }
            }
        }

    } else {

        C <- matrix(1, nrow = p, ncol = 1)

    }

    z[,!transformed_flag] <- log(z[,!transformed_flag])
    
    relevant_columns <- features_output
    expanded_features <- z[, relevant_columns, drop = FALSE]
    
    colnames(expanded_features) <-
        paste0(
            "h", rep(eta:tau, times = sapply(levels[(eta + 1):(tau + 1)], length)),
            "_x", relevant_columns
        )
    
    return(list(data = expanded_features, C = C))
}

generate_beta_vector <- function(length) {

    values <- seq(-2, 2, by = 0.25)
    values <- values[values != 0]

    if (length == 1) {
        output <- c(sample(values, 1, replace = TRUE))
    } else {
        selected_values <- sample(values, length - 1, replace = TRUE)
        last_value <- -sum(selected_values)

        output <- c(selected_values, last_value)
    }

    return(output)
}

sample_with_exclusion <- function(level, exclude, n) {
    candidates <- setdiff(level, exclude)
    if (length(candidates) == 0) {
        return(integer(0))
    } else if (length(candidates) < n) {
        return(candidates)
    } else {
        return(sample(candidates, n))
    }
}

sample_non_zero_parent <- function(tree_info, tau, n_per_level = 2, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)

    selected <- c()
    exclude <- c()

    for (level in rev(2:(tau+2))) {
        current_level <- tree_info$levels[[level]]
        sampled <- sample_with_exclusion(current_level, exclude, n_per_level)
        selected <- c(selected, sampled)
        
        for (i in sampled) {
            exclude <- c(exclude, tree_info$children[[i]])
        }
    }

    return(selected)
}