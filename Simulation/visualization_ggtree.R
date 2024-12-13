
# -*- coding: utf-8 -*-

rm(list = ls())

#
library(ggtree)
library(ape)

set.seed(123)

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

p <- 10

#
tree <- ape::rtree(p)
tree <- ape::rcoal(p)
# tree <- ape::read.tree("path_to_reference_tree.newick")

#
ggtree::ggtree(tree, layout = "rectangular", color = "black", size = 0.5, linetype = 1)
ggtree::ggtree(tree, layout = "rectangular", branch.length="none", color = "black", size = 0.5, linetype = 1)

ggtree::ggtree(tree, layout = "dendrogram", color = "black", size = 0.5, linetype = 1)
ggtree::ggtree(tree, layout = "dendrogram", branch.length="none", color = "black", size = 0.5, linetype = 1)

ggtree::ggtree(tree, layout = "slanted", color = "black", size = 0.5, linetype = 1)
ggtree::ggtree(tree, layout = "slanted", branch.length="none", color = "black", size = 0.5, linetype = 1)

ggtree::ggtree(tree, layout = "ellipse", color = "black", size = 0.5, linetype = 1)
ggtree::ggtree(tree, layout = "ellipse", branch.length="none", color = "black", size = 0.5, linetype = 1)

ggtree::ggtree(tree, layout = "roundrect", color = "black", size = 0.5, linetype = 1)
ggtree::ggtree(tree, layout = "roundrect", branch.length="none", color = "black", size = 0.5, linetype = 1)

ggtree::ggtree(tree, layout = "circular", color = "black", size = 0.5, linetype = 1)
ggtree::ggtree(tree, layout = "circular", branch.length="none", color = "black", size = 0.5, linetype = 1)

ggtree::ggtree(tree, layout = "fan", color = "black", size = 0.5, linetype = 1)
ggtree::ggtree(tree, layout = "fan", branch.length="none", color = "black", size = 0.5, linetype = 1)

#
g <- ggtree::ggtree(tree, layout = "circular", branch.length="none", color = "black", size = 0.5, linetype = 1, alpha = 0.3)

g +
    ggtree::geom_nodepoint(color = "blue", size = 4, alpha = 0.9, shape = 16) +
    ggtree::geom_tippoint(color = "red", size = 4, alpha = 0.9, shape = 16)

tree_info <- create_tree_structure(tips = 1:p, edges = tree$edge)
tree_info$levels

height_length <- length(tree_info$levels)

g$data$level <- NA
for (h in seq_len(height_length)) {
    nodes <- tree_info$levels[[h]]
    g$data$level[g$data$node %in% nodes] <- h
}
g$data$level <- factor(g$data$level, levels = 1:height_length)



palette_function <- colorRampPalette(RColorBrewer::brewer.pal(n = min(height_length, 8), name = "Set1"))
colors <- palette_function(height_length)

#
g +
    ggtree::geom_tippoint(aes(color = level), size = 4, alpha = 0.9, shape = 16) +
    ggtree::geom_nodepoint(aes(color = level), size = 4, alpha = 0.9, shape = 16) +
    ggtree::scale_color_manual(values = colors) +
    ggtree::geom_text(aes(label=node), hjust = -0.5, vjust = -0.5, size = 4) +
    ggtree::theme(legend.position = "right") +
    ggplot2::labs(color = "Level")

g +
    ggtree::geom_tippoint(aes(color = level), size = 4, alpha = 0.9, shape = 16) +
    ggtree::geom_nodepoint(aes(color = level), size = 4, alpha = 0.9, shape = 16) +
    ggtree::scale_color_manual(values = colors) +
    ggrepel::geom_text_repel(aes(label = node), size = 4, box.padding = 0.3, point.padding = 0.5) +
    ggtree::theme(legend.position = "right") +
    ggplot2::labs(color = "Level")

#
g$data$label_text <- NA

g$data$label_text[g$data$node == 1] <- "A"
g$data$label_text[g$data$node == 2] <- "B"
g$data$label_text[g$data$node == 3] <- "C"


g$data$label_text[g$data$node == 4] <- "D"

g$data$label_text[is.na(g$data$label_text)] <- ""

g +
    ggtree::geom_tippoint(aes(color = level), size = 4, alpha = 0.9, shape = 16) +
    ggtree::geom_nodepoint(aes(color = level), size = 4, alpha = 0.9, shape = 16) +
    ggtree::scale_color_manual(values = colors) +
    ggrepel::geom_text_repel(aes(label = label_text), size = 4, box.padding = 0.3, point.padding = 0.5) +
    ggtree::theme(legend.position = "right") +
    ggplot2::labs(color = "Level")

#
g$data$level_label <- paste("Group", g$data$level)

g$data$level_label <- factor(
    g$data$level_label,
    levels = paste("Group", 1:height_length)
)

g +
    ggtree::geom_tippoint(aes(color = level_label), size = 4, alpha = 0.9, shape = 16) +
    ggtree::geom_nodepoint(aes(color = level_label), size = 4, alpha = 0.9, shape = 16) +
    ggtree::scale_color_manual(
        values = colors,
        drop = FALSE
    ) +
    ggrepel::geom_text_repel(aes(label = label_text), size = 4, box.padding = 0.3, point.padding = 0.5) +
    ggtree::theme(legend.position = "right") +
    ggplot2::labs(color = "Level")

#
g +
    ggtree::geom_tippoint(aes(fill = level_label), size = 4, color = "black", stroke = 0.5, alpha = 1, shape = 21) +
    ggtree::geom_nodepoint(aes(fill = level_label), size = 4, color = "black", stroke = 0.5, alpha = 1, shape = 21) +
    ggtree::scale_fill_manual(
        values = colors,
        drop = FALSE
    ) +
    ggrepel::geom_text_repel(aes(label = label_text), size = 4, box.padding = 0.3, point.padding = 0.5) +
    ggtree::theme(legend.position = "right") +
    ggplot2::labs(fill = "Level")

#
g$data$highlight <- ifelse(g$data$node %in% c(1, 6), "highlight", "normal")

g +
    ggtree::geom_tippoint(
        data = g$data[g$data$highlight == "highlight", ],
        aes(fill = level_label),
        color = "red",
        size = 6,
        shape = 21,
        stroke = 2,
        alpha = 0.9,
        show.legend = FALSE
    ) +
    ggtree::geom_nodepoint(
        data = g$data[g$data$highlight == "highlight", ],
        aes(fill = level_label),
        color = "red",
        size = 6,
        shape = 21,
        stroke = 2,
        alpha = 0.9,
        show.legend = FALSE
    ) +
    ggtree::geom_tippoint(aes(fill = level_label), size = 4, color = "black", stroke = 0.5, alpha = 1, shape = 21) +
    ggtree::geom_nodepoint(aes(fill = level_label), size = 4, color = "black", stroke = 0.5, alpha = 1, shape = 21) +
    ggtree::scale_fill_manual(
        values = colors,
        drop = FALSE
    ) +
    ggrepel::geom_text_repel(aes(label = label_text), size = 4, box.padding = 0.3, point.padding = 0.5) +
    ggtree::theme(legend.position = "right") +
    ggplot2::labs(fill = "Level")

#
g$data$alpha_value <- 0.3
g$data$alpha_value[g$data$node == 1] <- 0.8
g$data$alpha_value[g$data$node == 6] <- 1

g +
    ggtree::geom_tippoint(aes(fill = level_label, alpha = alpha_value), color = "black", size = 4, shape = 21, stroke = 0.5) +
    ggtree::geom_nodepoint(aes(fill = level_label, alpha = alpha_value), color = "black", size = 4, shape = 21, stroke = 0.5) +
    ggplot2::scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
    ggtree::scale_fill_manual(
        values = colors,
        drop = FALSE
    ) +
    ggrepel::geom_text_repel(aes(label = label_text), size = 4, box.padding = 0.3, point.padding = 0.5) +
    ggtree::theme(legend.position = "right") +
    ggplot2::labs(fill = "Level")






