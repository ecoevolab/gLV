# =============================================================================
# Graph generation for controls simulation, different magnitude method.
# =============================================================================
#   Author: Manuel Rivera
#   Date:   2024-06-17
#   Description: Generates control method interaction graphs for slides.
#
# ========================================================================

library(igraph)
library(ggplot2)
library(ggraph)
library(ggtext)

# Generate common parameters for all methods
n_species <- 5
k         <- sample(n_species, 1)          # sample(1:n) == sample(n)
p_noint   <- runif(1, min = 0, max = 0.9)
p_neg     <- 1

# Interaction counts
total     <- (n_species - 1)^2
num_noint <- floor(p_noint * total)
num_negs  <- floor(p_neg   * (total - num_noint))
num_pos   <- total - num_noint - num_negs

# Build matrix
node_names <- paste0("S", seq_len(n_species))
M <- matrix(NA, n_species, n_species, dimnames = list(node_names, node_names))

# Precompute mask once instead of calling row()/col() twice
off_diag_off_k        <- col(M) != k & row(M) != col(M)
M[off_diag_off_k]     <- sample(c(
    rep(0,  num_noint),
    -runif(num_negs, 1, 2),
     runif(num_pos,  0, 1)
))
M[, k] <- -runif(n_species, 0, 1)
diag(M) <- -0.5
M        <- round(M, 5)

#---------------- Save Plots ----------------
make_grid_layout <- function(g, min_spacing = 1.5) {
    n <- igraph::vcount(g)
    
    # Compute grid dimensions
    ncols <- ceiling(sqrt(n))
    nrows <- ceiling(n / ncols)
    # Generate grid points
    x <- rep(seq(0, (ncols - 1) * min_spacing, by = min_spacing), length.out = n)
    y <- rep(seq(0, (nrows - 1) * min_spacing, by = min_spacing), each = ncols)[1:n]
    # Build layout data frame (required by ggraph)
    layout <- data.frame(
        x     = x,
        y     = y,
        name  = V(g)$name
    )
    
    return(layout)
}


graph_eg <- function(mat, highlight_node = NULL, highlight_color = "orange", highlight_edges = NULL) {
    g <- igraph::graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE)
    V(g)$name     <- rownames(mat)
    V(g)$key_node <- ifelse(V(g)$name %in% highlight_node, "highlight", "normal")
    E(g)$group    <- "normal"

    # Vectorized: pass all from/to pairs at once instead of looping
    if (!is.null(highlight_edges)) {
        eids <- igraph::get_edge_ids(g, t(rbind(highlight_edges$from, highlight_edges$to)))
        E(g)$group[eids[eids > 0]] <- "highlight"
    }

    layout <- make_grid_layout(g, min_spacing = 2.0)

    ggraph::ggraph(g, layout = "manual", x = layout$x, y = layout$y) +
        geom_edge_arc(
            aes(color = group, width = group),
            arrow      = arrow(length = unit(3, "mm")),
            end_cap    = circle(6, "mm"),
            strength   = 0.35,
            angle_calc = "along"
        ) +
        scale_edge_color_manual(values = c("normal" = "lightblue", "highlight" = "darkred"), guide = "none") +
        scale_edge_width_manual(values  = c("normal" = 0.6,         "highlight" = 1.2),       guide = "none") +
        geom_node_point(aes(fill = key_node), shape = 21, size = 12, color = "grey25", stroke = 1.5) +
        scale_fill_manual(values = c("highlight" = highlight_color, "normal" = "steelblue"), guide = "none") +
        geom_node_text(aes(label = name), color = "white", size = 4, fontface = "bold") +
        theme_graph(base_family = "sans") +
        labs(
            title    = "Ejemplo de red de interacción",
            subtitle = "Método de magnitudes negativas"
        ) +
        theme(
            plot.title    = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
        ) +
        annotate(
            "richtext",
            x = Inf, y = Inf,
           label = "<span style='color:darkred'>&#9632;</span> U(-1, 0)<br>
         <span style='color:steelblue'>&#9632;</span> U(-2, -1)",
            hjust = 1.05, vjust = 1.5, size = 3,
            fill  = "#FFFDD0", label.color = "grey80",
            label.r = unit(0.15, "lines")
        )
}

edges_to_highlight <- data.frame(
    from = paste0("S", k),   # whom_rows
    to   =  paste0("S", seq_len(n_species)[-k])
)

# generate the plot
save_dir = '/mnt/data/sur/users/mrivera/thesis_plots'
p = graph_eg(M, highlight_node =  paste0("S", k), highlight_color = "orange", highlight_edges = edges_to_highlight)
print(p)
ggsave(paste0(save_dir, "/magnitude_graph.png"), plot = p, width = 8, height = 6, dpi = 300)
# ggsave(paste0(save_dir, "/directed_graph.pdf"), plot = p, width = 8, height = 6)  # vector format

