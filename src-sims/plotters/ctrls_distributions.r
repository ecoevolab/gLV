# 
# Description: Sctipt to generate plots for all available controls methods

#-------------------------------------
# Section: Load data
#-------------------------------------
library(dplyr)
library(parallel)

process_row <- function(i, exp, sim_params) {
    row <- sim_params[i, ]
    arrow::read_feather(
        file.path(exp, paste0("ExtSummaries/ExtSummary_", row$id, ".feather")) ) |>
        rename(dissimilarity = dissimilarity_bc) |>
        mutate(
            is_keystone = row_number() == as.numeric(row$key),
            p_noint     = as.numeric(row$p_noint),
            specie      = row_number()
        ) |>
        select(-pop_initial, -p_noint, -time_stability)
}

#-------------------------------------
# Section: Metric distributions
#-------------------------------------
library(ggplot2)
library(tidyr)
library(ggpubr)

# Centralize palette, variables, labels
KEYSTONE_COLORS  <- c("FALSE" = "#4575b4", "TRUE" = "#d73027")
KEYSTONE_LABELS  <- c("FALSE" = "No esencial", "TRUE" = "Esencial")
METRICS          <- c("keystoneness", "dissimilarity", "prop_extinctions")
METRIC_LABELS <- c(keystoneness     = "Keystoneness",
                   dissimilarity    = "Disimilitud BC",
                   prop_extinctions = "Prop. Extinciones")

theme_custom <- theme_minimal(base_size = 12) +
    theme(
        legend.position  = "none",
        plot.title       = element_text(face = "bold",   size = 11, hjust = 0.5),
        plot.subtitle    = element_text(face = "italic", size = 9,  hjust = 0.5),
        axis.title       = element_text(size = 10),
        strip.text       = element_text(face = "bold",   size = 10),
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 1.5)
    )

# Shared helper: pivot + factor
pivot_metrics <- function(data, metrics = METRICS) {
    data |>
        pivot_longer(cols = all_of(metrics), names_to = "variable", values_to = "value") |>
        mutate(variable = factor(variable, levels = metrics))
}

# Shared helper: common scales
scales_keystone <- function() {
    list(
        scale_fill_manual(values  = KEYSTONE_COLORS),
        scale_color_manual(values = KEYSTONE_COLORS),
        scale_y_continuous(trans  = "log1p"),
        facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = METRIC_LABELS))
    )
}

make_boxplot <- function(data, title, subtitle, metrics = METRICS) {
    plot_data <- pivot_metrics(data, metrics)
    plot_data |>
        ggplot(aes(x = is_keystone, y = value, fill = is_keystone)) +
        geom_jitter(aes(color = is_keystone), width = 0.15, size = 0.3, alpha = 0.15, show.legend = FALSE) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA) +
        # Move stars higher
        stat_compare_means(method = "wilcox.test", label = "p.signif",label.x = 1.5, label.y.npc = 0.99) +
        scale_x_discrete(labels = KEYSTONE_LABELS) +
        scales_keystone() +
        labs(title = title, subtitle = subtitle, x = "Estado de especie", y = "Valor de la métrica (escala log1p)") +
        theme_custom
}

#-------------------------------------
# Section: Metrics vs Relative Abundance for keystone species
#-------------------------------------
library(ggtext)  

met_relabu <- function(data, threshold = 1e-6, title, subtitle, metrics = METRICS) {
    # Generate data for plotting
    plot_data <- data |>
        filter(is_keystone) |>
        select(rel_pop_initial, all_of(metrics)) |>
        pivot_longer(
            cols      = all_of(metrics),
            names_to  = "variable",
            values_to = "value"
        ) |>
        mutate(
            variable      = factor(variable, levels = metrics),
            low_abundance = rel_pop_initial < threshold
        )
    # Generate counts of groups that didnt pass the thershold
    counts <- plot_data |>
        summarise(n = n(), .by = c(variable, low_abundance)) |>
        mutate(color = ifelse(low_abundance, "#d73027", "#4575b4")) |>
        summarise(
            label = paste(
                paste0("<span style='color:", color, "'>&#9679;</span> ", n),
                collapse = "<br>"
            ),
            x_pos = Inf, y_pos = Inf,
            .by = variable
        )
    # Plotting
    plot_data |>
        ggplot(aes(x = rel_pop_initial, y = value, color = low_abundance)) +
        geom_point(alpha = 0.4, size = 0.8) +
        geom_richtext(
            data    = counts,
            aes(x = x_pos, y = y_pos, label = label),
            hjust   = 1.1, vjust = 1.5, size = 3,
            fill    = "white", color = "grey30",
            label.r = unit(0.25, "lines"),
            label.padding = unit(c(0.3, 0.4, 0.3, 0.4), "lines"),
            inherit.aes = FALSE
        ) +
        scale_color_manual(
            values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
            labels = c("TRUE" = paste0("< ", threshold), "FALSE" = paste0("\u2265 ", threshold)),
            name   = "Abundancia relativa"
        ) +
        scale_x_continuous(trans = "log1p") +
        scale_y_continuous(trans = "log1p") +
        facet_wrap(~ variable, scales = "free_y", labeller = labeller(variable = METRIC_LABELS)) +
        labs(title = title, subtitle = subtitle,
            x = "Abundancia relativa (escala log1p)",
            y = "Valor de la métrica (escala log1p)") +
        theme_custom +
        theme(legend.position = "right")
}

# --------------------------------
ctrl_dir = c('/mnt/data/sur/users/mrivera/Data/clean_controls/Negative_magnitude_d1e776c23c74',
            '/mnt/data/sur/users/mrivera/Data/Controls/Cascade_keystone',
            '/mnt/data/sur/users/mrivera/Data/Controls/Boosted_keystone')

save_dir = '/mnt/data/sur/users/mrivera/thesis_plots'
methods <- c("Magnitud negativa", "Efecto cascada", "Boosting de columna")

process_experiment <- function(exp, met, i, save_dir = save_dir, metrics = METRICS) {
    plot_met  <- paste0("ctrl", i)
    sim_params <- data.table::fread(
        file.path(exp, "simulation_params.tsv"),
        select = c("id", "key", "p_noint")
    )
    full_summary <- mclapply(1:nrow(sim_params), process_row, exp = exp, sim_params = sim_params, mc.cores = detectCores() - 1) |>
        do.call(what = rbind)

    save_plot <- function(plot, filename, width = 12, height = 5) {
        path = file.path(save_dir, paste0(plot_met, filename))
        ggsave(path, plot = plot, width = width, height = height)
        cat('>> Image saved at: ', path,'\n')
    }
    variables_to_plot <- c("keystoneness", "dissimilarity", "prop_extinctions")
    # Plot 1: Metric distributions
    make_boxplot(
        data     = full_summary |> select(is_keystone, all_of(variables_to_plot)),
        title    = "Métricas por grupo: todas las especies",
        subtitle = met
    ) |> save_plot("_metrics_distribution.png")

    # Plot 2: Surviving species only
    make_boxplot(
        data     = full_summary |> select(is_keystone, rel_pop_initial, all_of(variables_to_plot)) |> filter(rel_pop_initial > 1e-6),
        title    = "Métricas por grupo: especies sobrevivientes",
        subtitle = met
    ) |> save_plot("_metrics_distribution_surviving.png")

    # Plot 3: Relative abundance vs metrics
    met_relabu(
        data     = full_summary |> select(is_keystone, rel_pop_initial, all_of(variables_to_plot)),
        title    = "Métricas por grupo vs abundancia relativa",
        subtitle = met
    ) |> save_plot("_relabu_vs_metrics.png")
}
i = 3
exp = ctrl_dir[i]
met = methods[i]
process_experiment(exp, met, i, save_dir = save_dir)
for (i in seq_along(ctrl_dir)) {
    exp <- as.character(ctrl_dir[[i]])  # [[ ]] drops names, as.character ensures plain string
    met <- as.character(methods[[i]])
    process_experiment(exp, met, i, save_dir = save_dir)
}


