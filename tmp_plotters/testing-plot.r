


info <- data.table::fread("/home/rivera/cluster/glv-research/Results/TEST01-D21-May/sims-summary.tsv")
data = data.table::fread("/home/rivera/cluster/glv-research/Data/TEST01-D21-May.tsv")

library(dplyr)

data_to_add = expand.grid(p_neg = 0, p_noint = 0, n_failed = 10, 
              total = 10, n_ran = 0, n_species = seq(20, 100, 20),
              prop = 0)

df = data  %>% 
      full_join(info, by = "id") %>% 
      group_by(p_neg, p_noint, n_species) %>% 
      summarise(n_failed = sum(na_ct, na.rm = TRUE), total = n(), .groups = "drop") %>% 
      mutate(n_ran = total - n_failed) %>%
      mutate(prop = n_ran/total) %>%
      dplyr::full_join(data_to_add, by = c("p_neg", "p_noint", "n_species", "prop", "n_failed", "total", "n_ran"))


library(ggplot2)
library(svglite)
library(pheatmap)
library(reshape2)


theme_blackbox <- function(border.size = 3){
  theme(axis.text = element_text(color="black"),
        axis.title = element_text(face="bold"),
        panel.background = element_rect(color = "black",
                                        size = border.size,
                                        fill = NA),
        panel.grid = element_blank())
}

p= ggplot(df, aes(x = p_noint, y = p_neg)) + 
  facet_wrap(~ n_species) + 
  scale_fill_gradientn(colours = c("#fee8c8", "#e34a33"), values = c(0,1)) +
  theme_blackbox() +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  geom_tile( aes(fill = prop), color = "black") 

# Add the highlighted row
p +  geom_rect(aes(xmin = -0.05, xmax = 1.05, ymin = 0.95, ymax = 1.05), fill = NA, color = "#00bfc4", size = 1) 

p2 = p + 
  annotate("segment", x = 0, xend = 0, y = 0, yend = 1, colour = "#0571b0", arrow = arrow(length = unit(0.3, "cm"))) +   # azul oscuro
  annotate("segment", x = 0, xend = 1, y = 0, yend = 0, colour = "#800080", arrow = arrow(length = unit(0.3, "cm")))     # morado intenso


ggsave(
  filename = "/home/rivera/thesis/na-count.pdf",
  plot = p2,
  width = 8,
  height = 6,
  units = "in"
)

