
#` Load parameters table
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/Exp03-D25M02.tsv")

#' The following function is an updated version of the code to generate parameters based on the `parameters_table`.  
#' More information can be found in the `Forge-gLV-Parameters.R` file.  
source("/mnt/atgc-d3/sur/users/mrivera/glv-research/GIT-gLV/Forge-gLV-Parameters.R")

#' The following code regenerates the parameters for each simulation with the function `regenerate`,  
#' calculates the proportions of observed null and negative interactions, returns the results as a matrix, 
#' and converts it to a data frame for column binding.  
system.time({
  tmp <-  lapply(1:nrow(params_table), function(i) {
    
    # cat("Starting simulation ", i,"...\n")
    p <- regenerate(params_table[i,]) # Regenerate the parameters
  
    # Get non-diagonal elements in one step
    non_diag_elements <- p$M[row(p$M) != col(p$M)]
    
    # Efficiently compute proportions using mean (avoids explicit sums)
    p_neg_obs <- mean(non_diag_elements < 0)
    p_nint_obs <- mean(non_diag_elements == 0)
    
    # Return named vector
    return(c(p_neg_obs = p_neg_obs, 
             p_neg_noint = p_nint_obs))
  })
})

res <- unlist(tmp)
tmp_df <- as.data.frame(matrix(res, ncol = 2, byrow = TRUE))
colnames(tmp_df) <- unique(names(res))

#' The following code selects only the column of interest from `params_table`,  
#' adds the counts of true negative and null interactions (`tmp_df`),  
#' and calculates the `Root Mean Squared Error` between Observed and Expected values.
#' 
#' In total, hm_df should have: length(unique(hm_df$p_neg)) * length(unique(hm_df$p_noint))
#' Because the data should be grouped by the combination of parameters.
library(dplyr)
hm_df <- params_table %>%
  select(p_noint, p_neg, n_species) %>% # select columns of interest
  mutate(tmp_df) %>%  # Negative and Null oserved probabilities 
  group_by(p_neg, p_noint, p_pos) %>% 
 summarise(
  Total_specs = sum(n_species), # Calculate total number of species
  Sims_done = n(),  # Count total simulations done with the combination of parameters
  RMSE_neg = sqrt(mean(SE_negs, na.rm = TRUE)), 
  RMSE_noint = sqrt(mean(SE_noint, na.rm = TRUE)),
  RMSE_pos = sqrt(mean(SE_pos, na.rm = TRUE)),
  Total_RMSE = RMSE_neg + RMSE_noint + RMSE_pos,
  .groups = "drop"
)

#' The following code is for creating the heatmap.
library(ggplot2)
library(plotly)

# Create heatmap
heatmap_plot1 <- ggplot(hm_df, aes(x = p_neg, y = p_noint, fill = Total_RMSE,
                                   text = paste("Total Sims:", Sims_done,
                                                "<br>Total_species:", Total_specs,
                                                "<br>p_neg:", p_neg,
                                                "<br>p_noint:", p_noint,
                                                "<br>p_noint:", p_pos,
                                                "<br>RMSE_neg:", RMSE_neg,
                                                "<br>RMSE_noint:", RMSE_noint,
                                                "<br>RMSE_pos:", RMSE_pos
                                                ))) +
  geom_tile(colour = "black") +  # Adjust 'size' to make the line thicker
  scale_fill_gradient(low = "white",  high = "steelblue") +
  labs(x = "p_neg", y = "p_noint", fill = "Total_RMSE") +
  theme_minimal()


# Convert to interactive plots
interactive_heatmap1 <- ggplotly(heatmap_plot1, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))

# Save as an interactive HTML file
library(htmlwidgets)
saveWidget(interactive_heatmap1, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/HM-D25M02-Params.html", selfcontained = FALSE)
