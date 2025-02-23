
#`-----------------------------Load master table------------#'
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv")

# params_table <- data.table::fread("Cluster/glv-research/Data/D13M02Y25.tsv")


index <- params_table[1,]
M <- matrix(NA, nrow = 5, ncol = 5)

#' The following function is an updated version of the code to generate parameters based on the `parameters_table`.  
#' More information can be found in the `Forge-gLV-Parameters.R` file.  

regenerate <- function(index) {
  
  n_species <- as.numeric(index[["n_species"]])
  
  #------------------Populations-----------------------------#
  set.seed(as.numeric(index[["Pop_seed"]]))
  x0 <- stats::runif(n_species, min = 0.1, max = 1)
  
  #------------------------Growth Rates---------------------#
  set.seed(as.numeric(index[["Growth_seed"]]))
  mu <- stats::runif(n_species, min = 0.001, max = 1)
  
  #--------------------Interactions-------------------------#
  
  #' Create a matrix full of ones and diagonal of 0
  M <- matrix(1, nrow = n_species, ncol = n_species)
  diag(M) <- -0.5
  
  #' Define proportions for zero and negative values
  p_noint <- as.numeric(index[["p_noint"]])
  p_neg <- as.numeric(index[["p_neg"]])
  set.seed(as.numeric(index[["A_seed"]]))
  
  #' Calculate off-diagonal indices and interaction counts
  indices_off_diag <- which(row(M) != col(M), arr.ind = TRUE)
  Noff_diag <- nrow(indices_off_diag)
  num_0 <- floor(p_noint * Noff_diag)
  num_negs <- floor(p_neg * (Noff_diag - num_0))
  num_pos <- Noff_diag - (num_0 + num_negs)
  
  # Reorder the off-diagonal `indices` randomly
  zero_index <- sample(Noff_diag, num_0)
  
  # Assign null interactions (set to zero)
  M[indices_off_diag[zero_index, 1], indices_off_diag[zero_index, 2]] <- 0
  
  
  # Assign negative interactions
  M[indices_off_diag[reorder_index[(num_0 + 1):(num_0 + num_negs)], ]] <- -runif(n = num_negs, min = 0, max = 1) 
  
  # Assign positive interactions
  M[indices_off_diag[reorder_index[(num_0 + num_negs + 1):Noff_diag], ]] <- runif(n = num_pos, min = 0, max = 1) 
  
  # Extract ID
  id <- index[["id"]]
  
  # Return parameters as a list
  params <- list(x0 = x0,
                 M = M,
                 mu = mu,
                 id = id,
                 n = n_species)
  
  return(params)
}

#' These lines are for testing purposes, dont need to be runned. 
x <- params_table[1,]
p <- regenerate(x)
# 
# non_diag_indices <- which(row(p$M) != col(p$M), arr.ind = TRUE)
# test_pos <- sum(p$M[non_diag_indices] > 0)/(p$n^2-p$n)
# test_neg <- sum(p$M[non_diag_indices] < 0)/(p$n^2-p$n)
# test_noint <- sum(p$M[non_diag_indices] == 0)/(p$n^2-p$n)

#' The following code regenerates the parameters for each simulation with the function `regenerate`,  
#' calculates the proportions of observed null and negative interactions, returns the results as a matrix, 
#' and converts it to a data frame for column binding.  
tmp <-  lapply(1:nrow(params_table), function(i) {
  
  p <- regenerate(i) # Regenerate the parameters
  # p$M = 2 - diag(3)

  # Get non-diagonal elements in one step
  non_diag_elements <- p$M[row(p$M) != col(p$M)]
  
  # Efficiently compute proportions using mean (avoids explicit sums)
  true_pos <- mean(non_diag_elements > 0)
  true_neg <- mean(non_diag_elements < 0)
  true_noint <- mean(non_diag_elements == 0)
  
  # Return named vector
  return(c(PropsPos = true_pos, PropsNeg = true_neg, PropsZer = true_noint))
}) 

tmp <- lapply(1:nrow(params_table), function(i) {
  print(params_table[i, ])  # Check the row
  
  p <- regenerate(params_table[i, ])  # Pass the row to regenerate
  print(p)  # Inspect the returned value
  
  # Additional debugging can be added inside regenerate()
})

























# Convert to data frame
tmp_df <- as.data.frame(t(tmp))

# 
#' The following code selects only the column of interest from `params_table`,  
#' adds the counts of true negative and null interactions (`tmp_df`),  
#' and calculates the `Root Mean Squared Error` between Observed and Expected values.
#' 
#' In total, hm_df should have: length(unique(hm_df$p_neg)) * length(unique(hm_df$p_noint))
#' Because the data should be grouped by the combination of parameters.
library(dplyr)
hm_df <- params_table %>%
  select(p_noint, p_neg) %>% # select columns of interest
  mutate( # Add columns
    tmp_df, # Add counts
    SE_neg = (true_negs - p_neg)^2 ,  # Compute Squared Error for negatives
    SE_noint = (true_noint - p_noint)^2 # Compute Squared Error for no interactions
  )  %>% 
  group_by(p_neg, p_noint) %>% 
 summarise(
  Sims_done = n(),  # Count total simulations done with the combination of parameters
  RMSE_neg = sqrt(mean(SE_neg, na.rm = TRUE)), 
  RMSE_noint = sqrt(mean(SE_noint, na.rm = TRUE)),
  # RMSE_combined = sqrt(sum(RMSE_neg, RMSE_noint)/2),
  .groups = "drop"
)

#' The following code is for creating the heatmap.
library(ggplot2)
library(plotly)

# Create heatmap
heatmap_plot1 <- ggplot(hm_df, aes(x = p_neg, y = p_noint, fill = RMSE_neg,
                                   text = paste("Total Sims:", Sims_done,
                                                "<br>p_neg:", p_neg,
                                                "<br>p_noint:", p_noint,
                                                "<br>RMSE:", RMSE_neg,
                                                ))) +
  geom_tile(colour = "black") +  # Adjust 'size' to make the line thicker
  scale_fill_gradient(low = "white",  high = "steelblue") +
  labs(x = "p_neg", y = "p_noint", fill = "RMSE_neg") +
  theme_minimal()
# Create heatmap
heatmap_plot2 <- ggplot(hm_df, aes(x = p_neg, y = p_noint, fill = RMSE_noint,
                                   text = paste("Total Sims:", Sims_done,
                                                "<br>p_neg:", p_neg,
                                                "<br>p_noint:", p_noint,
                                                "<br>RMSE:", RMSE_noint,
                                   ))) +
  geom_tile(colour = "black") +  # Adjust 'size' to make the line thicker
  scale_fill_gradient(low = "white",  high = "tomato") +
  labs(x = "p_neg", y = "p_noint", fill = "RMSE_noint") +
  theme_minimal()

# Convert to interactive plots
interactive_heatmap1 <- ggplotly(heatmap_plot1, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
interactive_heatmap2 <- ggplotly(heatmap_plot2, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))

#-----------------------------------#
# Combine both interactive plots in a single figure
combined_plot <- subplot(interactive_heatmap1, interactive_heatmap2, nrows = 1, shareY = TRUE, titleX = TRUE) %>%
  layout(
    annotations = list(
      list(x = 0.2, y = 1, text = "RMSE negative int", showarrow = FALSE, xref='paper', yref='paper', font=list(size=16)),
      list(x = 0.8, y = 1, text = "RMSE null int", showarrow = FALSE, xref='paper', yref='paper', font=list(size=16))
    ),
    xaxis = list(title = "Negative Probability"),
    xaxis2 = list(title = "Negative Probability"),
    yaxis = list(title = "No Interaction Probability")
  )

# Save as an interactive HTML file
library(htmlwidgets)
saveWidget(combined_plot, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/heatmap_plot.html", selfcontained = FALSE)
