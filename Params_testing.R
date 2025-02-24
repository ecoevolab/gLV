
#`-----------------------------Load master table------------#'
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv")

# params_table <- data.table::fread("Cluster/glv-research/Data/D13M02Y25.tsv")


index <- params_table[1,]

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
  
  # Create a matrix full of NAs
  M <- rep(NA, times = n_species^2)
  
  # Assign diagonal elements as -0.5
  diagonal <- seq(from = 1, to = n_species^2, by = n_species + 1)
  M[diagonal] <- -0.5
  
  # Define proportions for zero and negative values
  p_noint <- as.numeric(index[["p_noint"]])
  p_neg <- as.numeric(index[["p_neg"]])
  set.seed(as.numeric(index[["A_seed"]]))
  
  # Calculate interaction counts
  #' For testing run the next: num_noint + num_negs + num_pos
  interactions_index <- which(is.na(M))  # Indices of off-diagonal elements
  num_noint <- floor(p_noint * (n_species^2 - n_species))  # Neutral interactions
  num_negs <- floor(p_neg * ((n_species^2 - n_species) - num_noint))  # Negative interactions
  num_pos <- (n_species^2 - n_species) - (num_noint + num_negs)  # Positive interactions
  
  # Assign neutral interactions (set to zero)
  #' sort(sample_noint)
  #' sort(interactions_index)
  sample_noint <- sample(interactions_index, size = num_noint, replace = FALSE)
  M[sample_noint] <- 0
  interactions_index <- interactions_index[!interactions_index %in% sample_noint]  # Remove from pool
  
  # Assign negative interactions (set to random negative values)
  if (num_negs > 1) {
    sample_negs <- sample(interactions_index, size = num_negs, replace = FALSE)
    M[sample_negs] <- -runif(n = num_negs, min = 0, max = 1)
    interactions_index <- interactions_index[!interactions_index %in% sample_negs]  # Remove from pool
  } 
  
  if (num_negs == 1) {
    M[interactions_index] <- -runif(n = num_negs, min = 0, max = 1)
    interactions_index <- c()
  }
  
  # Assign positive interactions to remaining positions
  if (length(interactions_index) != 0) {
    M[interactions_index] <- runif(n = num_pos, min = 0, max = 1)
  } 
  
  # Reshape the vector into a matrix and round it
  M_df <- matrix(round(M, digits = 5), nrow = n_species, ncol = n_species)
  
  # Extract ID
  id <- index[["id"]]
  
  # Return parameters as a list
  params <- list(x0 = x0,
                 M = M_df,
                 mu = mu,
                 id = id,
                 n = n_species)
  
  return(params)
}

#' The following code regenerates the parameters for each simulation with the function `regenerate`,  
#' calculates the proportions of observed null and negative interactions, returns the results as a matrix, 
#' and converts it to a data frame for column binding.  
tmp <-  lapply(1:nrow(params_table), function(i) {
  
  cat("Starting simulation ", i,"...\n")
  p <- regenerate(params_table[i,]) # Regenerate the parameters

  # Get non-diagonal elements in one step
  non_diag_elements <- p$M[row(p$M) != col(p$M)]
  
  # Efficiently compute proportions using mean (avoids explicit sums)
  true_pos <- mean(non_diag_elements > 0)
  true_neg <- mean(non_diag_elements < 0)
  true_noint <- mean(non_diag_elements == 0)
  
  # Return named vector
  return(c(PropsPos = true_pos, PropsNeg = true_neg, PropsZer = true_noint))
})

res <- unlist(tmp)
tmp_df <- as.data.frame(matrix(res, ncol = 3, byrow = TRUE))
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
  mutate( # Add columns
    tmp_df, # Add counts
    SE_negs = (PropsNeg - p_neg)^2,
    SE_noint = (PropsZer - p_noint)^2
  )  %>% 
  group_by(p_neg, p_noint) %>% 
 summarise(
  Total_specs = sum(n_species), # Calculate total number of species
  Sims_done = n(),  # Count total simulations done with the combination of parameters
  RMSE_neg = sqrt(mean(SE_negs, na.rm = TRUE)), 
  RMSE_noint = sqrt(mean(SE_noint, na.rm = TRUE)),
  RMSE_total = sqrt(mean(RMSE_noint, RMSE_neg)),
  .groups = "drop"
)

#' The following code is for creating the heatmap.
library(ggplot2)
library(plotly)

# Create heatmap
heatmap_plot1 <- ggplot(hm_df, aes(x = p_neg, y = p_noint, fill = RMSE_total,
                                   text = paste("Total Sims:", Sims_done,
                                                "<br>Total_species:", Total_specs,
                                                "<br>p_neg:", p_neg,
                                                "<br>p_noint:", p_noint,
                                                "<br>RMS_neg:", RMSE_neg,
                                                "<br>RMSE_noint:", RMSE_noint
                                                ))) +
  geom_tile(colour = "black") +  # Adjust 'size' to make the line thicker
  scale_fill_gradient(low = "white",  high = "steelblue") +
  labs(x = "p_neg", y = "p_noint", fill = "RMSE_total") +
  theme_minimal()


# Convert to interactive plots
interactive_heatmap1 <- ggplotly(heatmap_plot1, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))

# Save as an interactive HTML file
library(htmlwidgets)
saveWidget(interactive_heatmap1, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/New-D23M02-updated01.html", selfcontained = FALSE)
