
#`-----------------------------Load master table------------#'
params_table <- data.table::fread("/mnt/atgc-d3/sur/users/mrivera/glv-research/Data/D13M02Y25.tsv")

#' Function to regenerate parameters.
regenerate <- function(index) {
  
  N_species <- as.numeric(index[["n_species"]])
  
  #------------------Populations-----------------------------#
  set.seed(as.numeric(index[["Pop_seed"]]))
  Pobl <- stats::runif(N_species, min = 0.1, max = 1)
  
  #------------------------Growth Rates---------------------#
  set.seed(as.numeric(index[["Growth_seed"]]))
  Grow <- stats::runif(N_species, min = 0.001, max = 1)
  
  #--------------------Interactions-------------------------#c
  set.seed(as.numeric(index[["A_seed"]]))
  
  # Probability of negative interaction and vector of 0's and 1's
  p_neg <- as.numeric(index[["p_neg"]])
  V_neg <- stats::rbinom(N_species * N_species, 1, p_neg)
  
  # Probability of null interaction and vector of 0's and 1's
  p_noint <- as.numeric(index[["p_noint"]])
  V_noint <- stats::rbinom(N_species * N_species, 1, 1 - p_noint)
  
  tmp <- V_noint * ifelse(V_neg != 0,
                          stats::runif(N_species * N_species, min = 0, max = 1),
                          -stats::runif(N_species * N_species, min = 0, max = 1)
  )
  
  inter <- matrix(tmp, nrow = N_species, ncol = N_species)
  
  # Set diagonal values
  diag(inter) <- -0.5
  
  # Extract ID
  id <- index[["id"]]
  
  # Return parameters as a list
  params <- list(x0 = Pobl,
                 M = inter,
                 mu = Grow,
                 ID = id)
  
  return(params)
}


#' The following code regenerates the parameters for each simulation with the function `regenerate`,  
#' calculates the proportions of observed null and negative interactions, returns the results as a matrix, 
#' and converts it to a data frame for column binding.  
tmp <- apply(params_table, 1, function(x) {
  params <- regenerate(x) # Regenerate the parameters
  true_negs <- (sum(params$M < 0))/ncol(params$M)^2 # Proportions of observed negative ineractions
  true_noint <- (sum(params$M == 0))/ncol(params$M)^2 # Proportions of observed null ineractions
  return(c(true_noint = true_noint, true_negs = true_negs ))
}) 

# Convert to data frame
tmp_df <- as.data.frame(t(tmp))


#' The following code selects only the column of interest from `params_table`,  
#' adds the counts of true negative and null interactions (`tmp_df`),  
#' and calculates the `Row-wise Root Mean Squared Error` between Observed and Expected values and the #' average of both parameters 
library(dplyr)
hm_df <- params_table %>%
  select(p_noint, p_neg) %>% # select columns of interest
  mutate( # Add columns
    tmp_df, # Add counts
    RMSE_neg = sqrt((true_negs - p_neg)^2) ,  # Compute RMSE for negatives
    RMSE_noint = sqrt((true_noint - p_noint)^2) 
  )  %>% 
  group_by(p_neg, p_noint) %>% 
summarise(
  Total_Sims = n(),  # Count total simulations
  Sum_RMSE_neg = sum(RMSE_neg, na.rm = TRUE), 
  Sum_RMSE_noint = sum(RMSE_noint, na.rm = TRUE),  
  .groups = "drop"
) %>%

head(hm_df)

  

#' The following code is for creating the heatmap.
library(ggplot2)
library(plotly)

# Create heatmap
heatmap_plot1 <- ggplot(hm_df, aes(x = p_neg, y = p_noint, fill = true_negs)) +
  geom_tile(colour = "black") +  # Adjust 'size' to make the line thicker
  scale_fill_gradient(low = "white",  high = "steelblue") +
  labs(x = "p_neg", y = "p_noint", fill = "true_negs") +
  theme_minimal()
# Create heatmap
heatmap_plot2 <- ggplot(hm_df, aes(x = p_neg, y = p_noint, fill = true_noint)) +
  geom_tile(colour = "black") +  # Adjust 'size' to make the line thicker
  scale_fill_gradient(low = "white",  high = "tomato") +
  labs(x = "p_neg", y = "p_noint", fill = "true_noint") +
  theme_minimal()

# Convert to interactive plots
interactive_heatmap1 <- ggplotly(heatmap_plot1, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))
interactive_heatmap2 <- ggplotly(heatmap_plot2, tooltip = "text") %>% layout(hoverlabel = list(bgcolor = "white"))

#-----------------------------------#
# Combine both interactive plots in a single figure
combined_plot <- subplot(interactive_heatmap1, interactive_heatmap2, nrows = 1, shareY = TRUE, titleX = TRUE) %>%
  layout(
    annotations = list(
      list(x = 0.2, y = 1, text = "True negative proportion", showarrow = FALSE, xref='paper', yref='paper', font=list(size=16)),
      list(x = 0.8, y = 1, text = "True null  proportion", showarrow = FALSE, xref='paper', yref='paper', font=list(size=16))
    ),
    xaxis = list(title = "Negative Probability"),
    xaxis2 = list(title = "Negative Probability"),
    yaxis = list(title = "No Interaction Probability")
  )

# Save as an interactive HTML file
library(htmlwidgets)
saveWidget(combined_plot, "/mnt/atgc-d3/sur/users/mrivera/glv-research/Graphs/heatmap_plot.html", selfcontained = FALSE)
