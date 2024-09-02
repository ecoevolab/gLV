# Load the data from the TSV file
data <- read.table("/home/rivera/Cluster/Scan/CPr_time.tsv", header = TRUE, sep = "\t")

# Order the data based on a specific column, for example, the first column
ordered_data <- data[order(data$Simulation_time.s), ]
