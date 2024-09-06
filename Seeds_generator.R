
library(random)
seeds <- randomNumbers(n=2000, min=1, max=10000)

seeds_path <- "./Seeds.tsv"
file.create(seeds_path) # Create file
write.table(seeds, file = seeds_path, sep = "\t", row.names = FALSE, col.names = FALSE) # Save


