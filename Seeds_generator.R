
library(random)
seeds <- randomNumbers(n=2000, min=1, max=10000)

wd <- "~/Documents/LAB_ECO/"
seeds_path <- paste(wd, "./Seeds.tsv", sep="")
file.create(seeds_path) # Create file
write.table(seeds, file = seeds_path, sep = "\t", row.names = FALSE, col.names = FALSE) # Save


