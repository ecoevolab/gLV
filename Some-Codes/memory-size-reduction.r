# 2025-10-02: File Compression Analysis
# Evaluates storage reduction by zipping simulation files
# Test case: 100 species specs (worst-case: 100 extinctions)
# Calculates: uncompressed size → compressed size → % improvement
library(data.table)
dat = readr::read_tsv("/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/parameters-sims.tsv")
library(dplyr)
specs100 = dat %>% filter(n_species == 100) %>% slice(1)

pattern = paste0("E_", specs100$id, "-*")
ext_dir = "/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/Post-exts"
mcmd <- sprintf('find "%s" -type f -name "%s" ', ext_dir, pattern)               # Command
src_files = system(mcmd, intern = TRUE)   

unzip_mb = sum(file.info(src_files)[, "size"]/(1024^3))   # size in GB (3)

zip("/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/100specs-test.zip", src_files)
zip_mb = sum(file.info("/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/100specs-test.zip")[, "size"]/(1024^3))
performance = 1 - (zip_mb/unzip_mb)
cat(sprintf(" >> Memory reduced by %.1f%% (%.2f GB → %.2f GB)\n", performance * 100, unzip_mb, zip_mb))

# Directory reduction
total_unzip = unzip_mb * nrow(dat)
total_zip = zip_mb * nrow(dat)
total_performance = 1- (total_zip/total_unzip)
cat(sprintf(" >> Total memory reduced by %.1f%% (%.2f GB → %.2f GB)\n", total_performance * 100, total_unzip, total_zip))