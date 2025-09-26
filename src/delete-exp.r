


id = "aaf5fb40-8dc1"

lapply(id, function(id){
    # Declare paths
    dat = file.path("/mnt/data/sur/users/mrivera/Data", paste0(id,".tsv"))
    res = file.path("/mnt/data/sur/users/mrivera/Experiments", id)
    # Remove files
    unlink(res, recursive = TRUE)
    print("Experiment files removed")
    unlink(dat)
    print("Data files removed")
})