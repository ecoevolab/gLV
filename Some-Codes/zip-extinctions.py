# 2025-10-02: Extinctions File Compression 
# Zips all simulation extinctions.
import pandas as pd
import os
from zipfile import ZipFile, ZIP_DEFLATED
from multiprocessing import Pool
import multiprocessing

def zip_worker(args):
    # Declare arguments
    spec_id, nsp, ext_dir, zip_dir = args
    # Extinctins of simulation spec_id
    src_files = [f'{ext_dir}E_{spec_id}-S{i}.feather' for i in range(1, nsp + 1)]
    zip_path = f"{zip_dir}/exts-{spec_id}.zip"
    # Generate Zip
    with ZipFile(zip_path, 'w', ZIP_DEFLATED, compresslevel=6) as zipf:
        for file in src_files:
            if os.path.exists(file):
                zipf.write(file, arcname=os.path.basename(file))
    print(f">> Finished zipping of extinctions for: {spec_id}")           
    return True

# Setup
ext_dir = "/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/Post-exts/"
zip_dir = "/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/PostExts-zips"
os.makedirs(zip_dir, exist_ok=True)

# Load data
data = pd.read_csv("/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/parameters-sims.tsv", sep='\t')
data_indexed = data.set_index('id')

# Prepare work
args_list = [
    (spec_id, data_indexed.loc[spec_id, 'n_species'], ext_dir, zip_dir)
    for spec_id in data_indexed.index
]


num_cores = multiprocessing.cpu_count()
with Pool(processes=num_cores-2) as pool:
    results = pool.map(zip_worker, args_list)

