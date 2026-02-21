import pandas as pd
import time

# Load data
data = pd.read_csv("/mnt/data/sur/users/mrivera/Train-sims/4379fd40-9f0a/parameters-sims.tsv", sep='\t')
data_indexed = data.set_index('id')

# Test with 1000 lookups
spec_ids = data['id'].head(1000).tolist()

# Method 1: WITHOUT indexing (filtering each time)
start = time.time()
for spec_id in spec_ids:
    nsp = data[data["id"] == spec_id]['n_species'].values[0]
without_index_time = time.time() - start

# Method 2: WITH indexing
start = time.time()
for spec_id in spec_ids:
    nsp = data_indexed.loc[spec_id, 'n_species']
with_index_time = time.time() - start

print(f"Without index: {without_index_time:.3f}s")
print(f"With index: {with_index_time:.3f}s")
print(f"Speedup: {without_index_time/with_index_time:.1f}x faster")