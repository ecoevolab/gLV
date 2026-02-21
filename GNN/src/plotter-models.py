

# Load dataframe
import pathlib
from pathlib import Path
dir = Path('/home/mriveraceron/glv-research/GNN-Logs/24Jan26')

import pandas as pd
df = pd.read_feather(f"{dir}/results_dataframe.feather")   
filtered_df = df[df['tr_NaN'] == False]
index = filtered_df.index

# Load loss
import numpy as np
#paths = [f"{dir}/loss_model-{i}.npy" for i in filtered_df.index]
#loss_list = [np.load(path) for path in paths]

loss_dict = {}

for idx, row in filtered_df.iterrows():
    model = f'{row["arch"]}-L{row["layers"]}-Ch{row["channels"]}'
    loss = np.load(f"{dir}/loss_model-{idx}.npy")
    loss_dict[model] = np.log1p(loss)

# Plot loss over epochs
import matplotlib.pyplot as plt

plt.clf()
for model, loss in loss_dict.items():
    n = 50
    x = np.arange(n)
    plt.plot(x, loss[:n], label=model)

plt.xlabel("Epoch")
plt.ylabel("Log(loss)")
plt.ylabel("Loss per model")
plt.legend()
plt.savefig(f'/home/mriveraceron/glv-research/plots/24Jan26-test.png', dpi=300, bbox_inches='tight')
