"""
Elapsed time comparissons plots
=================================================
Purpose:
    Generate model at n epochs to compare running time with normal method.
    Compare if values that model returns are the same

Author: Manuel Rivera
Date:   March 13, 2026
"""
#-------------------------------
# Section: Load data
#-------------------------------
import os 
import numpy as np

dir_path = '/home/mriveraceron/glv-research/Results/KBoost_v2_testing'

# Load DataLoader method data
KEY_MAP = {
    "idx_max_true":  "max_idx_true",
    "idx_max_pred":  "max_idx_pred",
    "metrics_true":  "values_true",
    "metrics_pred":  "values_pred",
    "loss_history":  "loss_history",
    "time_elapsed":  "elpased",
}

def load_results(path):
    data = np.load(path)
    return {var: data[key] for var, key in KEY_MAP.items() if key in data}

# Load both
dl     = load_results(f'{dir_path}/DataLoader_testing/model_results.npz')
normal = load_results(f'{dir_path}/NormalMethod_testing/model_results.npz')
#-------------------------------
from scipy.stats import pearsonr
from scipy.stats import spearmanr

correlationP, pvalue = pearsonr(dl["idx_max_true"], normal["idx_max_true"])
correlationP, pvalue = pearsonr(dl["idx_max_pred"], normal["idx_max_pred"])
correlationP, pvalue = pearsonr(dl["loss_history"], normal["loss_history"])
#-------------------------------
# Section: Loss plot
#-------------------------------
# Section: Plot loss over time
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 5))
plt.plot(loss_history)
plt.title("Loss over time")
plt.xlabel("Epoch")
plt.ylabel("Loss")
plt.grid(True)
plt.savefig(f'{result_path}/DataLoss_plot.png')

#-------------------------------
# Section: Expected vs predicted maximum node
#-------------------------------
import matplotlib.pyplot as plt

# Calculate accuracy
accuracy_DL = np.mean(dl["idx_max_pred"] == dl["idx_max_true"])
accuracy_normal = np.mean(normal["idx_max_pred"] == normal["idx_max_true"])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# --- Left plot ---
accuracy = np.mean(dl["idx_max_true"] == normal["idx_max_true"])
ax1.scatter(x=dl["idx_max_true"], y=normal["idx_max_true"], color='steelblue', s=80)
ax1.text(0.95, 0.95, f'Accuracy: {accuracy:.2f}', transform=ax1.transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax1.set_title("True maximum node index")
ax1.set_xlabel("DataLoader method")
ax1.set_ylabel("Normal method")
ax1.legend()

# --- Right plot ---
accuracy2 = np.mean(dl["idx_max_pred"] == normal["idx_max_pred"])
ax2.scatter(x=dl["idx_max_pred"], y=normal["idx_max_pred"], color='steelblue', s=80)
ax2.text(0.95, 0.95, f'Accuracy: {accuracy2:.2f}', transform=ax2.transAxes,ha='right', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
ax2.set_title("Predicted maximum node index")
ax2.set_xlabel("DataLoader method")
ax2.set_ylabel("Normal method")
ax2.legend()

plt.suptitle("Max Node Index: True vs Predicted", fontsize=14)
plt.tight_layout()
plt.savefig(f"{dir_path}/Indexes_plot.png", dpi=150)

#-------------------------------
# Section: Expected vs predicted values
#-------------------------------
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Generate correlations
correlationP, pvalue = pearsonr(dl["metrics_true"].flatten(), normal["metrics_true"])

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# --- Left plot ---
ax1.scatter(x=dl["metrics_true"], y=normal["metrics_true"], color='steelblue', s=80)
ax1.text(0.95, 0.95, 
    f'Pearson Correlation: {correlationP.item():.4f}',
    transform=ax1.transAxes, ha='right', va='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)
ax1.set_title("True Values")
ax1.set_xlabel("DataLoader method")
ax1.set_ylabel("Normal method")
ax1.legend()

# --- Right plot ---
correlationP, pvalue = pearsonr(dl["metrics_pred"].flatten(), normal["metrics_pred"])
ax2.scatter(x=dl["metrics_pred"], y=normal["metrics_pred"], color='steelblue', s=80)
ax2.text(0.95, 0.95, 
    f'Pearson Correlation: {correlationP.item():.4f}',
    transform=ax2.transAxes, ha='right', va='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)
ax2.set_title("Predicted Values")
ax2.set_xlabel("DataLoader method")
ax2.set_ylabel("Normal method")
ax2.legend()

plt.suptitle("Values: True vs Predicted", fontsize=14)
plt.tight_layout()
plt.savefig(f"{dir_path}/Values_plot.png", dpi=150)


#-------------------------------
# Section: Time plot
#-------------------------------
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8, 5))
epochs = np.arange(1, len(dl["time_elapsed"]) + 1)

ax.plot(epochs, dl["time_elapsed"],     label="DataLoader", marker='o', markersize=3)
ax.plot(epochs, normal["time_elapsed"], label="Normal",     marker='o', markersize=3)
ax.set_title("Elapsed Time per Epoch")
ax.set_xlabel("Epoch")
ax.set_ylabel("Time (s)")
ax.legend()
plt.tight_layout()
plt.savefig(f"{dir_path}/elapsed_time.png", dpi=150)

diff = sum(dl["time_elapsed"]) - sum(normal["time_elapsed"])
print(f'>> The difference of elapsed time: {diff}')