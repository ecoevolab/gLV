"""
Hexbin plot for values
22-March-2026
Description:
    Script to generate plots for predicted and expected values, while avoiding overplotting some datapoints.

"""

# Load data
import os
import numpy as np

dir = '/home/mriveraceron/glv-research/tuning_results/91074c4e25b4/Variant_11'
data = np.load(f'{dir}/Variant_11-values.npz')

mt  = data['values_true']
mp  = data['values_pred']

# Commpute correlations
from scipy.stats import pearsonr
from scipy.stats import spearmanr

correlationP, _ = pearsonr(mt.flatten(), mp.flatten())
correlationS, _ = spearmanr(mt.flatten(), mp.flatten())

#-------------------
# Hexbin plot
#-------------------
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
hb = ax.hexbin(x=mp.flatten(), y=mt.flatten(), gridsize=30, cmap="Blues")
plt.colorbar(hb, ax=ax, label="Count")

ax.plot([0, 1], [0, 1], 'r--', linewidth=1.5, label="Ideal")
ax.set_xlabel("Predicted values")
ax.set_ylabel("Expected values")
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_title("Keystoneness values")
ax.legend(loc="lower right")
ax.text(0.95, 0.95,
    f'Pearson Correlation: {correlationP:.4f}\nSpearman Correlation: {correlationS:.4f}',
    transform=ax.transAxes, ha='right', va='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)
ax.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(f'{dir}/Hexbin_plot.png', dpi=150)

#-------------------
# Density plot
#-------------------
import seaborn as sns
import matplotlib.pyplot as plt

g = sns.jointplot(x=mp.flatten(), y=mt.flatten(), kind="kde", fill=True, cmap="Blues")

g.set_axis_labels("Predicted values", "Expected values")
g.figure.suptitle("Keystoneness values", y=1.02)

# Add correlations to the main axis
g.ax_joint.text(0.95, 0.95,
    f'Pearson Correlation: {correlationP:.4f}\nSpearman Correlation: {correlationS:.4f}',
    transform=g.ax_joint.transAxes, ha='right', va='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)

# Diagonal reference line
g.ax_joint.plot([0, 1], [0, 1], 'r--', linewidth=1.5, label="Ideal")
g.ax_joint.set_xlim(0, 1)
g.ax_joint.set_ylim(0, 1)
g.ax_joint.legend(loc="lower right")
g.ax_joint.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig(f'{dir}/KDE_plot.png', dpi=150)