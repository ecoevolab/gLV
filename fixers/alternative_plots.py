import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr

# Generate values plot
dir = '/home/mriveraceron/glv-research/tuning_results/GraphConv_sample_size'
path = f'{dir}/Variant_5/metric-values.npz'
data = np.load(path)

# Assign vectors
y_true = data['values_true']
y_pred = data['values_pred']

# Generate correlations
correlationP, _ = pearsonr(y_pred, y_true)
correlationS, _ = spearmanr(y_pred, y_true)

eps = 1
plt.clf() # 
fig, ax = plt.subplots()
hb = ax.hexbin(
    np.log10(y_true + eps), np.log10(y_pred + eps),
    gridsize=50,
    cmap="Blues",
    mincnt=1,
    bins="log"  # log color scale so dense cluster doesn't wash out the rest
)
plt.colorbar(hb, ax=ax, label="log10(count)")

# Regression line
m, b = np.polyfit(np.log10(y_true.flatten() + eps), np.log10(y_pred.flatten() + eps), 1)
x = np.linspace(np.log10(y_true.min() + eps), np.log10(y_true.max() + eps), 100)
ax.plot(x, m * x + b, color="red", linewidth=2)
ax.text(
    0.98, 0.98,
    f"Pearson: {correlationP.item():.4f}\nSpearman: {correlationS.item():.4f}",
    transform=ax.transAxes,
    fontsize=10,
    verticalalignment="top",
    horizontalalignment="right",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
)
ax.set_xlabel("log10(True values + 1)")
ax.set_ylabel("log10(Predicted values + 1)")
ax.set_title("Keystoneness Values")
plt.savefig(f'{dir}/Variant_5/Hexbin_plot.png', dpi=150, bbox_inches="tight")
plt.show()