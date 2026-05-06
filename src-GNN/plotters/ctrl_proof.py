"""
    Control plots

    Description: 
    Script to generate plots of our control method for testing wether a GNN was capable of predicting keystoneness.
    This approach was just for a proof of concept.
"""
# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Load data
data_dir = '/home/mriveraceron/glv-research/Results/Boosted_keystone/Filtered_AllFeats_V3'
data = np.load(f'{data_dir}/metric_values.npz')
save_dir = '/home/mriveraceron/glv-research/thesis_plots'


# Plotting settings
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'DejaVu Serif'],
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 15,
    'axes.linewidth': 1.2,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'legend.fontsize': 11,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8',
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

#------------------- Plotting loss ------------------- 
loss_train = data['loss']
plt.clf
# Number of epochs
epochs = range(1, len(loss_train) + 1)
fig, ax = plt.subplots(figsize=(8, 8))
ax.plot(epochs, loss_train, color='#2166ac',linewidth=1.8,label='Pérdida de entrenamiento')
# Axes labels & title
ax.set_xlabel('Época', labelpad=8)
ax.set_ylabel('Pérdida', labelpad=8)
# Grid — subtle, behind data
ax.grid(True, which='major', linestyle='--', linewidth=0.6, alpha=0.5)
ax.grid(True, which='minor', linestyle=':', linewidth=0.4, alpha=0.3)
ax.set_axisbelow(True)
# Tick formatting
ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
# Legend
ax.legend(loc='upper right')
# Clean spines (keep bottom & left only)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_title('Evolución de la pérdida durante el entrenamiento', pad=12)
ax.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig(f'{save_dir}/ctrl_proof_loss.png', dpi=150)

# ----------- Values plot -----------
# Plotting predictions
mt = data['mt']
mp = data['mp']
# Compute correlations
r_p, _ = pearsonr(mt, mp)
r_s, _ = spearmanr(mt, mp)
x = np.log1p(mp)
y = np.log1p(mt)

# Hexbin
plt.close('all')
fig, ax = plt.subplots(figsize=(5, 5))
hb = ax.hexbin(
    x, y,
    gridsize=45,
    cmap='YlOrRd',
    mincnt=1,
    bins='log',
    linewidths=0.15,
    edgecolors='none',      # cleaner look at high density
)

# Colorbar
cb = fig.colorbar(hb, ax=ax, pad=0.03, shrink=0.82, aspect=22)
cb.set_label('Recuento (escala log)', fontsize=11, labelpad=8)
cb.ax.tick_params(labelsize=10, direction='in')
cb.outline.set_linewidth(0.8)

# Identity line (y = x)
lims = [-0.02, 1.05]
ax.plot(lims, lims, color='0.3', linewidth=1.2, linestyle='--', alpha=0.85, zorder=5, label='Correlación perfecta')

# Labels, title
ax.set_xlabel('Keystoneness predicho', labelpad=8)
ax.set_ylabel('Keystoneness esperado', labelpad=8)
fig.suptitle('GraphConv — valores de keystoneness', fontsize=15, fontfamily='serif')
ax.set_title('Prueba de concepto', fontsize=11, color='0.4', pad=6)
# Figure limits
ax.set_xlim(lims)
ax.set_ylim(lims)
# Ticks
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

# Correlation annotation box
stats_text = (
    f'$r$ = {r_p:.3f}$\;$(Pearson)\n'
    f'$\\rho$ = {r_s:.3f}$\;$(Spearman)'
)
ax.text(
    0.05, 0.95, stats_text,
    transform=ax.transAxes,
    ha='left', va='top',
    fontsize=10,
    fontfamily='serif',
    bbox=dict(
        boxstyle='round,pad=0.5,rounding_size=0.4',
        facecolor='white',
        edgecolor='0.75',
        linewidth=0.8,
        alpha=0.9,
    ),
    zorder=10,
)

# Legend for identity line
ax.legend(loc='lower right', fontsize=10, framealpha=0.9, edgecolor='0.75')

# Clean spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(f'{save_dir}/ctrl_proof_values.png')

# ----------- Confussion matrix -----------
