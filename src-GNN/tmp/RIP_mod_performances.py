"""
    This sscript is to generate plots comparing some metrics of different models.
    For each model, the best 2 variants will be chosen and correlations will be plotted
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os


# Your data — replace with actual values
def performance_plot(table, title):
    # --- Data ---
    table = table.dropna().copy()                      # .copy() avoids SettingWithCopyWarning
    table = table[table['pearson_corr'] >= 0].reset_index(drop=True)
    table = table.sort_values(by="pearson_corr", ascending=False)         # sort so best model is at top
    variants = table['model_id']
    x = np.arange(len(variants))
    height = 0.25
    # --- Data-driven bar definitions ---
    metrics = [
        ('pearson_corr',  'Pearson ($r$)',       'steelblue'),
        #('spearman_corr', 'Spearman ($ρ$)',      'tomato'),
        ('accuracy_idx',  'Accuracy ($acc$)',    'mediumseagreen'),
    ]
    offsets = [height, 0, -height]
    fig, ax = plt.subplots(figsize=(8, max(4, len(variants) * 0.6)))  # dynamic height
    ax.set_xlim(0, table[['pearson_corr', #'spearman_corr', 
                          'accuracy_idx']].max().max() + 0.1)
    # Move legend outside the plot to avoid the empty space
    ax.legend(framealpha=0.7, loc='lower right', bbox_to_anchor=(1.0, 0.0))
    for (col, label, color), offset in zip(metrics, offsets):
        bars = ax.barh(x + offset, table[col], height, label=label, color=color, alpha=0.85)
        for bar in bars:
            w = bar.get_width()
            ax.text(w + 0.005, bar.get_y() + bar.get_height() / 2 + offset, f'{w:.2f}', va='center', fontsize=8)
    # --- Labels ---
    ax.set_yticks(x)                                   # removed duplicate set_yticks
    ax.set_yticklabels(variants, fontweight='bold')
    ax.set_xlabel('Value', fontweight='bold')
    ax.set_ylabel('Model variant', fontweight='bold')
    ax.set_title(title, fontweight='bold', pad=12)
    ax.set_xlim(0, 1.15)                               # slightly more room for labels
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.legend(framealpha=0.7, loc='lower right')
    plt.tight_layout()
    return fig  


# Results dir
RES_DIR = '/home/mriveraceron/glv-research/plots/RIP'
os.makedirs(RES_DIR, exist_ok=True)

# Load models tables
DIR = '/home/mriveraceron/glv-research/tuning_results/91074c4e25b4'
graph_table = pd.read_csv(f'{DIR}/tuning_results.csv')
p = performance_plot(graph_table, title = 'GraphConv performance metrics by variant')
plt.savefig(f'{RES_DIR}/GraphConv_performance.jpg', dpi=300, bbox_inches='tight')

# SAGE model
DIR = '/home/mriveraceron/glv-research/tuning_results/SAGE-hpo'
sage_table = pd.read_csv(f'{DIR}/tuning_results.csv')
p = performance_plot(sage_table, title = 'SAGE performance metrics by variant')
plt.savefig(f'{RES_DIR}/SAGE_performance.jpg', dpi=300, bbox_inches='tight')

# Random forest model
DIR = '/home/mriveraceron/glv-research/Results/RForest'
random_table = pd.read_csv(f'{DIR}/rf_results.csv')
random_table = random_table.rename(columns={'node_acc': 'accuracy_idx', 'name':'model_id'})
p = performance_plot(random_table, title = 'RF performance metrics by variant')
plt.savefig(f'{RES_DIR}/RF_performance.jpg', dpi=300, bbox_inches='tight')


