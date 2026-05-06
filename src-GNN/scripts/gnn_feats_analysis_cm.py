"""
Confussion matrix for feature analysis.
---------------------------------------

GraphConv best model was tested with or without features.
Confussion matrix will be generated based on its results.
"""

# Imports
import torch
import numpy as np
import pandas as pd

# Experiment directory
exp_dir = '/home/mriveraceron/glv-research/Results/feats_analysis'
models = {'features_run', 'dummy_run'}

records = []
for mod in models:
    data = np.load(f'{exp_dir}/{mod}/model_preds.npz')
    idxt        = data['idxt_eval']
    idxp        = data['idxp_eval']
    graph_nodes = data['nodes_eval']
    # Confusion matrix
    tp = sum((idxt == idxp).astype(int))
    fp = sum((idxt != idxp).astype(int))
    fn = fp
    tn = sum((graph_nodes - 1)) - fp
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    ppv = tp / (tp + fp)
    print(f'Model {mod} had an accuracy of {accuracy} and precision of {ppv}')
    records.append({'model': mod, 'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn})

df_cm = pd.DataFrame(
    [[tp, fp], [fn, tn]],
    index   = ['Expected_Positive',  'Expected_Negative'],
    columns = ['Predicted_Positive','Predicted_Negative']
)