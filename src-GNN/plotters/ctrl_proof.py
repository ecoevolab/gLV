"""
    Control plots

    Description: 
    Script to generate plots of our control method for testing wether a GNN was capable of predicting keystoneness.
    This approach was just for a proof of concept.
"""
# Imports
import numpy as np

data_dir = '/home/mriveraceron/glv-research/Results/Boosted_keystone/Filtered_AllFeats_V3'
loss = np.load(f'{data_dir}/loss_history.npy')
mp = np.load(f'{data_dir}/values_pred.npy')
mt = np.load(f'{data_dir}/values_true.npy')
idxt = np.load(f'{data_dir}/max_idx_true.npy')
idxp =  np.load(f'{data_dir}/max_idx_pred.npy')
np.savez(f'{data_dir}/metric_values.npz',
    idxt  = idxt,
    idxp  = idxp,
    mt   = mt,
    mp   = mp,
    loss  = loss
)

# Load data
data = np.load(f'{data_dir}/metric_values.npz')


# Plotting loss


# Plotting predictions


# Confussion matrix
