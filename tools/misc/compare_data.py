"""
Author: Gokhan Oztarhan
Created date: 20/11/2022
Last modified: 25/05/2024
"""

import numpy as np
import pandas as pd

cols = [
    'm_r', 'kappa', 'a', 'nelec', 'nup', 'ndn', 'norb', 'nbasis', 'ncent',
    'nctype', 'gndot_v0', 'gndot_rho', 'gndot_s', 'gndot_k', 'gauss_sigma', 
    'scalek', 'nopt_iter', 'add_diag', 'p_var', 
    'etrial', 'eshift', 'gauss_sigma_best',
    'nelec_inside', 'uni_den_mean', 'uni_den_std',
    'ss_corr_den', 'ss_corr_pairden', 'edge_pol_den', 'edge_pol_pairden',
    'tot_E_1st', 'tot_E', 'tot_E_err', 'pop_err', 
    'stdev_1st', 'stdev', 'stdev_err', 
    'pot_E', 'pot_E_err', 'int_E', 'int_E_err', 
    'jf_kin_E', 'jf_kin_E_err', 'pb_kin_E', 'pb_kin_E_err', 
    'acceptance_1st', 'acceptance', 'delta', 'Tcorr_1st', 'Tcorr',
    'overlap_dmc', 'nwalk_eff', 'nwalk_eff_f', 'nwalk_eff_fs',
    'wts', 'wts_f', 'wts_fs', 'n_walkers_last', 'n_walkers_global',
    'n_walkers',
]

df_new = pd.read_csv('data_SI_new.csv', header=0, index_col=0)
df_old = pd.read_csv('data_SI_old.csv', header=0, index_col=0)

sort_cols = ['run_dir', 'flk_type', 'info']

df_new = df_new.sort_values(sort_cols,  ignore_index=True)
df_old = df_old.sort_values(sort_cols,  ignore_index=True)

new = df_new[cols].to_numpy()
old = df_old[cols].to_numpy()

new[np.isnan(new)] = 0
old[np.isnan(old)] = 0

diff = np.sqrt((new - old)**2).sum()

print(diff)


