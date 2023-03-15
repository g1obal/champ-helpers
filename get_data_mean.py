"""
get_data_mean

Calculate the mean of mpi2 runs with different random seeds.

Author: Gokhan Oztarhan
Created date: 04/07/2022
Last modified: 15/03/2023
"""

import sys
import os
import time
from copy import deepcopy
import pickle

import numpy as np
import pandas as pd

from champio.outputparser import OutputParser


ROOT_DIR = '.'
OUTPUT_FILE_NAME = 'output_file'

DATA_MEAN_DIR = 'runs_mean'
DATA_MEAN_FILE_NAME = 'mean_file.pkl'

# Manually set m_r, kappa and a (overwrites run title info)
M_R = None
KAPPA = None
A = None

# SI conversion units
EUNIT = 'meV' # meV, eV, J
LUNIT = 'nm' # A, nm, m

# Weight features in order to calculate the mean, 
# these features are multiplied with each other.
# If empty list, all weights are equal to 1.
WEIGHT_FEATURES = ['n_walkers', 'n_cpu']

# Features
FEATURES = [
    ['tot_E_1st', None],
    ['tot_E', 'tot_E_err'],
    ['tot_E', 'pop_err'],
    ['stdev_1st', None],
    ['stdev', 'stdev_err'],
    ['pot_E', 'pot_E_err'],
    ['int_E', 'int_E_err'],
    ['jf_kin_E', 'jf_kin_E_err'],
    ['pb_kin_E', 'pb_kin_E_err'],
    ['acceptance_1st', None],
    ['acceptance', None],
    ['Tcorr_1st', None],
    ['Tcorr', None],
    ['overlap_dmc', None],
    ['nwalk_eff', None],
    ['nwalk_eff_f', None],
    ['nwalk_eff_fs', None],
    ['wts', None],
    ['wts_f', None],
    ['wts_fs', None],
    ['n_walkers_last', None],
    ['n_walkers_global', None],
    ['n_walkers', None],
    ['n_cpu', None],
    ['n_block_eq', None],
    ['n_block', None],
    ['n_steps', None],
    ['den2d_t', None],
    ['den2d_s', None],
    ['den2d_nelec_calc', None],
    ['pairden_dd', None],
    ['pairden_dt', None],
    ['pairden_du', None],
    ['pairden_ud', None],
    ['pairden_ut', None],
    ['pairden_uu', None],
    ['pairden_t', None],
    ['pairden_s', None],
]

# Matching directories, using sorted() to guarantee os.walk() order.
# Notice that using sorted on os.walk makes os.walk to return all values at once
# instead of returning the iteration elements one by one since the path list
# needs to be sorted before iterating on it.
PATH = {}
for root, dirs, files in sorted(os.walk(ROOT_DIR)):
    if OUTPUT_FILE_NAME in files:
        run_dir = os.path.split(root)[-1]
        if run_dir in PATH:
            PATH[run_dir].append(root)
        else:
            PATH[run_dir] = [root]

# Create mean directory
if not os.path.exists(DATA_MEAN_DIR):
    os.mkdir(DATA_MEAN_DIR)


def get_data_mean():
    tic = time.time()

    # Form an empty dataframe using parser features
    parser = OutputParser()
    
    # Dataframe for writing mean data
    df = pd.DataFrame(columns=parser.features + ['n_data_mean'])
    df_SI = pd.DataFrame(columns=parser.features + ['n_data_mean'])
    
    # Dataframe for writing all data
    df_all = pd.DataFrame(columns=parser.features)
    df_all_SI = pd.DataFrame(columns=parser.features)
    
    for run_dir in PATH:
        print('run_dir = %s' %run_dir)
        
        parser_list = []
        for path in PATH[run_dir]:
            # Parse DMC output file and add to parser list
            parser_list.append(
                OutputParser(
                    OUTPUT_FILE_NAME, 
                    path=path,
                    m_r=M_R,
                    kappa=KAPPA,
                    a=A,
                )
            )  
            parser_list[-1].parse()
            print('OutputParser done: %s' %path)
        
        if parser_list:
            # Calculate mean values
            parser_mean = deepcopy(parser_list[-1])
            for _features in FEATURES:
                parser_mean = get_mean(*_features, parser_list, parser_mean)
            parser_mean._ratio_int_kin_dmc()
            parser_mean = set_ss_corrs(parser_mean)
            parser_mean = get_time_mean(parser_list, parser_mean)
            
            # Set variables of parser_mean
            parser_mean.fname = DATA_MEAN_FILE_NAME
            parser_mean.parent_path = DATA_MEAN_DIR
            parser_mean.path = \
                os.path.join(parser_mean.parent_path, parser_mean.run_dir)
            
            # Create run_dir for mean file
            if not os.path.exists(parser_mean.path):
                os.mkdir(parser_mean.path)
            
            # Save all variables as dictionary to a pickle file
            pkl_file = open(
                os.path.join(parser_mean.path, parser_mean.fname), 'wb'
            )
            pickle.dump(
                parser_mean.__dict__, pkl_file, pickle.HIGHEST_PROTOCOL
            )
            pkl_file.close()
            
            # Append to dataframe
            df_temp = parser_mean.dataframe()
            df_temp['n_data_mean'] = len(parser_list)
            df = pd.concat([df, df_temp], ignore_index=True, axis=0)
            
            # Convert to SI units 
            parser_mean.to_SI(eunit=EUNIT, lunit=LUNIT)
            
            # Append to dataframe
            df_temp = parser_mean.dataframe()
            df_temp['n_data_mean'] = len(parser_list)
            df_SI = pd.concat([df_SI, df_temp], ignore_index=True, axis=0)
            
            # Append all runs into df_all and df_all_SI
            for _parser in parser_list:
                df_all = pd.concat(
                    [df_all, _parser.dataframe()], 
                    ignore_index=True, axis=0
                )
                _parser.to_SI(eunit=EUNIT, lunit=LUNIT)
                df_all_SI = pd.concat(
                    [df_all_SI, _parser.dataframe()], 
                    ignore_index=True, axis=0
                )
            
        print('')
        
    # Save dataframes to csv
    df.to_csv('data_au.csv', header=True, float_format='% .6f')
    df_SI.to_csv('data_SI.csv', header=True, float_format='% .6f')
    df_all.to_csv('data_all_au.csv', header=True, float_format='% .6f')
    df_all_SI.to_csv('data_all_SI.csv', header=True, float_format='% .6f')
    
    toc = time.time() 
    print("Execution time, get_data_mean = %.3f s" %(toc-tic))
        
        
def get_mean(feature, feature_err, parser_list, parser_mean):
    # Weighted variance calculation sources;
    # https://en.wikipedia.org/wiki/Variance
    # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
    
    weights = []
    norm = 0
    for i, parser in enumerate(parser_list):
        weight = 1.0
        for weight_feature in WEIGHT_FEATURES:
            weight *= parser.__dict__[weight_feature]
        
        norm += weight    
        weights.append(weight)
    
    mean = 0
    for i, parser in enumerate(parser_list):
        mean +=  weights[i] * parser.__dict__[feature]
    
    mean /= norm
    parser_mean.__dict__[feature] = mean

    if feature_err is not None:
        Not_NaN = []
        for parser in parser_list:
            Not_NaN.append(~np.isnan(parser.__dict__[feature_err]))
                 
        var = 0
        if all(Not_NaN):
            for i, parser in enumerate(parser_list):
                var += weights[i]**2 * parser.__dict__[feature_err]**2
            var /= norm**2
        else:
            for i, parser in enumerate(parser_list):
                var += weights[i] * (mean - parser.__dict__[feature])**2
            var /= norm
        
        parser_mean.__dict__[feature_err] = np.sqrt(var)
        
    return parser_mean
    
    
def set_ss_corrs(parser_mean):
    # Integral of electron density
    dx = parser_mean.den2d_t[1,0,0] - parser_mean.den2d_t[0,0,0]
    dy = parser_mean.den2d_t[0,1,1] - parser_mean.den2d_t[0,0,1]
    parser_mean.den2d_nelec_calc = np.trapz(
        np.trapz(parser_mean.den2d_t[:,:,2], dx=dx, axis=0), dx=dy, axis=0
    )
    
    # Calculate real space spin-spin correlation for density
    parser_mean.ss_corr_den = \
        parser_mean._ss_corr(parser_mean.den2d_t, parser_mean.den2d_s)

    # Calculate real space spin-spin correlation for pair density
    parser_mean.ss_corr_pairden = \
        parser_mean._ss_corr(parser_mean.pairden_t, parser_mean.pairden_s)
    
    # Calculate edge polarization for density
    parser_mean.edge_pol_den = parser_mean._edge_pol(parser_mean.den2d_s)
    
    # Calculate edge polarization for pair density
    parser_mean.edge_pol_pairden = parser_mean._edge_pol(parser_mean.pairden_s)
        
    return parser_mean
    
    
def get_time_mean(parser_list, parser_mean):
    seconds_mean = 0
    norm = 0
    for parser in parser_list:
        weight = 1.0
        for weight_feature in WEIGHT_FEATURES:
            weight *= parser.__dict__[weight_feature]
            
        cpu_time = parser.tot_cpu_time.split()
        hours = int(cpu_time[0][:-1])
        mins = int(cpu_time[1][:-1])
        secs = int(cpu_time[2][:-1])
        
        seconds_mean += weight * (secs + mins * 60 + hours * 3600)
        norm += weight
        
    seconds_mean /= norm
    
    hours, remaining = divmod(seconds_mean, 3600)
    mins, secs = divmod(remaining, 60)
    
    parser_mean.tot_cpu_time = '%dh %dm %ds' %(hours, mins, secs)  
        
    return parser_mean
    

if __name__ == '__main__':
    get_data_mean()

