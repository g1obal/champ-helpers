"""
get_data_ens

Calculate the ensemble averages of extrapolated data.

Author: Gokhan Oztarhan
Created date: 02/08/2022
Last modified: 15/03/2023
"""

import sys
import os
import time
from copy import deepcopy
import pickle

import numpy as np
import pandas as pd

from champio.auconverter import AUConverter
from champio.outputparser import OutputParser


ROOT_DIR = '.'
OUTPUT_FILE_NAME = 'output_file'

DATA_EXT_DIR = '../ext'
DATA_EXT_FILE_NAME = 'ext_file.pkl'

DATA_ENS_DIR = 'runs_ens'
DATA_ENS_FILE_NAME = 'ens_file.pkl'

# Temperatures (K)
TEMP = [4]

# Manually set m_r, kappa and a (overwrites run title info)
M_R = None
KAPPA = None
A = None

# SI conversion units
EUNIT = 'meV' # meV, eV, J
LUNIT = 'nm' # A, nm, m

# Features
FEATURES = [
    ['gauss_sigma', None],
    ['etrial', None],
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
    ['delta', None],
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

# Load pickle files
groups = {}
for root, dirs, files in sorted(os.walk(DATA_EXT_DIR)):
    if DATA_EXT_FILE_NAME in files:
        with open(os.path.join(root, DATA_EXT_FILE_NAME), 'rb') as pkl_file_ext:
            data = pickle.load(pkl_file_ext)
        
        parser = OutputParser()
        parser.__dict__.update(data)
        
        if M_R is not None:
            parser.m_r = M_R
        if KAPPA is not None:
            parser.kappa = KAPPA
        if A is not None:
            parser.a = A
        
        # Get group names of matching runs
        group_name = parser.run_dir.split('_')[0]
        if group_name in groups:
            groups[group_name].append(deepcopy(parser))
        else:
            groups[group_name] = [deepcopy(parser)]

# Create ensemble directory
if not os.path.exists(DATA_ENS_DIR):
    os.mkdir(DATA_ENS_DIR)


def get_data_ens():
    tic = time.time()
    
    # Form an empty dataframe using parser features
    parser = OutputParser()
    
    # Dataframe for writing mean data
    df = pd.DataFrame(columns=parser.features + ['n_data_mean'])
    df_SI = pd.DataFrame(columns=parser.features + ['n_data_mean'])
    
    # Dataframe for writing all data
    df_all = pd.DataFrame(columns=parser.features)
    df_all_SI = pd.DataFrame(columns=parser.features)
    
    for temp in TEMP:
        print('temp = %s K' %temp)
        
        for group_name in groups:
            print('group_name = %s' %group_name)
            
            # Deep copy needed since at the end of the loop 
            # _parser.to_SI() edits all parser values.
            parser_list = deepcopy(groups[group_name])
            
            if parser_list:
                # Calculate mean values
                parser_mean = deepcopy(parser_list[-1])
                
                # Boltzmann constant
                kb_SI = 1.380649e-23 # in J K^-1
                kb = AUConverter(parser_mean.m_r, parser_mean.kappa).\
                    energy_to_au(kb_SI, 'J')
                beta = 1 / (kb * temp)
                
                shift = calculate_exponent_shift(beta, parser_list)
                
                for _features in FEATURES:
                    parser_mean = get_mean(
                        *_features, beta, shift, parser_list, parser_mean
                    )
                parser_mean._ratio_int_kin_dmc()
                parser_mean = set_ss_corrs(parser_mean)
                parser_mean = get_time_mean(
                    beta, shift, parser_list, parser_mean
                )
                
                # Set variables of parser_mean
                parser_mean.fname = DATA_ENS_FILE_NAME
                parser_mean.parent_path = DATA_ENS_DIR
                parser_mean.run_dir = \
                    '_'.join(
                        parser_mean.run_dir.split('_')[:-1] + [str(temp) + 'K']
                    )
                parser_mean.path = \
                    os.path.join(parser_mean.parent_path, parser_mean.run_dir)
                parser_mean.info = 'T = %s K' %temp
                
                # Create run_dir for ensemble file
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
        
    # Add ensemble data to df_all in order to plot them together
    df_all = pd.concat([df_all, df], ignore_index=True, axis=0)
    df_all_SI = pd.concat([df_all_SI, df_SI], ignore_index=True, axis=0)
        
    # Save dataframes to csv
    df.to_csv('data_au.csv', header=True, float_format='% .6f')
    df_SI.to_csv('data_SI.csv', header=True, float_format='% .6f')
    df_all.to_csv('data_all_au.csv', header=True, float_format='% .6f')
    df_all_SI.to_csv('data_all_SI.csv', header=True, float_format='% .6f')
    
    toc = time.time() 
    print("Execution time, get_data_ens = %.3f s" %(toc-tic))
        
        
def get_mean(feature, feature_err, beta, shift, parser_list, parser_mean):
    # Weighted variance calculation sources;
    # https://en.wikipedia.org/wiki/Variance
    # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean

    weights = []
    norm = 0
    for i, parser in enumerate(parser_list):
        weight = np.exp(- beta * parser.__dict__['tot_E'] - shift)
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
    
    
def get_time_mean(beta, shift, parser_list, parser_mean):
    try:
        seconds_mean = 0
        norm = 0
        for parser in parser_list:
            weight = np.exp(- beta * parser.__dict__['tot_E'] - shift)
                
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
    except Exception:
        parser_mean.tot_cpu_time = np.nan
            
    return parser_mean


def calculate_exponent_shift(beta, parser_list):
    # Calculate the shift in exponents if exponent is greater than 150
    # in order to prevent overflow
    exponents = []
    for i, parser in enumerate(parser_list):
        exponents.append(- beta * parser.__dict__['tot_E'])
        
    bools = []
    for i, parser in enumerate(parser_list):
        bools.append(abs(exponents[i]) >= 150)
        
    if any(bools):
        ind = np.argmax(np.abs(exponents))
        shift = exponents[ind] - np.sign(exponents[ind]) * 150
    else:
        shift = 0
        
    return shift
    

if __name__ == '__main__':
    get_data_ens()


