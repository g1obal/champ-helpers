"""
CHAMP Data Reader and Extrapolated Estimator Calculator

Author: Gokhan Oztarhan
Created date: 24/06/2022
Last modified: 06/03/2023
"""

import sys
import os
import time
from copy import deepcopy
import pickle

import numpy as np
import pandas as pd

from champio.outputparser import OutputParser

VMC_ROOT_DIR = '../vmc/den'
DMC_ROOT_DIR = '../dmc/mpi2'
OUTPUT_FILE_NAME = 'output_file'

MPI2_MEAN = False
DATA_MEAN_FILE_NAME = 'mean_file.pkl'

EXTRAPOLATED_FILE_NAME = 'runs_ext.pkl'

# Manually set m_r, kappa and a (overwrites run title info)
M_R = None
KAPPA = None
A = None

# SI conversion units
EUNIT = 'meV' # meV, eV, J
LUNIT = 'nm' # A, nm, m

# PATHS
PATH_VMC = {}
PATH_DMC = {}

# VMC Directories 
for root, dirs, files in os.walk(VMC_ROOT_DIR):
    if OUTPUT_FILE_NAME in files:
        run_dir = os.path.split(root)[-1]
        PATH_VMC[run_dir] = root

# DMC Directories 
if MPI2_MEAN:
    search_fname = DATA_MEAN_FILE_NAME
else:
    search_fname = OUTPUT_FILE_NAME
for root, dirs, files in os.walk(DMC_ROOT_DIR):
    if search_fname in files:
        run_dir = os.path.split(root)[-1]
        PATH_DMC[run_dir] = root


def get_data_ext():
    tic = time.time()
    
    # Open pickle file
    pkl_file = open(EXTRAPOLATED_FILE_NAME, 'wb')
    
    # Form an empty dataframe using parser features
    parser = OutputParser()
    df = pd.DataFrame(columns=parser.features)
    df_SI = pd.DataFrame(columns=parser.features)
    
    for run_dir in PATH_DMC:
        try:
            parser_vmc = get_parser_vmc(run_dir)
            parser_dmc = get_parser_dmc(run_dir)
            
            # Path info
            print('vmc path: %s\n' %parser_vmc.path \
                + 'dmc path: %s' %parser_dmc.path)
            
            # Calculate den2d related extrapolated estimators
            parser_dmc = den2d_ext(parser_vmc, parser_dmc)
            
            # Calculate pairden related extrapolated estimators
            parser_dmc = pairden_ext(parser_vmc, parser_dmc)
            
            # Save all variables as dictionary to a pickle file
            pickle.dump(parser_dmc.__dict__, pkl_file, pickle.HIGHEST_PROTOCOL)
            pkl_file.flush()
            
            # Append to dataframe
            df = pd.concat(
                [df, parser_dmc.dataframe()], ignore_index=True, axis=0
            )
            
            # Convert to SI units 
            parser_dmc.to_SI(eunit=EUNIT, lunit=LUNIT)
            
            # Append to dataframe
            df_SI = pd.concat(
                [df_SI, parser_dmc.dataframe()], ignore_index=True, axis=0
            )
            
            print('Done.\n')
        except OSError as err:
            print('%s in files.\n' %type(err).__name__)
            
    # Save dataframes to csv
    df.to_csv('data_au.csv', header=True, float_format='% .6f')
    df_SI.to_csv('data_SI.csv', header=True, float_format='% .6f')
    
    # Close pickle file
    pkl_file.close()
    
    toc = time.time() 
    print("Execution time, get_data_ext = %.3f s" %(toc-tic))


def get_parser_vmc(run_dir):
    parser_vmc = OutputParser(
        OUTPUT_FILE_NAME,
        path=PATH_VMC[run_dir],
        m_r=M_R,
        kappa=KAPPA,
        a=A,
        calculate_ss_corr=False,
        calculate_edge_pol=False,
    )
    parser_vmc.parse()
    
    return parser_vmc


def get_parser_dmc(run_dir):
    if MPI2_MEAN:
        pkl_mean_fname = os.path.join(PATH_DMC[run_dir], DATA_MEAN_FILE_NAME)
        with open(pkl_mean_fname, 'rb') as pkl_mean_file:
            data = pickle.load(pkl_mean_file)
        
        # Initialize DMC parser using empty parser
        # Update it with the loaded dictionary
        parser_dmc = OutputParser()
        parser_dmc.__dict__.update(data)
        
        if parser_dmc.path not in PATH_DMC[run_dir]:
            print('Warning: parser_dmc.path not in PATH_DMC[run_dir]')
            print('parser_dmc.path (read from pkl file) = %s' %parser_dmc.path)
        
        parser_dmc.path = PATH_DMC[run_dir]
        
        if M_R is not None:
            parser_dmc.m_r = M_R
        if KAPPA is not None:
            parser_dmc.kappa = KAPPA
        if A is not None:
            parser_dmc.a = A        
    else:
        parser_dmc = OutputParser(
            OUTPUT_FILE_NAME,
            path=PATH_DMC[run_dir],
            m_r=M_R,
            kappa=KAPPA,
            a=A,
            calculate_ss_corr=False,
            calculate_edge_pol=False,
        )  
        parser_dmc.parse()

    return parser_dmc


def den2d_ext(parser_vmc, parser_dmc):
    try:
        # Extrapolated estimator for den2d        
        parser_dmc.den2d_t[:,:,2] = 2 * parser_dmc.den2d_t[:,:,2] \
                                    - parser_vmc.den2d_t[:,:,2]
        parser_dmc.den2d_s[:,:,2] = 2 * parser_dmc.den2d_s[:,:,2] \
                                    - parser_vmc.den2d_s[:,:,2]

        # Integral of total electron density
        dx = parser_dmc.den2d_t[1,0,0] - parser_dmc.den2d_t[0,0,0]
        dy = parser_dmc.den2d_t[0,1,1] - parser_dmc.den2d_t[0,0,1]
        parser_dmc.den2d_nelec_calc = np.trapz(
            np.trapz(parser_dmc.den2d_t[:,:,2], dx=dx, axis=0), dx=dy, axis=0
        )
        
        # Normalization of spin density:
        # Integral of absolute value of spin density should be equal to 
        # integral of total electron density if spins are perfectly polarized.
        # While calculating extrapolated estimator for spin density,
        # if VMC spin density is close to metallic phase and DMC spin density
        # is close to antiferromagnetic phase, resulting extrapolated spin 
        # density would be an antiferromagnetic phase however the normalization 
        # of it would be corrupted (since 2*DMC - VMC). It should be corrected 
        # using integral of total electron density.
        abs_den2d_s_sum = np.trapz(
            np.trapz(np.abs(parser_dmc.den2d_s[:,:,2]), dx=dx, axis=0), 
            dx=dy, axis=0
        )
        
        # Normalize spin density
        if abs_den2d_s_sum > parser_dmc.den2d_nelec_calc:
            parser_dmc.den2d_s[:,:,2] *= \
                parser_dmc.den2d_nelec_calc / abs_den2d_s_sum
        
        # Calculate real space spin-spin correlation for density 
        parser_dmc.ss_corr_den = \
            parser_dmc._ss_corr(parser_dmc.den2d_t, parser_dmc.den2d_s)
            
        # Calculate edge polarization for density
        parser_dmc.edge_pol_den = parser_dmc._edge_pol(parser_dmc.den2d_s)
    
    except (TypeError, AttributeError) as err:
        print('%s in den2d_ext.' %type(err).__name__)
        parser_dmc.den2d_t = np.nan
        parser_dmc.den2d_s = np.nan
        parser_dmc.den2d_nelec_calc = np.nan
        parser_dmc.ss_corr_den = np.nan
        parser_dmc.edge_pol_den = np.nan
    
    return parser_dmc


def pairden_ext(parser_vmc, parser_dmc):
    try:
        # Extrapolated estimator for pairden
        parser_dmc.pairden_dd[:,:,2] = 2 * parser_dmc.pairden_dd[:,:,2] \
                                       - parser_vmc.pairden_dd[:,:,2]
        parser_dmc.pairden_dt[:,:,2] = 2 * parser_dmc.pairden_dt[:,:,2] \
                                       - parser_vmc.pairden_dt[:,:,2]
        parser_dmc.pairden_du[:,:,2] = 2 * parser_dmc.pairden_du[:,:,2] \
                                       - parser_vmc.pairden_du[:,:,2]
        parser_dmc.pairden_ud[:,:,2] = 2 * parser_dmc.pairden_ud[:,:,2] \
                                       - parser_vmc.pairden_ud[:,:,2]
        parser_dmc.pairden_ut[:,:,2] = 2 * parser_dmc.pairden_ut[:,:,2] \
                                       - parser_vmc.pairden_ut[:,:,2]
        parser_dmc.pairden_uu[:,:,2] = 2 * parser_dmc.pairden_uu[:,:,2] \
                                       - parser_vmc.pairden_uu[:,:,2]
        parser_dmc.pairden_t[:,:,2] = 2 * parser_dmc.pairden_t[:,:,2] \
                                      - parser_vmc.pairden_t[:,:,2]                             
        parser_dmc.pairden_s[:,:,2] = 2 * parser_dmc.pairden_s[:,:,2] \
                                      - parser_vmc.pairden_s[:,:,2]
             
        # Integration variables (for normalization)
        dx = parser_dmc.pairden_dd[1,0,0] - parser_dmc.pairden_dd[0,0,0]
        dy = parser_dmc.pairden_dd[0,1,1] - parser_dmc.pairden_dd[0,0,1]    
               
        # Normalization of spin density:
        # Integral of absolute value of spin density should be equal to 
        # integral of total electron density if spins are perfectly polarized.
        # While calculating extrapolated estimator for spin density,
        # if VMC spin density is close to metallic phase and DMC spin density
        # is close to antiferromagnetic phase, resulting extrapolated spin 
        # density would be an antiferromagnetic phase however the normalization 
        # of it would be corrupted (since 2*DMC - VMC). It should be corrected 
        # using integral of total electron density.
        pairden_t_sum = np.trapz(
            np.trapz(parser_dmc.pairden_t[:,:,2], dx=dx, axis=0), dx=dy, axis=0
        )
        abs_pairden_s_sum = np.trapz(
            np.trapz(np.abs(parser_dmc.pairden_s[:,:,2]), dx=dx, axis=0), 
            dx=dy, axis=0
        )
        
        # Normalize spin density
        if abs_pairden_s_sum > pairden_t_sum:
            parser_dmc.pairden_s[:,:,2] *= pairden_t_sum / abs_pairden_s_sum
                    
        # Calculate real space spin-spin correlation for pair density 
        parser_dmc.ss_corr_pairden = \
            parser_dmc._ss_corr(parser_dmc.pairden_t, parser_dmc.pairden_s)
            
        # Calculate edge polarization for pair density
        parser_dmc.edge_pol_pairden = parser_dmc._edge_pol(parser_dmc.pairden_s)

    except (TypeError, AttributeError) as err:
        print('%s in pairden_ext.' %type(err).__name__)
        parser_dmc.pairden_dd = np.nan
        parser_dmc.pairden_dt = np.nan
        parser_dmc.pairden_du = np.nan
        parser_dmc.pairden_ud = np.nan
        parser_dmc.pairden_ut = np.nan
        parser_dmc.pairden_uu = np.nan
        parser_dmc.pairden_t = np.nan
        parser_dmc.pairden_s = np.nan
        parser_dmc.ss_corr_pairden = np.nan
        parser_dmc.edge_pol_pairden = np.nan
        
    return parser_dmc


if __name__ == '__main__':
    get_data_ext()

