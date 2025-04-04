"""
edit_input

Edits currently present inputs using InputEditor.

Author: Gokhan Oztarhan
Created date: 30/01/2022
Last modified: 27/05/2024
"""

import os
import sys

import numpy as np

from champio.outputparser import OutputParser
from champio.inputeditor import InputEditor


IRN = 'auto' # random seed
ETRIAL = 'auto'
NSTEP = 1000
NBLK = 100
NBLKEQ = 5
N_WALKERS = 5
NCONF_NEW = 100 # save n_conf_new final walkers per cpu
ISITE = 1 # 0: read walkers from file, 1: generate walkers in sites.f90
DELTA = 'auto'
TAU = 0.1
IFIXE = -3 # 0: no output, -1: 2d density, -3: 2d and pair density
XMAX = None
NOPT_ITER = 0
ADD_DIAG = None
P_VAR = None
OPT_MODE = None # 0: all, 1: width+pos, 2: width+jastrow, 3: pos+jastrow
                # 4: only pos, 5: only width, 6: only jastrow

UPDATE_POS = True
UPDATE_GAUSS_WIDTH = True
UPDATE_JASTROW = True

GENERATE_MC_CONFIGS = False
MC_CONFIGS_PICK_RANDOM = False
N_MC_CONFIGS = 120

ROOT_DIR = '.'
INPUT_FILE_NAME = 'input_file'
OUTPUT_FILE_NAME = 'output_file'
NEW_INPUT_FILE_NAME = 'input_file'

CLEAN_RUN_DIR = True
EXCLUDE_FILES = [
    INPUT_FILE_NAME,
    NEW_INPUT_FILE_NAME,
    'orb_dot_coef',
    'mc_configs',
] 


def edit_input():
    for root, dirs, files in os.walk(ROOT_DIR):
        dirs.sort()
        if INPUT_FILE_NAME in files:
            print('Directory: %s' %root)
        
            parser = OutputParser(
                OUTPUT_FILE_NAME, path=root, 
                parse_density=False, 
                calculate_ss_corr=False, 
                calculate_edge_pol=False,
                calculate_U_onsite=False,
                calculate_nelec_inside_uni=False,
            )
            parser.parse()
            
            editor = InputEditor(INPUT_FILE_NAME, path=root)
            editor.edit_irn(IRN)
            editor.edit_etrial(
                parser.run_mode, parser.tot_E, parser.tot_E_err, 
                parser.pop_err, ETRIAL
            )
            editor.edit_nstep(NSTEP)
            editor.edit_nblk(NBLK)
            editor.edit_nblkeq(NBLKEQ)
            editor.edit_nconf(N_WALKERS)
            editor.edit_nconf_new(NCONF_NEW)
            editor.edit_isite(ISITE)
            editor.edit_delta(parser.acceptance, DELTA)
            editor.edit_tau(TAU)
            editor.edit_ifixe(IFIXE)
            editor.edit_xmax(XMAX)
            editor.edit_nopt_iter(NOPT_ITER)
            editor.edit_add_diag(ADD_DIAG)
            editor.edit_p_var(P_VAR)
            editor.edit_opt_mode(parser.nbasis, parser.constraint, OPT_MODE)
            
            pos, width, jast = False, False, False
            if parser.opt_type == 'all':
                pos, width, jast = True, True, True
            elif parser.opt_type == 'gwidth+pos':
                pos, width, jast = True, True, False
            elif parser.opt_type == 'gwidth+jastrow':
                pos, width, jast = False, True, True
            elif parser.opt_type == 'pos+jastrow':
                pos, width, jast = True, False, True
            elif parser.opt_type == 'pos':
                pos, width, jast = True, False, False
            elif parser.opt_type == 'gwidth':
                pos, width, jast = False, True, False
            elif parser.opt_type == 'jastrow':
                pos, width, jast = False, False, True

            if UPDATE_POS and pos:
                editor.update_pos(parser.opt_pos())
            if UPDATE_GAUSS_WIDTH and width:
                editor.update_gauss_width(parser.opt_gauss_width())
            if UPDATE_JASTROW and jast:
                editor.update_jastrow(parser.nctype, parser.opt_jastrow())
            
            if GENERATE_MC_CONFIGS:
                editor.generate_mc_configs(
                    root, files, N_MC_CONFIGS, MC_CONFIGS_PICK_RANDOM
                )
                
            editor.commit(NEW_INPUT_FILE_NAME)
            
            print('Done InputEditor.\n')
            
            if CLEAN_RUN_DIR:
                remove_files(root, files, EXCLUDE_FILES)
                
                
def remove_files(root, files, exclude):
    for f in files:
        if f not in exclude:
            os.remove(os.path.join(root, f))
                     

if __name__ == '__main__':
    raise_cont = False
    if ISITE is not None:
        if ISITE and GENERATE_MC_CONFIGS or \
        not ISITE and not GENERATE_MC_CONFIGS:
            raise_cont = True
    else:
        if GENERATE_MC_CONFIGS:
            raise_cont = True
    
    if raise_cont:
        print('ISITE = %s\nGENERATE_MC_CONFIGS = %s' \
            %(ISITE, GENERATE_MC_CONFIGS))
        cont = input('Continue? [y/n]\n')
        if cont == 'y':
            print('')
        elif cont == 'n':
            sys.exit('Stopped.')
        else:
            sys.exit('Stopped: Please enter y or n.')
    
    edit_input()
    

