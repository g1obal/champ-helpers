"""
edit_input

Edits currently present inputs using InputEditor.

Author: Gokhan Oztarhan
Created date: 30/01/2022
Last modified: 23/09/2022
"""

import os

import numpy as np

from champio.outputparser import OutputParser
from champio.inputeditor import InputEditor


IRN = 'auto' # random seed
ETRIAL = 'auto'
NSTEP = 100
NBLK = 10
NBLKEQ = 2
N_WALKERS = 5
DELTA = 'auto'
TAU = 0.1
IFIXE = -3 # 0: no output, -1: 2d density, -3: 2d and pair density
XMAX = None
NOPT_ITER = 0
OPT_MODE = 0 # 0: both, 1: only width, 2: only jastrow

UPDATE_GAUSS_WIDTH = True
UPDATE_JASTROW = True

ROOT_DIR = '.'
INPUT_FILE_NAME = 'input_file'
OUTPUT_FILE_NAME = 'output_file'
NEW_INPUT_FILE_NAME = 'input_file'

CLEAN_RUN_DIR = True
EXCLUDE_FILES = [
    INPUT_FILE_NAME,
    NEW_INPUT_FILE_NAME,
    'orb_dot_coef',
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
            editor.edit_delta(parser.acceptance, DELTA)
            editor.edit_tau(TAU)
            editor.edit_ifixe(IFIXE)
            editor.edit_xmax(XMAX)
            editor.edit_nopt_iter(NOPT_ITER)
            editor.edit_opt_mode(OPT_MODE)
            
            if UPDATE_GAUSS_WIDTH:
                editor.update_gauss_width(parser.opt_gauss_width())
            if UPDATE_JASTROW:
                editor.update_jastrow(parser.opt_jastrow())
                
            editor.commit(NEW_INPUT_FILE_NAME)
            
            print('Done InputEditor.\n')
            
            if CLEAN_RUN_DIR:
                remove_files(root, files, EXCLUDE_FILES)
                
                
def remove_files(path, files, exclude):
    for f in files:
        if f not in exclude:
            os.remove(os.path.join(path, f))
                     

if __name__ == '__main__':
    edit_input()
    
    

    
