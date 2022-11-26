"""
clear_dmc_dirs

Removes VMC density files from DMC runs directories

Author: Gokhan Oztarhan
Created date: 27/06/2022
Last modified: 13/08/2022
"""

import os


ROOT_DIR = '.'
INPUT_FILE_NAME = 'input_file'
OUTPUT_FILE_NAME = 'output_file'

RM_LIST = [
    'den2d_t_vmc',
    'den2d_d_vmc',
    'den2d_u_vmc',
    
    'pairden_ut_vmc',
    'pairden_ud_vmc',
    'pairden_uu_vmc',
    
    'pairden_dt_vmc',
    'pairden_dd_vmc',
    'pairden_du_vmc',
    
    'pot_t_vmc',
    'pot_d_vmc',
    'pot_u_vmc',
]


def clear_dmc_dirs():
    for root, dirs, files in os.walk(ROOT_DIR):
        dirs.sort()
        if OUTPUT_FILE_NAME in files and 'dmc' in root:
            if 'pairden_ut_dmc' in files:
                print(root)
            
                for f in RM_LIST:
                    try:
                        os.remove(os.path.join(root, f))
                    except Exception as err:
                        print(type(err).__name__, f)
            
            print('')
                     

if __name__ == '__main__':
    clear_dmc_dirs()
    
    

    
