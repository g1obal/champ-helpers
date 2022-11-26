"""
CHAMP Data Reader

Author: Gokhan Oztarhan
Created date: 10/05/2021
Last modified: 23/11/2022
"""

import os
import time

import pandas as pd

from champio.outputparser import OutputParser


ROOT_DIR = '.'
OUTPUT_FILE_NAME = 'output_file'

# Manually set m_r, kappa and a (overwrites run title info)
M_R = None
KAPPA = None
A = None

# SI conversion units
EUNIT = 'meV' # meV, eV, J
LUNIT = 'nm' # A, nm, m


def get_data():
    tic = time.time()
    
    # Form an empty dataframe using parser features
    parser = OutputParser()
    df = pd.DataFrame(columns=parser.features)
    df_SI = pd.DataFrame(columns=parser.features)

    for root, dirs, files in os.walk(ROOT_DIR):
        dirs.sort()
        if OUTPUT_FILE_NAME in files:    
            # Parse output file        
            parser = OutputParser(
                OUTPUT_FILE_NAME, 
                path=root,
                m_r=M_R,
                kappa=KAPPA,
                a=A,
            )  
            parser.parse()  
            
            # Append to dataframe
            df = pd.concat(
                [df, parser.dataframe()], ignore_index=True, axis=0
            )
            # Convert to SI units 
            parser.to_SI(eunit=EUNIT, lunit=LUNIT)
            
            # Append to dataframe
            df_SI = pd.concat(
                [df_SI, parser.dataframe()], ignore_index=True, axis=0
            )
            
            # Print info
            print('Done OutputParser: %s' %os.path.split(root)[-1])
    
    # Save dataframes to csv
    df.to_csv('data_au.csv', header=True, float_format='% .6f')
    df_SI.to_csv('data_SI.csv', header=True, float_format='% .6f')

    toc = time.time() 
    print("Execution time, get_data = %.3f s" %(toc-tic))

if __name__ == '__main__':
    get_data()

 
