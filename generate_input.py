"""
generate_input

Generating inputs for CHAMP program.

Author: Gokhan Oztarhan
Created date: 09/06/2019
Last modified: 24/05/2024
"""

import os
import sys
import shutil
import time
import datetime
import logging

import numpy as np

from champio.logger import set_logger, unset_logger
from champio.inputgenerator import InputGenerator


logger = logging.getLogger(__name__)

GENERATE_INPUT_MULTIPLE = False

PLOT = 0
PLOT_FNAME = 'model.png'
PLOT_DPI = 300
PLOT_INTERACTIVE = 0

CONFIG = {
    # [file]
    'champ_input_file_name': 'input_file',

    # [title]
    'info': None, 

    # [au]
    'm_r': 0.067, # Electron mass ratio, m_eff / m_e
    'kappa': 12.4,  # Dielectric constant, varepsilon / varepsilon_0

    # [units]
    'eunit': 'meV',
    'lunit': 'nm',

    # [potential]
    'gndot_v0': -20,
    'gndot_rho': 20,
    'gndot_s': 1.40, # unitless
    'gndot_k': 0, # unit: [ENERGY] / [LENGTH]^2

    # [basis]
    'gauss_sigma': 12, # width guess for Gaussian basis

    # [lattice]
    'a': 50.0,
    'n_side': 4,
    'width': 4, # for nanoribbon
    'bc': 'xy', # for nanoribbon
    'lat_type': 'honeycomb',
    'flk_type': 1, # 0: hexagonal_zigzag, 
                   # 1: hexagonal_armchair, 
                   # 2: triangular_zigzag, 
                   # 3: triangular_armchair, 
                   # 4: nanoribbon

    # [electron]
    'total_charge': 0, # set total number of electrons, nelec
                          # None or 0 for charge neutral system
    'Sz': None, # total spin; to calculate the number of up and down electrons
                # None to arrange nup and ndn according to Lieb's theorem
    'spin_order': 'AFM', # electrons are located in a spin order
                         # AFM: antiferromagnetic, FM: ferromagnetic
    'spin_order_direction': 1, # the direction in which electrons are located
                               # 0: add electrons from outside to inside, 
                               #    additional electrons from inside to outside.   
                               # 1: add electrons from inside to outside, 
                               #    additional electrons from outside to inside.  

    # [orbital]
    'orb_dot_coef': 0, # 0: gauss, 1: orbitals read from orb_dot_coef file 
                       # >0: overwrites spin_order

    # [random seed champ]
    'irn': 'auto', # irn: a string of 16 digit integer
                   # 'auto': automatically set random seed using os.urandom

    # [run]
    'nstep': 100, # number of steps per block
    'nblk': 10,
    'nblkeq': 3,

    # [dmc]
    'etrial': 0.1,
    'n_walkers': 5,
    'tau': 0.1,
    
    # [jastrow]
    'scalek': 0.2, # 0.2 is default
    'nctype_of_edges': 1, # 1: same iwtype for all dots
                          # 2: dots at the edges have different iwtype
    'nctype2_include_corners': 1, # 0: do not include corner sites
                                  # 1: include corner sites
    
    # [optional champ]
    'ifixe': 0, #  0: do not write 2d density, 
                # -1: write out 2d density
                # -3: write out full pair density and 2d density
    'xfix_pos': None, # fixed electron position, example: [0.2, 0.5]
                      # None to find automatically
    'xfix_angle': 60, # fixed electron symmetry, default is 60
    
    # [opt]
    'opt_mode': 0, # 0: all, 1: width+pos, 2: width+jastrow, 3: pos+jastrow
                   # 4: only pos, 5: only width, 6: only jastrow
    'opt_constraint': 0,
    'nopt_iter': 20,
    'nblk_max': 25, # max. number of blocks which can be increased during opt.
    'add_diag': 1e-4, # CHAMP uses abs(add_diag)
                      # negative sign for fixed add_diag
                      # positive sign for optimization of add_diag
    'p_var': 0.1, # 0: energy, 1:variance
    'tol_energy': 1e-12, # energy tolerance to finish optimization
    'iopt': '00002', # last digit 2 is newton, 
                     # 01002 also a good choice, 31101 is linear (bad choice)
    'ipr_opt': -2, # -2: minimal output, 2: maximum output
}


def generate_input():   
    tic = time.time()
    
    # Initialize logger to print only to console
    set_logger(0, 1)
    
    # Initial info of the script
    now = datetime.datetime.now()
    string = 'generate_input\n' + now.strftime('%Y-%m-%d %H:%M:%S\n\n')
    logger.info(string)  
    
    # Initialize InputGenerator
    inputgenerator = InputGenerator(CONFIG)
    
    # Set all variables for input file
    inputgenerator.set()
    
    # Writing settings to input file
    inputgenerator.write()
    
    # Plotting the figures
    if PLOT:
        inputgenerator.plot(
            plot_fname=PLOT_FNAME, plot_dpi=PLOT_DPI, 
            plot_interactive=PLOT_INTERACTIVE
        )
    
    # End of program
    toc = time.time()
    string = '\nExecution time = %.3f s\n' %(toc-tic)
    logger.info(string)
    
    # Close all handlers of logger
    unset_logger()
    
    
def generate_input_multiple(): 
    tic = time.time()
    
    # Initialize logger to print nothing
    set_logger(0, 0)
    
    # Initial info of the script
    now = datetime.datetime.now()
    string = 'generate_input_multiple\n' + now.strftime('%Y-%m-%d %H:%M:%S\n\n')
    print(string, end='')  
    
    # Set dynamical variables
    total_charge = 0
    Sz = None
    spin_order = 'AFM'
    spin_order_direction = 1
    
    info = 'gauss'
    orb_dot_coef = 0 # 0: gauss, 1: orbitals read from orb_dot_coef file 
                     # >0: overwrites spin_order

    gndot_v0_listoflist = [
        [-32.82],
        [-29.35],
        [-27.01],
        [-25.28],
        [-22.39],
    ]
    gndot_rho_list = [12.5, 15, 17.5, 20, 25]
    gndot_s = 1.40 # unitless
    
    # sigma guess for Gaussian basis
    gauss_sigma_list = [
        [8.35],
        [9.27],
        [10.13],
        [10.97],
        [13.11],
    ]

    # np.linspace(1.3e-4, 6.19e-4, 10) for rho20    
    gndot_k_list = [
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
    ]
    
    run_dir_naming_elements = [
        'rho',
        's',
        'v0',
        'gsigma',
        'k',
        'total-charge',
        'Sz',
        'nctype-edge',
        'info',
    ] 
    
    root_dir = 'run_dir'
    if not os.path.exists(root_dir): 
        os.mkdir(root_dir)
    
    for i, gndot_v0_list in enumerate(gndot_v0_listoflist):
        for gndot_v0 in gndot_v0_list:
            for gauss_sigma in gauss_sigma_list[i]:
                for gndot_k in gndot_k_list[i]:
                    var_dict_dynamic = {
                        'gndot_v0': gndot_v0,
                        'gndot_rho': gndot_rho_list[i],
                        'gndot_s': gndot_s,
                        'gndot_k': gndot_k,
                        'gauss_sigma': gauss_sigma,
                        
                        'info': info,
                        
                        'total_charge': total_charge,
                        'Sz': Sz,
                        'spin_order': spin_order,
                        'spin_order_direction': spin_order_direction,
                        
                        'orb_dot_coef': orb_dot_coef,
                    }
                    
                    inputgenerator = InputGenerator(CONFIG)
                    inputgenerator.__dict__.update(var_dict_dynamic)
                    inputgenerator.set()
                    Sz_calc = 0.5 * (inputgenerator.nup - inputgenerator.ndn)
                    
                    run_dir = ''
                    if 'rho' in run_dir_naming_elements:
                        run_dir += 'rho%s' %gndot_rho_list[i] + '_'
                    if 's' in run_dir_naming_elements:
                        run_dir += 's%s' %gndot_s + '_'
                    if 'v0' in run_dir_naming_elements:
                        run_dir += 'v0%s' %gndot_v0 + '_' 
                    if 'gsigma' in run_dir_naming_elements:
                        run_dir += 'gsigma%s' %gauss_sigma + '_'
                    if 'k' in run_dir_naming_elements:
                        run_dir += 'k%.5e' %gndot_k + '_'
                    if 'total-charge' in run_dir_naming_elements:
                        run_dir += 'total-charge%s' %total_charge + '_'
                    if 'Sz' in run_dir_naming_elements:
                        run_dir += 'Sz%.1f' %Sz_calc + '_'
                    if 'nctype-edge' in run_dir_naming_elements:
                        run_dir += \
                            'nctype-edge%i' %CONFIG['nctype_of_edges'] + '_'
                    if 'info' in run_dir_naming_elements:
                        run_dir += info
                    if run_dir[-1] == '_':
                        run_dir = run_dir[:-1]
                    run_dir = run_dir.replace('+', '')
            
                    os.mkdir(os.path.join(root_dir, run_dir))
                    
                    if orb_dot_coef != 0:
                        if 'gauss' not in info:     
                            copy_destination = os.path.join(
                                root_dir, run_dir, 'orb_dot_coef'
                            )
                            orb_dot_coef_dir = os.path.join(
                                info, 'orb_dot_coef'
                            )
                            shutil.copy2(orb_dot_coef_dir, copy_destination)
                        else:
                            sys.exit('info is wrong: %s' \
                                %os.path.join(root_dir, run_dir))
                    elif orb_dot_coef == 0:
                        if 'gauss' not in info:
                            sys.exit('info is wrong: %s' \
                                %os.path.join(root_dir, run_dir))        
                    
                    inputgenerator.champ_input_file_name = os.path.join(
                       root_dir, run_dir, inputgenerator.champ_input_file_name
                    )
                    inputgenerator.write() 
 
    # End of program
    toc = time.time()
    string = 'Execution time = %.3f s\n' %(toc-tic)
    print(string, end='')
    
    # Close all handlers of logger
    unset_logger()  
    

if __name__ == "__main__":      
    if GENERATE_INPUT_MULTIPLE:
        generate_input_multiple()
    else:
        generate_input()


