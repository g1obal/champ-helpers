"""
CHAMP Plot Optimization Steps

Author: Gokhan Oztarhan
Created date: 10/06/2021
Last modified: 25/05/2024
"""

import os
import sys
from copy import deepcopy
from itertools import cycle, product
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator, FixedLocator, MaxNLocator

from champio.outputparser import OutputParser


GROUPBY = ['rho'] # None for no groupby, list of strings for groupby
FILTER = None # None for no filter, list of strings for filter
to_SI = True

TOT_E = True
STDEV = True
TCORR = True
ACCEPTANCE = True
GAUSS_SIGMA = False

DPI = 200
PLOT_FORMAT = 'png'
FONTSIZE_XYLABEL = 16
LABELSIZE_TICK_PARAMS = 12
FONTSIZE_LEGEND = 10

ROOT_DIR = '.'
OUTPUT_FILE_NAME = 'output_file'
PLOT_DIR = 'data_plots_opt'
if GROUPBY is not None:
    PLOT_DIR = '_'.join([PLOT_DIR] + GROUPBY)
if FILTER is not None:
    PLOT_DIR = '_'.join([PLOT_DIR] + FILTER)
    
# Matching directories 
PATH = {}
for root, dirs, files in sorted(os.walk(ROOT_DIR)):
    if OUTPUT_FILE_NAME in files:
        run_dir = os.path.split(root)[-1]
        if FILTER is not None:
            if not all([s in run_dir for s in FILTER]):
                continue
        if GROUPBY is not None:
            run_dir_split = run_dir.split('_')
            s = []
            for group in GROUPBY:
                ind = [i for i, s in enumerate(run_dir_split) if group in s][0]
                s.append(run_dir_split[ind])
            s = '_'.join(s)
            if FILTER is not None:
                s = '_'.join([s] + FILTER)
            if s in PATH:
                PATH[s].append(root)
            else:
                PATH[s] = [root]
        else:
            if 'all' in PATH:
                PATH['all'].append(root)
            else:
                PATH['all'] = [root]

warnings.filterwarnings('ignore') # Suppress warnings!

plt.ioff()
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
    'text.latex.preamble': r'\usepackage{amsmath}'
})

if not os.path.exists(PLOT_DIR):
    os.mkdir(PLOT_DIR)


def plot_opt_steps(): 
    for groupname in PATH:
        parser_list = []
        for path in PATH[groupname]:
            parser_list.append(OutputParser(OUTPUT_FILE_NAME, path=path))
            parser_list[-1].parse()
            if to_SI:
                parser_list[-1].to_SI()
        
        if TOT_E:
            plot(
                groupname, 
                parser_list, 'tot_E_all', errfeature='tot_E_err_all',
                ylabel='Total Energy',
                title='group = %s, to\_SI = %s' %(groupname, to_SI),
            )
        if STDEV:
            plot(
                groupname, 
                parser_list, 'stdev_all', errfeature='stdev_err_all',
                ylabel='Standard Deviation of Total Energy',
                title='group = %s, to\_SI = %s' %(groupname, to_SI),
            )
        if TCORR:
            plot(
                groupname, 
                parser_list, 'Tcorr_all',
                ylabel='Tcorr',
                title='group = %s, to\_SI = %s' %(groupname, to_SI),
            )
        if ACCEPTANCE:
            plot(
                groupname, 
                parser_list, 'acceptance_all',
                ylabel='Acceptance',
                title='group = %s, to\_SI = %s' %(groupname, to_SI),
            )
        if GAUSS_SIGMA:
            plot(
                groupname, 
                parser_list, 'gauss_sigma_all',
                ylabel='Gauss Sigma (Au)',
                title='group = %s, to\_SI = %s' %(groupname, to_SI),
            )
            

def plot(
    groupname,
    parser_list, 
    yfeature, errfeature=None, 
    ylabel=None,
    ytype=float,
    markersize=3,
    elinewidth=1, 
    capsize=2,
    title=None,
):
    if ylabel is None:
        ylabel = yfeature.replace('_', '\_')

    # Initialize line style cycler
    line_cycler = _linestyle()
    
    # Initialize figure
    fig = plt.figure(figsize=plt.figaspect(1.0))
    
    # Initialize ax list
    ax = []
    ax.append(fig.add_subplot(1, 1, 1))
    
    # ax[0]
    for parser in parser_list:
        style = next(line_cycler)
        linestyle = style['linestyle']
        color = style['color']
        legend = parser.run_dir.replace('_', '\_')
    
        y = parser.__dict__[yfeature]
        x = range(y.shape[0])
    
        if errfeature is not None:
            yerr = parser.__dict__[errfeature]
        else:
            yerr = None

        ax[-1].errorbar(
            x, y, yerr, label=legend,
            color=color, linestyle=linestyle[1], 
            fmt=linestyle[0], markersize=markersize,
            elinewidth=elinewidth, capsize=capsize
        )

    ax[-1].xaxis.set_major_locator(MaxNLocator(integer=True)) 
    ax[-1].tick_params(labelsize=LABELSIZE_TICK_PARAMS)
    
    ax[-1].set_xlabel('Optimization Iteration', fontsize=FONTSIZE_XYLABEL)
    ax[-1].set_ylabel(ylabel.replace('_', '\_'), fontsize=FONTSIZE_XYLABEL)
    
    ax[-1].legend(
        bbox_to_anchor =(1.0, 1.0), ncol = 1, fontsize=FONTSIZE_LEGEND
    )
    ax[-1].set_title(title)

    filename = os.path.join(
        PLOT_DIR, yfeature + '_' + groupname + '.' + PLOT_FORMAT
    )
    
    fig.savefig(filename, dpi=DPI, format=PLOT_FORMAT, bbox_inches='tight') 
    plt.close(fig)
    
    print('Done plot: %s %s' %(groupname, yfeature))

    
def _linestyle(): 
    linestyle = { 
        'linestyle': [['o', '-'], ['^', '--'], ['s', '-.'], ['h', ':']], 
        'color': ['k', 'r', 'y', 'g', 'c', 'b', 'm'],
    }    
    linestyle = _grid(linestyle)
    return cycle(linestyle)


def _grid(params):
    return [dict(zip(params, i)) for i in product(*params.values())]
    
            
if __name__ == '__main__':
    plot_opt_steps()
        
 
