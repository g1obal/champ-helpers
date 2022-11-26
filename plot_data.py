"""
Plot data

Author: Gokhan Oztarhan
Created date: 12/10/2021
Last modified: 20/11/2022
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


DATA_FILE_NAME = 'data_SI.csv'
SORT_COLS = ['gndot_rho', 'info']

GROUPBY_COLS = ['info']
GS_GROUPBY_COLS = ['gndot_rho'] # for ground state plot

# Drop rows if run_dir or info columns include the desired string
DROP_ROWS = []

PLOT_DIR = 'data_plots'

PLOT_V0 = False
PLOT_SS_CORR = True
PLOT_EDGE_POL = True

DPI = 200
PLOT_FORMAT = 'png'
FONTSIZE_XYLABEL = 16
LABELSIZE_TICK_PARAMS = 12
FONTSIZE_LEGEND = 10
DRAW_LEGEND = True
DRAW_TITLE = False

# Set LEGEND_ORDER to None for plotting all groups regardless of order,
# as well as the groups that are not in this list
LEGEND_ORDER = [
    'tb',
    'mfh_U_t',
    'mfh_U_2t',
    'mfh_U_3t',
    'mfh_U_5t',
    'mfh_U_20t',
    'mfh_U_100t',
    'dft',
    'gauss',
    'T = 4 K',
]
LEGEND = {
    'tb': 'Tight-binding',
    'mfh_U_t': 'Hubbard $U = t$',
    'mfh_U_2t': 'Hubbard $U = 2t$',
    'mfh_U_3t': 'Hubbard $U = 3t$',
    'mfh_U_5t': 'Hubbard $U = 5t$',
    'mfh_U_20t': 'Hubbard $U = 20t$',
    'mfh_U_100t': 'Hubbard $U = 100t$',
    'dft': 'DFT',
    'gauss': 'Gaussian',
    'T = 4 K': 'T = 4 K',
}

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
   

def plot_data():
    # Load data and sort
    df = pd.read_csv(DATA_FILE_NAME, header=0, index_col=0)
    df = df.sort_values(SORT_COLS, ignore_index=True)
    
    # Drop rows
    for value in DROP_ROWS:
        mask = df[['run_dir', 'info']].apply(
            lambda s: s.str.contains(value)
        ).any(axis=1)
        df = df[~mask]
    
    # Check whether all runs are DMC
    run_mode_dmc = (df['run_mode'] == 'dmc').all(axis=0)
    
    # Check whether gndot_k is equal to zero for all runs
    gndot_k_zero = (df['gndot_k'] == 0).all(axis=0)
    
    # Calculate some features
    df['Tcorr_n_steps'] = df['Tcorr'] / df['n_steps']
    df['tot_E_shifted'] = df['tot_E'] - df['eshift']
    df['wts_n_walkers_global'] = df['wts'] / df['n_walkers_global']
    df['wts_f_n_walkers_global'] = df['wts_f'] / df['n_walkers_global']
    
    # Group data      
    df_grouped = df.groupby(GROUPBY_COLS)
  
    plot(
        df_grouped, 'gndot_rho', 'tot_E', errfeature='tot_E_err',
        xlabel='$\\rho$ (nm)', ylabel='Total Energy (meV)',
        title='Total Energy vs $\\rho$',
        groups=LEGEND_ORDER,
    )
    plot(
        df_grouped, 'gndot_rho', 'stdev', errfeature='stdev_err',
        xlabel='$\\rho$ (nm)', ylabel='Standard Deviation of Total Energy',
        title='Standard Deviation of Total Energy vs $\\rho$',
        groups=LEGEND_ORDER,
    )
    plot(
        df_grouped, 'gndot_rho', 'acceptance',
        xlabel='$\\rho$ (nm)', ylabel='Acceptance',
        title='Acceptance vs $\\rho$',
        groups=LEGEND_ORDER,
    )
    plot(
        df_grouped, 'gndot_rho', 'Tcorr_n_steps',
        xlabel='$\\rho$ (nm)', ylabel='Tcorr / n\_steps',
        title='Tcorr / n\_steps Ratio vs $\\rho$',
        groups=LEGEND_ORDER,
    )
    
    if PLOT_SS_CORR:
        plot(
            df_grouped, 'gndot_rho', 'ss_corr_den',
            xlabel='$\\rho$ (nm)', ylabel='$g$ (s.s. corr. - density)',
            title='Spin-spin Corr. for Density vs $\\rho$',
            groups=LEGEND_ORDER,
        )
        plot(
            df_grouped, 'gndot_rho', 'ss_corr_pairden',
            xlabel='$\\rho$ (nm)', ylabel='$g$ (s.s. corr. - pair density)',
            title='Spin-spin Corr. for Pair Density vs $\\rho$',
            groups=LEGEND_ORDER,
        )
        
    if PLOT_EDGE_POL:
        plot(
            df_grouped, 'gndot_rho', 'edge_pol_den',
            xlabel='$\\rho$ (nm)', ylabel='$p$ (edge pol. - density)',
            title='Edge Polarization for Density vs $\\rho$',
            groups=LEGEND_ORDER,
        )
        plot(
            df_grouped, 'gndot_rho', 'edge_pol_pairden',
            xlabel='$\\rho$ (nm)', ylabel='$p$ (edge pol. - pair density)',
            title='Edge Polarization for Pair Density vs $\\rho$',
            groups=LEGEND_ORDER,
        )
    
    if run_mode_dmc:
        plot(
            df_grouped, 'gndot_rho', 'tot_E', errfeature='pop_err',
            xlabel='$\\rho$ (nm)', ylabel='Total Energy (meV)',
            title='Total Energy with Population Error vs $\\rho$',
            groups=LEGEND_ORDER,
            fname='tot_E_pop_err'
        )
        plot(
            df_grouped, 'gndot_rho', 'overlap_dmc',
            xlabel='$\\rho$ (nm)', ylabel='Approx. Normalized Overlap',
            title='Approx. Normalized Overlap vs $\\rho$',
            groups=LEGEND_ORDER,
            fname='dmc_overlap'
        )
        plot(
            df_grouped, 'gndot_rho', 'nwalk_eff',
            xlabel='$\\rho$ (nm)', 
            ylabel='nwalk\_eff',
            title='nwalk\_eff vs $\\rho$',
            groups=LEGEND_ORDER,
            fname='dmc_nwalk_eff'
        )
        plot(
            df_grouped, 'gndot_rho', 'nwalk_eff_f',
            xlabel='$\\rho$ (nm)', 
            ylabel='nwalk\_eff\_f',
            title='nwalk\_eff\_f vs $\\rho$',
            groups=LEGEND_ORDER,
            fname='dmc_nwalk_eff_f'
        )
        plot(
            df_grouped, 'gndot_rho', 'nwalk_eff_fs',
            xlabel='$\\rho$ (nm)', 
            ylabel='nwalk\_eff\_fs',
            title='nwalk\_eff\_fs vs $\\rho$',
            groups=LEGEND_ORDER,
            fname='dmc_nwalk_eff_fs'
        )
        plot(
            df_grouped, 'gndot_rho', 'wts_n_walkers_global',
            xlabel='$\\rho$ (nm)', 
            ylabel='wts / n\_walkers\_global',
            title='wts / n\_walkers\_global vs $\\rho$',
            groups=LEGEND_ORDER,
            fname='dmc_wts'
        )
        plot(
            df_grouped, 'gndot_rho', 'wts_f_n_walkers_global',
            xlabel='$\\rho$ (nm)', 
            ylabel='wts\_f / n\_walkers\_global',
            title='wts\_f / n\_walkers\_global vs $\\rho$',
            groups=LEGEND_ORDER,
            fname='dmc_wts_f'
        )
    
    if not gndot_k_zero:
        plot(
            df_grouped, 'gndot_rho', 'tot_E_shifted', errfeature='tot_E_err',
            xlabel='$\\rho$ (nm)', ylabel='Total Energy Shifted (meV)',
            title='Total Energy Shifted vs $\\rho$',
            groups=LEGEND_ORDER,
        )
        
    if run_mode_dmc and not gndot_k_zero:
        plot(
            df_grouped, 'gndot_rho', 'tot_E_shifted', errfeature='pop_err',
            xlabel='$\\rho$ (nm)', ylabel='Total Energy Shifted (meV)',
            title='Total Energy Shifted vs $\\rho$',
            groups=LEGEND_ORDER,
            fname='tot_E_shifted_pop_err'
        )
        
    if PLOT_V0:
        plot(df_grouped, 'gndot_rho', 'gndot_v0', groups=LEGEND_ORDER) 
    
    # Ground states for all rho values
    df_gs = df.loc[df.groupby(GS_GROUPBY_COLS)['tot_E'].idxmin()]       
    df_grouped = df_gs.groupby(GROUPBY_COLS)
    if PLOT_SS_CORR:
        plot(
            df_grouped, 'gndot_rho', 'ss_corr_den',
            xlabel='$\\rho$ (nm)', ylabel='$g$ (s.s. corr. - density - GS)',
            title='Spin-spin Corr. for Density vs $\\rho$',
            groups=LEGEND_ORDER,
            connect_group_points=False,
            markersize=6,
            df_background_line=df_gs,
            fname='gs_ss_corr_den'
        )
        plot(
            df_grouped, 'gndot_rho', 'ss_corr_pairden',
            xlabel='$\\rho$ (nm)', 
            ylabel='$g$ (s.s. corr. - pair density - GS)',
            title='Spin-spin Corr. for Pair Density vs $\\rho$',
            groups=LEGEND_ORDER,
            connect_group_points=False,
            markersize=6,
            df_background_line=df_gs,
            fname='gs_ss_corr_pairden'
        )
        
    if PLOT_EDGE_POL:
        plot(
            df_grouped, 'gndot_rho', 'edge_pol_den',
            xlabel='$\\rho$ (nm)', ylabel='$p$ (edge pol. - density - GS)',
            title='Edge Polarization for Density vs $\\rho$',
            groups=LEGEND_ORDER,
            connect_group_points=False,
            markersize=6,
            df_background_line=df_gs,
            fname='gs_edge_pol_den'
        )
        plot(
            df_grouped, 'gndot_rho', 'edge_pol_pairden',
            xlabel='$\\rho$ (nm)', 
            ylabel='$p$ (edge pol. - pair density - GS)',
            title='Edge Polarization for Pair Density vs $\\rho$',
            groups=LEGEND_ORDER,
            connect_group_points=False,
            markersize=6,
            df_background_line=df_gs,
            fname='gs_edge_pol_pairden'
        )


def plot(
    df_grouped, 
    xfeature, yfeature, errfeature=None, 
    xlabel=None, ylabel=None, 
    xtype=float, ytype=float,
    title=None,
    groups=None,
    connect_group_points=True,
    markersize=3,
    elinewidth=1, 
    capsize=2,
    df_background_line=None,
    fname=None
):     
    if xlabel is None:
        xlabel = xfeature.replace('_', '\_')
    if ylabel is None:
        ylabel = yfeature.replace('_', '\_')
        
    if groups is not None:
        groups = [
            groupname for groupname in groups \
            if groupname in df_grouped.groups
        ]
    else:
        groups = df_grouped.groups

    # Initialize line style cycler
    line_cycler = _linestyle()
    
    # Initialize figure
    fig = plt.figure(figsize=plt.figaspect(1.0))
    
    # Initialize ax list
    ax = []
    ax.append(fig.add_subplot(1, 1, 1))
    
    # Draw background line
    if df_background_line is not None:
        x = get_feature(df_background_line, xfeature, xtype)
        y = get_feature(df_background_line, yfeature, ytype)
        ax[-1].plot(x, y, '--k', alpha=0.5)
    
    for groupname in groups:     
        df = df_grouped.get_group(groupname)
    
        style = next(line_cycler)
        linestyle = style['linestyle']
        color = style['color']
        legend = LEGEND.get(groupname, groupname)
        if isinstance(legend, tuple):
            legend = ', '.join(legend).replace('_', '\_')
        
        x = get_feature(df, xfeature, xtype)
        y = get_feature(df, yfeature, ytype)
            
        if errfeature is not None:
            yerr = df[errfeature]
        else:
            yerr = None
            
        if not connect_group_points:
            linestyle[1] = 'none'
        
        ax[-1].errorbar(
            x, y, yerr, label=legend,
            color=color, linestyle=linestyle[1], 
            fmt=linestyle[0], markersize=markersize,
            elinewidth=elinewidth, capsize=capsize
        )
    
    # Get all x values if group points are not connected
    if df_background_line is not None and not connect_group_points:
        x = get_feature(df_background_line, xfeature, xtype)

    # x-axis locators
    if xtype == str:
        ax[-1].tick_params(axis='x', rotation=90)
    else:
        ax[-1].xaxis.set_major_locator(FixedLocator(x)) 
        #ax[-1].xaxis.set_major_locator(FixedLocator([10, 15, 20, 25, 30, 35])) 
        #ax[-1].xaxis.set_major_locator(LinearLocator(5)) 
        #ax[-1].xaxis.set_major_locator(MaxNLocator(nbins=6)) 
        #ax[-1].tick_params(axis='x', rotation=45)
        #ax[-1].ticklabel_format(axis='x', style='scientific',scilimits=(0,0))
        #ax[-1].tick_params(axis='x', labelsize=LABELSIZE_TICK_PARAMS)
    
    # y-axis locators
    #ax[-1].yaxis.set_major_locator(LinearLocator(5))
    #ax[-1].yaxis.set_major_locator(MaxNLocator(nbins=6)) 
    #ax[-1].tick_params(axis='y', labelsize=LABELSIZE_TICK_PARAMS)
    
    ax[-1].tick_params(labelsize=LABELSIZE_TICK_PARAMS)
    
    ax[-1].set_xlabel(xlabel, fontsize=FONTSIZE_XYLABEL)
    ax[-1].set_ylabel(ylabel, fontsize=FONTSIZE_XYLABEL)
        
    if DRAW_LEGEND:
        ax[-1].legend(
            loc='best',
            bbox_to_anchor =(1.0, 1.0), 
            ncol = 1, 
            fontsize=FONTSIZE_LEGEND
        )
        
    if DRAW_TITLE and title is not None:
        ax[-1].set_title(title)
    
    if fname is None:
        filename = os.path.join(PLOT_DIR, yfeature + '.' + PLOT_FORMAT)
    else:
        filename = os.path.join(PLOT_DIR, fname + '.' + PLOT_FORMAT)
        
    fig.savefig(filename, dpi=DPI, format=PLOT_FORMAT, bbox_inches='tight') 
    plt.close(fig)
    
    print('Done plot: %s' %filename)
    
    
def get_feature(df, feature, _type):
    if _type == str:
        x = df[feature].apply(lambda string: string.replace('_', '\_'))
    else:
        x = df[feature]
    return x
            

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
    plot_data()
        

