"""
Plot energy data for all x values on different subplots

Author: Gokhan Oztarhan
Created date: 12/04/2024
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


DATA_FILE_NAME = 'data_SI.csv'

X_AXIS = 'gndot_v0'
X_LABEL = '$V_{0}$ (meV)'
X_TITLE = '$V_{0}$'

SORT_COLS = ['gndot_k', 'gndot_v0', 'info']

GROUPBY_COLS = ['gndot_k', 'gndot_v0']
GROUPBY_COLS_LABEL = ['$k$', '$V_{0}$']

GROUPBY_COLS_INSIDE = ['info']

# Drop rows if run_dir or info columns include the desired string
DROP_ROWS = []

PLOT_DIR = 'data_plots_energy'

DPI = 600
PLOT_FORMAT = 'png'
FONTSIZE_XYLABEL = 10
LABELSIZE_TICK_PARAMS = 10
FONTSIZE_LEGEND = 10
FONTSIZE_TEXT = 10
SCALE_FACTOR = 0.0045

# Set LEGEND_ORDER to None for plotting all groups regardless of order,
# as well as the groups that are not in this list
LEGEND_ORDER = None
LEGEND = {
    'tb': 'TB',
    'mfh_U_t': 'MFH $U = t$',
    'mfh_U_2t': 'MFH $U = 2t$',
    'mfh_U_3t': 'MFH $U = 3t$',
    'mfh_U_5t': 'MFH $U = 5t$',
    'mfh_U_20t': 'MFH $U = 20t$',
    'mfh_U_100t': 'MFH $U = 100t$',
    'dft': 'DFT',
    'gauss': 'Gauss',
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
    
    # Calculate some features
    df['Tcorr_n_steps'] = df['Tcorr'] / df['n_steps']
    df['tot_E_shifted'] = df['tot_E'] - df['eshift']
    df['wts_n_walkers_global'] = df['wts'] / df['n_walkers_global']
    df['wts_f_n_walkers_global'] = df['wts_f'] / df['n_walkers_global']
    df['total_charge'] = df['ncent'] - df['nelec']
    df['Sz'] = 0.5 * (df['nup'] - df['ndn'])
    
    # Sort values
    df = df.sort_values(SORT_COLS, ignore_index=True)
    
    # Drop rows
    for value in DROP_ROWS:
        mask = df[['run_dir', 'info']].apply(
            lambda s: s.str.contains(value)
        ).any(axis=1)
        df = df[~mask]
    
    # Check whether all runs are DMC
    run_mode_dmc = (df['run_mode'] == 'dmc').all(axis=0)
    
    # Group data      
    df_grouped = df.groupby(GROUPBY_COLS)
  
    plot(
        df_grouped, X_AXIS, 'tot_E', errfeature='tot_E_err', 
        errfeature_dmc='pop_err',
        xlabel=X_LABEL, ylabel='Total Energy (meV)',
        title='Total Energy vs %s' %X_TITLE,
        groups=LEGEND_ORDER,
        run_mode_dmc=run_mode_dmc,
    )


def plot(
    df_grouped, 
    xfeature, yfeature, errfeature=None, errfeature_dmc=None,
    xlabel=None, ylabel=None, 
    xtype=float, ytype=float,
    title=None,
    groups=None,
    markersize=3,
    elinewidth=1, 
    capsize=3,
    capthick=1,
    run_mode_dmc=False,
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
    
    # Initialize ax grid values
    n_tot_ax = len(groups)
    n_ax_x = int(np.rint(np.sqrt(n_tot_ax)))
    n_ax_y = int(np.rint(np.sqrt(n_tot_ax)))
    if n_ax_x * n_ax_y < n_tot_ax:
        n_ax_y += 1
        
    # Initialize figure
    fig = plt.figure(figsize=plt.figaspect(n_ax_x / n_ax_y))
    
    # Initialize ax list
    ax = []
    
    i_ax = 1
    for groupname in groups: 
        # Add ax
        ax.append(fig.add_subplot(n_ax_x, n_ax_y, i_ax))
        i_ax += 1
        
        # Scale multiplier for fonts, markers and linewidths
        mult = np.sqrt(ax[-1].bbox.height * ax[-1].bbox.width) * SCALE_FACTOR
    
        df_grouped_inside = \
            df_grouped.get_group(groupname).groupby(GROUPBY_COLS_INSIDE)
            
        # Initialize line style cycler
        line_cycler = _linestyle()
        
        for groupname_inside in df_grouped_inside.groups:
            df = df_grouped_inside.get_group(groupname_inside)
        
            style = next(line_cycler)
            linestyle = style['linestyle']
            color = style['color']
            legend = LEGEND.get(groupname_inside, groupname_inside)
            if isinstance(legend, tuple):
                legend = ', '.join(str(lgn) for lgn in legend)\
                    .replace('_', '\_')
            
            x = get_feature(df, xfeature, xtype)
            y = get_feature(df, yfeature, ytype)
            
            if errfeature is not None:
                yerr = df[errfeature]
                # If all elements of yerr are NaN, an error is raised.
                if np.isnan(yerr).all():
                    yerr = None
                    
                if run_mode_dmc and errfeature_dmc is not None:
                    yerr_dmc = df[errfeature_dmc]
                    if np.isnan(yerr_dmc).all():
                        yerr_dmc = None
            else:
                yerr = None
                yerr_dmc = None
            
            if run_mode_dmc and errfeature_dmc is not None:
                ax[-1].errorbar(
                    x, y, yerr_dmc,
                    color='gray', linestyle='', 
                    fmt='', markersize=0,
                    elinewidth=elinewidth*mult, capsize=capsize*mult,
                    capthick=capthick*mult
                )
            
            ax[-1].errorbar(
                x, y, yerr, label=legend,
                color=color, linestyle='', 
                fmt=linestyle[0], markersize=markersize*mult,
                elinewidth=elinewidth*mult, capsize=capsize*mult,
                capthick=capthick*mult
            )
            
            ax[-1].text(
                x + abs(x) * 0.05, y, str(y.iloc[0]), 
                fontsize=FONTSIZE_TEXT*mult,
                verticalalignment='center', horizontalalignment='left',
                bbox=dict(facecolor='none', edgecolor='none', pad=0.5)
            )

            # x-axis locators
            ax[-1].set_xticks([])
            
            # y-axis locators
            ax[-1].set_yticks([])
            
            #ax[-1].tick_params(labelsize=LABELSIZE_TICK_PARAMS / font_mult)
            
            # limits
            xpoint = x.iloc[0]
            xpoint_abs = abs(xpoint)
            if xpoint_abs == 0:
                xpoint_abs = 1
            ax[-1].set_xlim(
                [xpoint - 0.10 * xpoint_abs, xpoint + 0.90 * xpoint_abs]
            )

            #ax[-1].set_xlabel(xlabel, fontsize=FONTSIZE_XYLABEL / font_mult)
            #ax[-1].set_ylabel(ylabel, fontsize=FONTSIZE_XYLABEL / font_mult)
            
            # Legend
            handles, labels = ax[-1].get_legend_handles_labels()
            handles = [hs[0] for hs in handles] # remove error bars
            legend = ax[-1].legend(
                handles, labels,
                loc='best',
                bbox_to_anchor =(1.0, 1.0), 
                ncol = 1, 
                fontsize=FONTSIZE_LEGEND * mult,
                handletextpad=0.1,
            )
            legend.get_frame().set_linewidth(0.375)
        
        # Group info text
        group_info = LEGEND.get(groupname, groupname)
        if isinstance(group_info, tuple):
            group_info = '\n'.join(
                GROUPBY_COLS_LABEL[it] + '=' + str(lgn).replace('_', '\_') \
                for it, lgn in enumerate(group_info)
            )
        else:
            group_info = \
                GROUPBY_COLS_LABEL[0] + '=' + str(group_info).replace('_', '\_')
        
        ax[-1].text(
            0.99, 0.01, group_info, transform=ax[-1].transAxes,
            fontsize=FONTSIZE_TEXT * mult,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(facecolor='none', edgecolor='none', pad=0.5)
        )
        
    # Axes border width
    for _ax in ax:
        for axis in ['top','bottom','left','right']:
            _ax.spines[axis].set_linewidth(0.375)
            _ax.tick_params(width=0.375, length=1.6875)
            
    # Space between subplots
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0, hspace=0)
    
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
        

