"""
Polyfit V0 Data

Author: Gokhan Oztarhan
Created date: 02/02/2023
Last modified: 07/02/2023
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


TOT_E_CUTOFF = 50
POLY_ORDER = 1
FIND_X_AT_Y_EQUALS = -0.15
BACKGROUND_LINE_Y = 0e0

DATA_FILE_NAME = 'data_SI.csv'
SORT_COLS = ['gndot_rho', 'info']
GROUPBY_COLS = ['gndot_rho']

PLOT_DIR = 'data_plots'

DPI = 200
PLOT_FORMAT = 'png'
FONTSIZE_XYLABEL = 16
LABELSIZE_TICK_PARAMS = 12
FONTSIZE_LEGEND = 10
DRAW_LEGEND = True

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


def polyfit_v0_data():
    # Load data and sort
    df = pd.read_csv(DATA_FILE_NAME, header=0, index_col=0)
    df = df.sort_values(SORT_COLS, ignore_index=True)
    
    # Filter data
    df = df[df['tot_E'] < TOT_E_CUTOFF]
    df = df[df['tot_E'] > -TOT_E_CUTOFF]
    
    # Group data
    df_grouped = df.groupby(GROUPBY_COLS)

    # Fit and plot
    for groupname in df_grouped.groups:
        df_rho = df_grouped.get_group(groupname)
        fit_plot(
            df_rho, 'gndot_v0', 'tot_E', errfeature='tot_E_err',
            xlabel='$V_{0}$ (meV)', ylabel='Total Energy (meV)',
            groupname=groupname,
            fname='rho' + str(groupname)[:5]
        )


def fit_plot(
    df_rho, 
    xfeature, yfeature, errfeature=None, 
    xlabel=None, ylabel=None, 
    xtype=float, ytype=float,
    markersize=3,
    elinewidth=1, 
    capsize=2,
    groupname=None,
    fname=None
):     
    if xlabel is None:
        xlabel = xfeature.replace('_', '\_')
    if ylabel is None:
        ylabel = yfeature.replace('_', '\_')
    
    # Initialize figure
    fig = plt.figure(figsize=plt.figaspect(1.0))
    
    # Initialize ax list
    ax = []
    ax.append(fig.add_subplot(1, 1, 1))
    
    #if isinstance(groupname, tuple):
    #    groupname = ', '.join(groupname).replace('_', '\_')
    legend = '$\\rho$ = ' + str(groupname)[:4]

    x = get_feature(df_rho, xfeature, xtype)
    y = get_feature(df_rho, yfeature, ytype)
    
    if errfeature is not None:
        yerr = df_rho[errfeature]
        # If all elements of yerr are NaN, an error is raised.
        if np.isnan(yerr).all():
            yerr = None
    else:
        yerr = None
    
    ax[-1].errorbar(
        x, y, yerr, label=legend,
        color='k', linestyle='none', 
        fmt='o', markersize=markersize,
        elinewidth=elinewidth, capsize=capsize
    )
    
    # Draw background line
    x_bg = np.linspace(x.min(), x.max(), 10)
    y_bg = np.full(x_bg.size, BACKGROUND_LINE_Y)
    ax[-1].plot(x_bg, y_bg, '--k', alpha=0.5)
    
    # Polynomial fit
    poly = np.polyfit(x.to_numpy(), y.to_numpy(), POLY_ORDER)
    x_poly = np.linspace(x.min(), x.max(), 100)
    y_poly = np.polyval(poly, x_poly)
    ax[-1].plot(x_poly, y_poly, '--r', alpha=0.35, label='polyfit')
    
    # Find the x value at which y value is around FIND_X_AT_Y_EQUALS.
    x_finder = np.linspace(x.min(), x.max(), 100000)
    y_finder = np.polyval(poly, x_finder)
    argmin = np.argmin(np.abs(y_finder - FIND_X_AT_Y_EQUALS))
    ax[-1].scatter(x_finder[argmin], y_finder[argmin], color='g', alpha=1)

    title = 'x = %.8f, y = %.8f' %(x_finder[argmin], y_finder[argmin])
    
    # x-axis locators
    if xtype == str:
        ax[-1].tick_params(axis='x', rotation=90)
    else:
        #ax[-1].xaxis.set_major_locator(FixedLocator(x)) 
        #ax[-1].xaxis.set_major_locator(FixedLocator([10, 15, 20, 25, 30, 35])) 
        ax[-1].xaxis.set_major_locator(LinearLocator(5)) 
        #ax[-1].xaxis.set_major_locator(MaxNLocator(nbins=6)) 
        #ax[-1].tick_params(axis='x', rotation=45)
        #ax[-1].ticklabel_format(axis='x', style='scientific', scilimits=(0,0))
        #ax[-1].tick_params(axis='x', labelsize=LABELSIZE_TICK_PARAMS)
    
    # y-axis locators
    #ax[-1].yaxis.set_major_locator(LinearLocator(5))
    #ax[-1].yaxis.set_major_locator(MaxNLocator(nbins=6)) 
    #ax[-1].tick_params(axis='y', labelsize=LABELSIZE_TICK_PARAMS)
    
    ax[-1].set_xlim([x.min(), x.max()])
    
    ax[-1].tick_params(labelsize=LABELSIZE_TICK_PARAMS)
    
    ax[-1].set_xlabel(xlabel, fontsize=FONTSIZE_XYLABEL)
    ax[-1].set_ylabel(ylabel, fontsize=FONTSIZE_XYLABEL)
        
    if DRAW_LEGEND:
        ax[-1].legend(loc='best', fontsize=FONTSIZE_LEGEND)
    
    ax[-1].set_title(title)
    
    if fname is None:
        filename = os.path.join(PLOT_DIR, yfeature + '.' + PLOT_FORMAT)
    else:
        filename = os.path.join(PLOT_DIR, fname + '.' + PLOT_FORMAT)
        
    fig.savefig(filename, dpi=DPI, format=PLOT_FORMAT, bbox_inches='tight') 
    plt.close(fig)
    
    print('Done fit_plot: %s' %filename)

    
def get_feature(df, feature, _type):
    if _type == str:
        x = df[feature].apply(lambda string: string.replace('_', '\_'))
    else:
        x = df[feature]
    return x


if __name__ == '__main__':
    polyfit_v0_data()
        

