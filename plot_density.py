"""
Plot Density

Author: Gokhan Oztarhan
Created date: 18/01/2022
Last modified: 03/06/2024
"""

import os
import sys
import time
import pickle
import multiprocessing

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

from champio.outputparser import OutputParser


N_CPU = 4

DEN2D = True
PAIRDEN = True
PAIRDEN_ALL = True

DRAW_SITE_POS = False

DRAW_TICKS = True
DRAW_TITLE = True

SHADING = 'nearest' # nearest or gouraud
PLOT_FORMAT = 'png'

DPI_DEN2D = 200
FONTSIZE_XYLABEL_DEN2D = 14
LABELSIZE_TICK_PARAMS_DEN2D = 10
FONTSIZE_TITLE_DEN2D = 14
FONTSIZE_INFO_DEN2D = 12

DPI_PAIRDEN = 200
FONTSIZE_XYLABEL_PAIRDEN = 14
LABELSIZE_TICK_PARAMS_PAIRDEN = 10
FONTSIZE_TITLE_PAIRDEN = 14
FONTSIZE_INFO_PAIRDEN = 12

DPI_PAIRDEN_ALL = 200
LABELSIZE_TICK_PARAMS_PAIRDEN_ALL = 6
FONTSIZE_TITLE_PAIRDEN_ALL = 10
FONTSIZE_INFO_PAIRDEN_ALL = 10

ROOT_DIR = '.'
PLOT_DIR = 'density_plots'
PAIRDEN_ALL_DIR = os.path.join(PLOT_DIR, 'pairden_all')

OUTPUT_FILE_NAME = 'output_file'
DATA_MEAN_FILE_NAME = 'mean_file.pkl'
DATA_EXT_FILE_NAME = 'ext_file.pkl'
DATA_ENS_FILE_NAME = 'ens_file.pkl'

plt.ioff()
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
    'text.latex.preamble': r'\usepackage{amsmath}'
})


def plot_density():
    tic = time.time()
    
    multiprocessing.set_start_method('spawn')
    
    if not os.path.exists(PLOT_DIR):
        os.mkdir(PLOT_DIR)
    
    if PAIRDEN_ALL:
        if not os.path.exists(PAIRDEN_ALL_DIR):
            os.mkdir(PAIRDEN_ALL_DIR)

    paths = []
    for root, dirs, files in sorted(os.walk(ROOT_DIR)):
        if not dirs:
            for f in files:
                if f in [OUTPUT_FILE_NAME, DATA_MEAN_FILE_NAME, 
                DATA_EXT_FILE_NAME, DATA_ENS_FILE_NAME]:
                    paths.append([root, f])

    chunks = get_chunks(len(paths), N_CPU)

    # Divide chunks to processes one by one (chunksize=1)
    # since I already calculated and distributed chunks 
    # to processes (in order to preserve list order).
    myPool = multiprocessing.Pool()
    result = myPool.map(
        parse_and_plot, [paths[i:j] for [i, j] in chunks], chunksize=1
    )
                         
    toc = time.time() 
    print("Execution time, plot_density = %.3f s" %(toc-tic))


def parse_and_plot(paths):
    for root, out_file in paths:
        try:
            if out_file == OUTPUT_FILE_NAME:
                # Parse output file for required variables
                outparser = OutputParser(
                    out_file, path=root,
                    parse_density=True, 
                    calculate_ss_corr=False,
                    calculate_edge_pol=False,
                    calculate_U_onsite=False,
                )
                outparser.parse()
            else:
                if out_file == DATA_MEAN_FILE_NAME:
                    pkl_info = 'dmc2mean'
                elif out_file == DATA_EXT_FILE_NAME:
                    pkl_info = 'extrapolated'
                elif out_file == DATA_ENS_FILE_NAME:
                    pkl_info = 'ensemble'
                
                # Load data
                with open(os.path.join(root, out_file), 'rb') as pkl_file:
                    pkl_data = pickle.load(pkl_file)
                outparser = OutputParser()
                outparser.__dict__.update(pkl_data)
                outparser.run_mode = pkl_info
        except Exception as err:
            print('%s, parse/load data: %s' %(type(err).__name__, root))
        
        if DEN2D:
            try:
                # Raise error if density parser fails
                if outparser.den2d_parse_err is not None:
                    raise(outparser.den2d_parse_err)
                    
                # Load, preprocess and plot
                plot_den2d(
                    outparser.den2d_t, 
                    outparser.den2d_s, 
                    outparser.den2d_nelec_calc,
                    outparser.pos,
                    root,
                    outparser.run_mode
                )
                
                # Print info
                print('Done den2d, %s: %s' \
                    %(outparser.run_mode, os.path.split(root)[-1]))
            except Exception as err:
                print('%s, den2d: %s' %(type(err).__name__, root))
        
        if PAIRDEN:                   
            try:
                # Raise error if density parser fails
                if outparser.pairden_parse_err is not None:
                    raise(outparser.pairden_parse_err)
                
                # Plot
                plot_pairden(
                    outparser.pairden_t, 
                    outparser.pairden_s, 
                    outparser.pairden_xfix,
                    root,
                    outparser.run_mode
                )
            
                # Print info
                print('Done pairden, %s: %s' \
                    %(outparser.run_mode, os.path.split(root)[-1]))
            except Exception as err:
                print('%s, pairden: %s' %(type(err).__name__, root))
        
        if PAIRDEN_ALL:                   
            try:
                # Raise error if density parser fails
                if outparser.pairden_parse_err is not None:
                    raise(outparser.pairden_parse_err)
                
                # Plot
                plot_pairden_all(
                    outparser.pairden_dd, 
                    outparser.pairden_dt, 
                    outparser.pairden_du, 
                    outparser.pairden_ud, 
                    outparser.pairden_ut, 
                    outparser.pairden_uu,
                    outparser.pairden_xfix,
                    root,
                    outparser.run_mode
                )
            
                # Print info
                print('Done pairden_all, %s: %s' \
                    %(outparser.run_mode, os.path.split(root)[-1]))
            except Exception as err:
                print('%s, pairden_all: %s' %(type(err).__name__, root))


def plot_den2d(den2d_t, den2d_s, nelec_calc, pos, path, run_mode):
    # Initialize figure
    fig = plt.figure(figsize=plt.figaspect(0.5))
    
    # Initialize GridSpec
    gs = gridspec.GridSpec(1, 2)
    
    # Initialize ax and surf list
    ax, surf = [], []
    ax.append(fig.add_subplot(gs[0, 0], adjustable='box', aspect=1.0))
    ax.append(fig.add_subplot(gs[0, 1], adjustable='box', aspect=1.0)) 
    surf.append(None)
    surf.append(None)   

    # _ax_plot arguments
    args = [
        [ax[0], den2d_t, 'jet', 'Electron Density'],
        [ax[1], den2d_s, 'seismic', 'Spin Density'],
    ]
    
    # Plot surfaces
    for i, _args in enumerate(args):
        ax[i], surf[i] = _ax_plot(*_args)
        if DRAW_SITE_POS:
            ax[i].scatter(
                pos[:,0], pos[:,1],
                marker='x', color='white',
                s=70, linewidth=2
            )
    
    # Colorbars and fonts
    for i, _ax in enumerate(ax):
        divider = make_axes_locatable(_ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(surf[i], cax=cax, ticks=LinearLocator(5))
        cbar.ax.tick_params(labelsize=LABELSIZE_TICK_PARAMS_DEN2D)
        if DRAW_TICKS:
            _ax.set_xlabel('x', fontsize=FONTSIZE_XYLABEL_DEN2D)
            _ax.set_ylabel('y', fontsize=FONTSIZE_XYLABEL_DEN2D)
            _ax.tick_params(labelsize=LABELSIZE_TICK_PARAMS_DEN2D)
        else:
            _ax.set_xticks([])
            _ax.set_yticks([]) 
        if DRAW_TITLE:
            _ax.set_title(args[i][-1], fontsize=FONTSIZE_TITLE_DEN2D)
        
    # Information ax
    #ax_info = fig.add_subplot(1,1,1, adjustable='box', aspect=1.0)
    #ax_info.clear()
    #ax_info.axis("off")
    info_str = os.path.split(path)[-1].replace('_', '\_') + '\n' \
               + '[%s], nelec\_calc = %.2f' %(run_mode, nelec_calc)
    fig.text(
        0.5, 0.90, info_str, ha='center', fontsize=FONTSIZE_INFO_DEN2D
    )
    
    # Add padding between subplots
    fig.tight_layout()
    
    # Save figure
    filename = '_'.join(['den2d', run_mode, os.path.split(path)[-1]])
    filename = os.path.join(PLOT_DIR, filename + '.' + PLOT_FORMAT)
    fig.savefig(
        filename, dpi=DPI_DEN2D, format=PLOT_FORMAT, bbox_inches='tight'
    ) 
    plt.close(fig)
 

def plot_pairden(pairden_t, pairden_s, xfix, path, run_mode):
    # Initialize figure
    fig = plt.figure(figsize=plt.figaspect(0.5))
    
    # Initialize GridSpec
    gs = gridspec.GridSpec(1, 2)
    
    # Initialize ax and surf list
    ax, surf = [], []
    ax.append(fig.add_subplot(gs[0, 0], adjustable='box', aspect=1.0))
    ax.append(fig.add_subplot(gs[0, 1], adjustable='box', aspect=1.0)) 
    surf.append(None)
    surf.append(None)   

    # _ax_plot arguments
    args = [
        [ax[0], pairden_t, 'jet', '(ut + dt) / 2'],
        [ax[1], pairden_s, 'seismic', '[ud + du - (uu + dd)] / 2'],
    ]
    
    # Plot surfaces and mark fixed point
    for i, _args in enumerate(args):
        ax[i], surf[i] = _ax_plot(*_args)
        
        if _args[-1] == '(ut + dt) / 2':
            mark_color = 'white'
        if _args[-1] == '[ud + du - (uu + dd)] / 2':
            mark_color = 'black'
            
        ax[i].scatter(
            xfix[0], xfix[1], 
            marker='x', color=mark_color, 
            s=70, linewidth=2
        )
    
    # Colorbars and fonts
    for i, _ax in enumerate(ax):
        divider = make_axes_locatable(_ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(surf[i], cax=cax, ticks=LinearLocator(5))
        cbar.ax.tick_params(labelsize=LABELSIZE_TICK_PARAMS_PAIRDEN) 
        if DRAW_TICKS:
            _ax.set_xlabel('x', fontsize=FONTSIZE_XYLABEL_PAIRDEN)
            _ax.set_ylabel('y', fontsize=FONTSIZE_XYLABEL_PAIRDEN)
            _ax.tick_params(labelsize=LABELSIZE_TICK_PARAMS_PAIRDEN)
        else:
            _ax.set_xticks([])
            _ax.set_yticks([])
        if DRAW_TITLE:
            _ax.set_title(args[i][-1], fontsize=FONTSIZE_TITLE_PAIRDEN)
        
    # Information ax
    #ax_info = fig.add_subplot(1,1,1, adjustable='box', aspect=1.0)
    #ax_info.clear()
    #ax_info.axis("off")
    info_str = os.path.split(path)[-1].replace('_', '\_') + '\n' \
               + '[%s], xfix = $(%.5f, %.5f), %.2f^{\\circ}$' \
               %(run_mode, xfix[0], xfix[1], xfix[2])
    fig.text(
        0.5, 0.915, info_str, ha='center', fontsize=FONTSIZE_INFO_PAIRDEN
    )
    
    # Add padding between subplots
    fig.tight_layout()
    
    # Save figure
    filename = '_'.join(['pairden', run_mode, os.path.split(path)[-1]])
    filename = os.path.join(PLOT_DIR, filename + '.' + PLOT_FORMAT)
    fig.savefig(
        filename, dpi=DPI_PAIRDEN, format=PLOT_FORMAT, bbox_inches='tight'
    ) 
    plt.close(fig) 
    
    
def plot_pairden_all(
    pairden_dd, pairden_dt, pairden_du, 
    pairden_ud, pairden_ut, pairden_uu, 
    xfix, path, run_mode
):
    # Initialize figure
    fig = plt.figure(figsize=plt.figaspect(0.5))
    
    # Initialize GridSpec
    gs = gridspec.GridSpec(2, 3)
    
    # Initialize ax and surf list
    ax = []
    ax.append(fig.add_subplot(gs[0, 0], adjustable='box', aspect=1.0))
    ax.append(fig.add_subplot(gs[0, 1], adjustable='box', aspect=1.0)) 
    ax.append(fig.add_subplot(gs[0, 2], adjustable='box', aspect=1.0))
    ax.append(fig.add_subplot(gs[1, 0], adjustable='box', aspect=1.0))
    ax.append(fig.add_subplot(gs[1, 1], adjustable='box', aspect=1.0))
    ax.append(fig.add_subplot(gs[1, 2], adjustable='box', aspect=1.0))
    surf = [None for i in range(len(ax))]
    
    # _ax_plot arguments
    args = [
        [ax[0], pairden_ut, 'jet', 'ut'],
        [ax[1], pairden_ud, 'jet', 'ud'],
        [ax[2], pairden_uu, 'jet', 'uu'],
        [ax[3], pairden_dt, 'jet', 'dt'],
        [ax[4], pairden_du, 'jet', 'du'],
        [ax[5], pairden_dd, 'jet', 'dd'],
    ]
    
    # Plot surfaces
    for i, _args in enumerate(args):
        ax[i], surf[i] = _ax_plot(*_args)
 
    # Colorbars and fonts
    for i, _ax in enumerate(ax):
        divider = make_axes_locatable(_ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(surf[i], cax=cax, ticks=LinearLocator(5))
        cbar.ax.tick_params(labelsize=LABELSIZE_TICK_PARAMS_PAIRDEN_ALL)
        _ax.set_xticks([])
        _ax.set_yticks([])
        _ax.xaxis.set_ticklabels([])
        _ax.yaxis.set_ticklabels([])
        if DRAW_TITLE:
            _ax.set_title(args[i][-1], fontsize=FONTSIZE_TITLE_PAIRDEN_ALL)

    # Information ax
    #ax_info = fig.add_subplot(1,1,1, adjustable='box', aspect=1.0)
    #ax_info.clear()
    #ax_info.axis("off")
    info_str = os.path.split(path)[-1].replace('_', '\_') + '\n' \
               + '[%s], xfix = $(%.5f, %.5f), %.2f^{\\circ}$' \
               %(run_mode, xfix[0], xfix[1], xfix[2])
    fig.text(
        0.5, 0.990, info_str, ha='center', fontsize=FONTSIZE_INFO_PAIRDEN_ALL
    )
    
    # Add padding between subplots
    fig.tight_layout()
    
    # Save figure
    filename = '_'.join(['pairden_all', run_mode, os.path.split(path)[-1]])
    filename = os.path.join(PAIRDEN_ALL_DIR, filename + '.' + PLOT_FORMAT)
    fig.savefig(
        filename, dpi=DPI_PAIRDEN_ALL, format=PLOT_FORMAT, bbox_inches='tight'
    ) 
    plt.close(fig)
    
    
def _ax_plot(ax, data, cmap, title):    
    xx = data[:,:,0]
    yy = data[:,:,1]
    zz = data[:,:,2]
    
    xlim = [xx.min().min(), xx.max().max()]
    ylim = [yy.min().min(), yy.max().max()]
    zlim = [zz.min().min(), zz.max().max()]
    
    if title == 'Spin Density' or title == '[ud + du - (uu + dd)] / 2':
        zmax = abs(zz).max().max()
        zlim = [-zmax, zmax]

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)  
    ax.xaxis.set_major_locator(LinearLocator(5))
    ax.yaxis.set_major_locator(LinearLocator(5))  

    surf = ax.pcolormesh(
        xx, yy, zz, vmin=zlim[0], vmax=zlim[1],
        cmap=cmap, shading=SHADING
    )    
    
    return ax, surf


def get_chunks(n_data, n_cpu):
    chunksize, extra = divmod(n_data, n_cpu)
    if extra:
       chunksize += 1
    chunks = []
    for i in range(0, n_data, chunksize):
        chunks.append([i, i + chunksize])
    chunks[-1][-1] = n_data
    return chunks


if __name__ == '__main__': 
    plot_density()


