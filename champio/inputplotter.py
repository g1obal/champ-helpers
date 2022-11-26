"""
PLotting module for InputGenerator

Author: Gokhan Oztarhan
Created date: 09/06/2019
Last modified: 16/10/2022
"""

import time
from copy import deepcopy
import logging

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.ticker import (
    FixedLocator, LinearLocator, MultipleLocator, FormatStrFormatter
)
import matplotlib.gridspec as gridspec
import networkx

from .potential import gndot


logger = logging.getLogger(__name__)
logging.getLogger('numexpr').setLevel(logging.WARNING) # mute np.numexpr info

plt.ioff()
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
    'text.latex.preamble': r'\usepackage{amsmath}'
})


def get_1D_pot(inputgenerator, res=100):
    tic = time.time() 
    
    xlim_1D = \
        [-inputgenerator.gndot_rho_au * 2.5, inputgenerator.gndot_rho_au * 2.5]
    
    x_1D = np.linspace(xlim_1D[0], xlim_1D[1], res)
    y_1D = gndot(x_1D, inputgenerator.gndot_v0_au, \
        inputgenerator.gndot_rho_au, inputgenerator.gndot_s)
     
    ylim_1D = [min([0, y_1D.min()]), max([0, y_1D.max()])]   
    
    toc = time.time()   
    string = 'get_1D_pot done. (%.3f s)\n' %(toc-tic)
    logger.info(string)
            
    return xlim_1D, ylim_1D, x_1D, y_1D
            

def get_3D_pot_single_point(inputgenerator, res=300):
    tic = time.time() 
    
    meshlim_3D = deepcopy(inputgenerator.gndot_rho_au * 2.5)
    xy_3D = np.linspace(-meshlim_3D, meshlim_3D, res)

    xx_3D, yy_3D = np.meshgrid(xy_3D, xy_3D) 
    dist = np.sqrt(xx_3D * xx_3D + yy_3D * yy_3D)
    zz_3D = gndot(dist, inputgenerator.gndot_v0_au, \
        inputgenerator.gndot_rho_au, inputgenerator.gndot_s)
    
    xlim_3D = [-meshlim_3D, meshlim_3D]
    ylim_3D = xlim_3D     
    zlim_3D = [min([0, zz_3D.min().min()]), max([0, zz_3D.max().max()])]

    toc = time.time()   
    string = 'get_3D_pot_single_point done. (%.3f s)\n' %(toc-tic)
    logger.info(string)

    return xlim_3D, ylim_3D, zlim_3D, xx_3D, yy_3D, zz_3D


def get_3D_pot_lattice(inputgenerator, res=300):
    tic = time.time()    
    
    # Get potential radius
    radius = inputgenerator.gndot_rho_au * 1.65
    
    # Calculate limits
    xmin = inputgenerator.pos[:,0].min() - radius
    ymin = inputgenerator.pos[:,1].min() - radius
    xmax = inputgenerator.pos[:,0].max() + radius
    ymax = inputgenerator.pos[:,1].max() + radius
    meshlim = max([np.abs(xmin), np.abs(ymin), ymax, xmax])
    
    # Calculate potential mesh points
    xy = np.linspace(-meshlim, meshlim, res)
    xx, yy = np.meshgrid(xy, xy)  
    zz = np.zeros([xx.shape[0], yy.shape[0]])
    
    # gndot potential  
    for i in range(inputgenerator.pos.shape[0]):
        dist = np.sqrt(
            (xx - inputgenerator.pos[i,0])**2 \
            + (yy - inputgenerator.pos[i,1])**2
        ) 
        zz += gndot(dist, inputgenerator.gndot_v0_au, \
            inputgenerator.gndot_rho_au, inputgenerator.gndot_s)
        
    # Add quadratic gate potential if defined
    if inputgenerator.gndot_k != 0:
        dist2 = xx**2 + yy**2
        zz += inputgenerator.gndot_k_au * dist2

    xlim = [xmin, xmax]
    ylim = [ymin, ymax]
    zlim = [min([0, zz.min().min()]), max([0, zz.max().max()])]        
        
    toc = time.time()    
    string = 'get_3D_pot_lattice done. (%.3f s)\n' %(toc-tic)
    logger.info(string)

    return xlim, ylim, zlim, xx, yy, zz

    
def lattice_network(pos, ind_NN):  
    """
    Lattice network
    Connect the lattice points with lines
    """
    lattice_net = np.zeros([pos.shape[0], pos.shape[0]])
    lattice_net[ind_NN[:,0],ind_NN[:,1]] = 1
    lattice_net[ind_NN[:,1],ind_NN[:,0]] = 1   
    lattice_net = networkx.Graph(lattice_net)
    for i in range(pos.shape[0]):
        lattice_net.add_node(i,pos=(pos[i,0], pos[i,1]))
        
    node_positions = networkx.get_node_attributes(lattice_net, 'pos')   
    
    return lattice_net, node_positions
     
#-------------------------------------------------------------------------------

def inputplotter(
    inputgenerator, plot_fname='model.png', plot_dpi=300, plot_interactive=0
):
    logger.info('\n')    
    
    # Colormap
    cmap = 'jet'
    shading = 'nearest' # shading = 'auto', 'nearest', 'gouraud'
    
    # Calculate 1D potential
    xlim_1D, ylim_1D, x_1D, y_1D = get_1D_pot(inputgenerator)

    # Calculate 3D potential for single lattice point
    [xlim_3D, ylim_3D, zlim_3D, 
     xx_3D, yy_3D, zz_3D] = get_3D_pot_single_point(inputgenerator)  
    
    # Calculate 3D potentials on lattice
    xlim, ylim, zlim, xx, yy, zz = get_3D_pot_lattice(inputgenerator)
    
    # Calculate lattice network (connect the lattice points with lines)
    [lattice_net, node_positions] = \
        lattice_network(inputgenerator.pos, inputgenerator.ind_NN)

    # Begin plot
    tic = time.time() 
    fig = plt.figure(figsize=plt.figaspect(0.50))

    # Initialize GridSpec
    gs = gridspec.GridSpec(6, 8) 

    # Initialize ax list
    ax = []
    
    # If an axis is required to add manually, instead of using gridspec,
    # uncomment the following. The figure is at the left bottom for left and
    # bottom are (0,0). This procedure gives user more control on the position
    # and dimensions of axes.
    #left, bottom, width, height = [0.25, 0.25, 0.7, 0.7]
    #ax.append(
    #    fig.add_axes(
    #        [left, bottom, width, height], projection='3d', zorder=100
    #    )
    #)

    # ax[0]: 1D potential
    ax.append(fig.add_subplot(gs[:2, 0:2], box_aspect=0.7))
    ax[-1].set_xlim(*xlim_1D)
    ax[-1].set_ylim(*ylim_1D)
    ax[-1].tick_params(labelsize=6) 
    ax[-1].xaxis.set_major_locator(
        FixedLocator([xlim_1D[0], -inputgenerator.gndot_rho_au, \
            inputgenerator.gndot_rho_au, xlim_1D[1]])
    )
    ax[-1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax[-1].yaxis.set_major_locator(
        FixedLocator([ylim_1D[0], inputgenerator.gndot_v0_au/2, ylim_1D[1]])
    )
    ax[-1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))   
    ax[-1].plot(x_1D,y_1D)
    
    # ax[1]: 3D projection of potential for single lattice point
    ax.append(fig.add_subplot(gs[2:4, 0:2], projection='3d'))
    ax[-1].patch.set_facecolor('none') # background is transparent
    ax[-1].patch.set_alpha(0.0) # background is transparent
    ax[-1].set_xlim(*xlim_3D)
    ax[-1].set_ylim(*ylim_3D)
    ax[-1].set_zlim(*zlim_3D)
    ax[-1].set_xticks([])
    ax[-1].set_yticks([])
    ax[-1].set_zticks([])
    # Set viewpoint.
    ax[-1].azim = 135
    ax[-1].elev = 27
    # The surface is made opaque by using antialiased=False.
    ax[-1].plot_surface(
        xx_3D, yy_3D, zz_3D, cmap=cmap, linewidth=0, antialiased=True
    )
    
    # ax[2]: 2D mesh of potential for single lattice point
    ax.append(fig.add_subplot(gs[4:, 0:2], adjustable='box', aspect=1.0))
    ax[-1].set_xlim(*xlim_3D)
    ax[-1].set_ylim(*ylim_3D)
    ax[-1].set_xlabel('')
    ax[-1].set_ylabel('')
    ax[-1].set_yticks([])
    ax[-1].set_xticks([])
    ax[-1].pcolormesh(xx_3D, yy_3D, zz_3D, cmap=cmap, shading=shading) 

    # ax[3]: 3D projection of potential on lattice 
    ax.append(fig.add_subplot(gs[0:3, 2:4], projection='3d'))
    ax[-1].patch.set_facecolor('none') # background is transparent
    ax[-1].patch.set_alpha(0.0) # background is transparent
    ax[-1].set_box_aspect((1, 1, 0.5))
    #ax[-1].axis('off')
    ax[-1].set_xlim(*xlim)
    ax[-1].set_ylim(*ylim)
    ax[-1].set_zlim(*zlim)
    ax[-1].azim = 45
    ax[-1].elev = 30
    #ax[-1].set_xlabel()
    #ax[-1].set_ylabel()
    #ax[-1].set_zlabel()
    ax[-1].set_zticks([])
    ax[-1].set_yticks([])
    ax[-1].set_xticks([])
    ax[-1].plot_surface(xx, yy, zz, cmap=cmap, linewidth=0, antialiased=True)
    
    # ax[4]: lattice network
    ax.append(fig.add_subplot(gs[3:, 2:4], adjustable='box', aspect=1.0))
    ax[-1].set_xlim(*xlim)
    ax[-1].set_ylim(*ylim)
    ax[-1].set_xlabel('x')
    ax[-1].set_ylabel('y')
    ax[-1].set_yticks([])
    ax[-1].set_xticks([])
    latticecolor = "black" # (#616161:gray) 
    latticealpha = 1.0 # 0.5 for printing lattice network on potentials
    latticewidth = 0.5
    networkx.draw_networkx(
        lattice_net, pos=node_positions, node_size=0, with_labels=False, 
        width=latticewidth, edge_color=latticecolor, alpha=latticealpha,
        ax=ax[-1]
    )
   
    # ax[5]: 2D mesh of potential on lattice
    ax.append(fig.add_subplot(gs[:, 4:], adjustable='box', aspect=1.0))
    ax[-1].set_xlim(*xlim)
    ax[-1].set_ylim(*ylim)
    ax[-1].set_xlabel('x')
    ax[-1].set_ylabel('y')
    ax[-1].xaxis.set_major_locator(LinearLocator(5))   
    ax[-1].yaxis.set_major_locator(LinearLocator(5))
    ax[-1].tick_params(labelsize=8) 
    surf = ax[-1].pcolormesh(
        xx, yy, zz, vmin=zlim[0], vmax=zlim[1], cmap=cmap,shading=shading
    )
    
    # Add colorbar
    # create an ax on the right side of ax[-1]. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax[-1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(surf, cax=cax, ticks=LinearLocator(5))
    cbar.ax.tick_params(labelsize=8) 
    #cbar.ax.tick_params(width=0.5, length=2.25, labelsize=8) 
    #cbar.outline.set_linewidth(0.5)
    #cbar.ax.yaxis.get_offset_text().set(size=5)  

    # Change ax border widths and tick params excluding mesh figure
    for _ax in ax[:5]:
        for axis in ['top','bottom','left','right']:
            _ax.spines[axis].set_linewidth(0.5)
            _ax.tick_params(width=0.5, length=2.25)  
    
    # Information ax
    ax_info = fig.add_subplot(1,1,1)
    ax_info.clear()
    ax_info.axis('off')
    info_str = '$V(\\textbf{R}) = V_{0} \\sum\\limits_{\\textbf{R}_{0}}' \
        + '\\exp[-(|\\textbf{R}-\\textbf{R}_{0}|^{2} / \\rho^{2})^{s}]' \
        + '+ k |\\textbf{R}|^{2}$'
    ax_info.text(0.5, 1.150, info_str, ha='center', fontsize=8)
    v0_SI = ' $(%.3f \\textrm{ %s})$' \
        %(inputgenerator.gndot_v0, inputgenerator.eunit)
    rho_SI = ' $(%.3f \\textrm{ %s})$' \
        %(inputgenerator.gndot_rho, inputgenerator.lunit)
    k_SI = ' $(\\textrm{%.5e} \\textrm{ %s} \\cdot \\textrm{%s}^{-2})$' \
        %(inputgenerator.gndot_k, inputgenerator.eunit, inputgenerator.lunit)
    info_str = '$V_{0}=%.5f$' %(inputgenerator.gndot_v0_au) + v0_SI + ',    ' \
        + '$\\rho=%.5f$' %(inputgenerator.gndot_rho_au) + rho_SI + ',    ' \
        + '$s=%.3f$' %(inputgenerator.gndot_s) + ',    ' \
        + '$k=\\textrm{%.5e}$' %(inputgenerator.gndot_k_au) + k_SI             
    ax_info.text(0.5, 1.080, info_str, ha='center', fontsize=8)
    a_SI = ' $(%.3f \\textrm{ %s})$' %(inputgenerator.a, inputgenerator.lunit)
    info_str = \
        '$a(\\textrm{dot-to-dot distance})=%.5f$' %(inputgenerator.a_au) \
        + a_SI \
        + ',    ' \
        +  '$\\textrm{ncent}=%.d$' %(inputgenerator.ncent)
    ax_info.text(0.5, 1.035, info_str, ha='center', fontsize=8)  
    
    # Spaces between subplots
    # tight_layout automatically adjusts subplot sizes and spaces between them
    # in a way that subplots do not overlap. If you want to adjust them even
    # more, use subplots_adjust.
    plt.tight_layout()
    #plt.subplots_adjust(wspace=0.01, hspace=0.5)
    
    # Save figure
    fig.savefig(
        plot_fname, 
        dpi=plot_dpi, 
        format='png', 
        bbox_inches='tight'
    ) 
    plt.close(fig)
    
    # Interactive 3D lattice plot
    if plot_interactive:
        fig_interactive = plt.figure()
        ax_interactive = fig_interactive.add_subplot(1, 1, 1, projection='3d')
        ax_interactive.set_xlim(*xlim)
        ax_interactive.set_ylim(*ylim)
        ax_interactive.set_zlim(*zlim)
        ax_interactive.azim = 45
        ax_interactive.elev = 50
        ax_interactive.set_xlabel('')
        ax_interactive.set_ylabel('')
        ax_interactive.set_yticks([])
        ax_interactive.set_xticks([])
        ax_interactive.zaxis.set_major_locator(LinearLocator(0))
        ax_interactive.plot_surface(
            xx, yy, zz, cmap=cmap, linewidth=0, antialiased=True
        )
        plt.show()
        plt.close(fig_interactive)
    
    toc = time.time()   
    string = 'inputplotter done. (%.3f s)\n' %(toc-tic)
    logger.info(string)
    
