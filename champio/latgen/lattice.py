"""
Lattice generator

Returns Lattice object.

Author: Gokhan Oztarhan
Created date: 17/11/2021
Last modified: 21/07/2022
"""

import time
import logging
import sys
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
import networkx

from .honeycomb import honeycomb


logger = logging.getLogger(__name__)


class Lattice():
    """
    Lattice
    
    Calculates the lattice points, positions, sub lattice indices,
    nearest neighbor indices, boundary conditions
    
    Parameters
    ----------
    lat_type : str
        lattice type
    flk_type : str
        flake type
    a : float, 
        dot-to-dot distance,
        lattice constant
    n_side : int
        number of points on one side of the flk_type
    width : int, optional
        multiple of 2
        vertical width of nanoribbon flk_type,
        vertical number of points
    bc : str, optional
        boundary condition for nanoribbon flk_type
    com_to_origin : bool, optional
        shift center of mass to origin
    debug : bool, optional
        debugging mode, plotting the lattice
        indicating the nearest neighbors    

    Returns
    -------
    Lattice object
    """                
    def __init__(
        self, lat_type, flk_type, a, n_side, 
        width=2, bc='', com_to_origin=False, debug=False, **kwargs
    ):
        self.lat_type = lat_type # Not used for now
        self.flk_type = flk_type
        self.a = a
        self.n_side = n_side
        self.width = width
        self.bc = bc
        self.com_to_origin = com_to_origin
        self.debug = debug
        self.kwargs = kwargs               
        self._honeycomb()
        
    def _honeycomb(self):    
        if self.width % 2 != 0: 
            self.width += 1
            
        self.pos, self.sub, self.m, self.n, \
        self.ind_A, self.ind_B, self.distances, \
        self.ind_NN, self.ind_NN_2nd, self.n_side = \
        honeycomb(
            self.flk_type, self.a, self.n_side, 
            width=self.width, bc=self.bc, 
            com_to_origin=self.com_to_origin,
            **self.kwargs
        )
                               
        self.n_total = self.pos.shape[0] 
        
        if self.debug:
            self._debug_plot()
        
    def set_spin(self, n_elec, spin_order='antiferromagnetic', Sz=None):
        ind_max = deepcopy(self.ind_A)
        ind_min = deepcopy(self.ind_B)
        if ind_min.shape[0] > ind_max.shape[0]:
            ind_max, ind_min = ind_min, ind_max
    
        if Sz is not None:
            n_up = int(n_elec // 2 + n_elec % 2 + Sz)
            n_dn = n_elec - n_up
        else:
            if spin_order == 'antiferromagnetic':
                n_up = ind_max.shape[0]
                n_dn = ind_min.shape[0]
            elif spin_order == 'ferromagnetic':
                n_up = n_elec
                n_dn = 0
    
        if spin_order == 'antiferromagnetic':
            ind_up = ind_max
            ind_dn = ind_min[:n_dn]
            if n_up != n_dn:
                ind_up = np.append(ind_up,ind_min[n_dn:])         
        elif spin_order == 'ferromagnetic':
            ind_up = np.arange(n_up)
            ind_dn = np.arange(n_up,n_elec)
            
        self.n_elec = n_elec
        self.n_up = n_up
        self.n_dn = n_dn
        self.ind_up = ind_up
        self.ind_dn = ind_dn  

    def _debug_plot(self):
        s = self.a * 0.05
        
        fig = plt.figure(figsize=plt.figaspect(1.0))
        
        ax = []
        ax.append(fig.add_subplot(1, 1, 1))
        
        ax[-1].scatter(self.pos[self.ind_A,0], self.pos[self.ind_A,1], c='r')
        ax[-1].scatter(self.pos[self.ind_B,0], self.pos[self.ind_B,1], c='b')
        
        for i in range(self.pos.shape[0]):
            string = '%d (%d %d)' %(i, self.m[i], self.n[i]) 
            ax[-1].text(self.pos[i,0] + s, self.pos[i,1] + s, string)
            
        lattice_net = np.zeros([self.pos.shape[0], self.pos.shape[0]])
        lattice_net[self.ind_NN[:,0],self.ind_NN[:,1]] = 1
        lattice_net[self.ind_NN[:,1],self.ind_NN[:,0]] = 1  
        lattice_net[self.ind_NN_2nd[:,0],self.ind_NN_2nd[:,1]] = 1
        lattice_net[self.ind_NN_2nd[:,1],self.ind_NN_2nd[:,0]] = 1
        lattice_net = networkx.Graph(lattice_net)
        for i in range(self.pos.shape[0]):
            lattice_net.add_node(i, pos=(self.pos[i,0], self.pos[i,1]))
        node_positions = networkx.get_node_attributes(lattice_net, 'pos')  
        
        networkx.draw_networkx(
            lattice_net, pos=node_positions, node_size=0, with_labels=False, 
            width=0.5, edge_color='#616161', alpha=0.5,
            ax=ax[-1]
        )
            
        plt.show()
        sys.exit()

        
