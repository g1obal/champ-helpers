"""
Lattice generator

Returns Lattice object.

Author: Gokhan Oztarhan
Created date: 17/11/2021
Last modified: 11/01/2023
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
        width=2, bc='', com_to_origin=False, **kwargs
    ):
        self.lat_type = lat_type # Not used for now
        self.flk_type = flk_type
        self.a = a
        self.n_side = n_side
        self.width = width
        self.bc = bc
        self.com_to_origin = com_to_origin
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
        
        self.n_site = self.pos.shape[0]
        
    def set(
        self, total_charge=None, Sz=None, 
        spin_order='AFM', spin_order_direction=1
    ):
        if total_charge is None:
            n_elec = self.n_site
        else:
            # there is a minus sign since electron charge is negative
            n_elec = self.n_site - total_charge
        
        # Check n_elec ranges. There should be at least 1 electron.
        # Total number of electrons should be less than or equal to 2 * n_site 
        # since there can be maximum 2 electrons (1 up, 1 dn) at a single site.
        if n_elec <= 0:
            sys.exit('\nn_elec should be greater than 0.')
        elif n_elec > 2 * self.n_site:
            sys.exit('\nn_elec should be less than or equal to 2 * n_site.')
            
        # Select the sublattice which has more sites
        ind_more = deepcopy(self.ind_A)
        ind_less = deepcopy(self.ind_B)
        if ind_less.shape[0] > ind_more.shape[0]:
            ind_more, ind_less = ind_less, ind_more
        
        # Set n_up and n_dn according to Lieb's theorem
        if Sz is None:
            # n_up and n_dn for a charge neutral system, n_elec = n_site
            n_up_cn = ind_more.shape[0]
            n_dn_cn = ind_less.shape[0]
            
            # Calculate the difference of n_up and n_dn
            # for a charge neutral system
            n_updn_diff = n_up_cn - n_dn_cn
            
            # Distribute the electrons to sublattice sites
            if n_elec <= self.n_site - n_updn_diff:
                n_up = int(n_elec // 2 + n_elec % 2)
                n_dn = n_elec - n_up
            elif n_elec <= self.n_site:
                n_add = n_updn_diff - (self.n_site - n_elec)
                n_up = n_dn_cn + n_add
                n_dn = n_dn_cn
            elif n_elec <= self.n_site + n_updn_diff:
                n_add = n_elec - self.n_site
                n_up = n_up_cn
                n_dn = n_dn_cn + n_add
            elif n_elec <= 2 * self.n_site:
                n_up = int(n_elec // 2 + n_elec % 2)
                n_dn = n_elec - n_up
            
        # Set n_up and n_dn according to Sz        
        else:
            n_up = int(n_elec // 2 + n_elec % 2 + Sz)
            n_dn = n_elec - n_up
            # Check Sz range
            if n_up > n_elec or n_up > self.n_site \
                or n_dn > n_elec or n_dn > self.n_site or n_dn > n_up:
                if n_elec >= self.n_site:
                    n_up_max = self.n_site
                    n_dn_max = n_elec - n_up_max
                else:
                    n_up_max = n_elec
                    n_dn_max = 0
                n_up_min = int(n_elec // 2 + n_elec % 2)
                n_dn_min = n_elec - n_up_min
                Sz_max = (n_up_max - n_dn_max) * 0.5
                Sz_min = (n_up_min - n_dn_min) * 0.5
                sys.exit(
                    '\nSz should be in range [%.1f, %.1f] for this system.' \
                    %(Sz_min, Sz_max)
                )

        # Set up and down indices
        if spin_order == 'AFM':
            # Sort indices according to their distances from center
            ind_draft_up = ind_more[
                np.argsort(np.sqrt((self.pos[ind_more,:]**2).sum(axis=1)))
            ]
            ind_draft_dn = ind_less[
                np.argsort(np.sqrt((self.pos[ind_less,:]**2).sum(axis=1)))
            ]
        
        elif spin_order == 'FM':
            # Sort indices according to their distances from center,
            # split them from y = - a / 8, in order to set up indices
            # to the sites in the upper half, and down indices to lower half.
            ind_sorted = np.argsort(np.sqrt((self.pos**2).sum(axis=1)))
            eps = self.a * 0.125
            ind_draft_up = ind_sorted[self.pos[ind_sorted,1] >= 0e0 - eps]
            ind_draft_dn = ind_sorted[self.pos[ind_sorted,1] < 0e0 - eps]
                
        # Add electrons from inside to outside, however add additional electrons
        # from outside to inside. For example, if there are additional up
        # electrons (Sz can be greater than 0 or total charge can be different 
        # than zero), these up electrons are located to outside lattice sites 
        # firstly, and continue to be located in outside to inside direction.
        if spin_order_direction:
            if n_up <= ind_draft_up.shape[0]:
                ind_up = ind_draft_up[:n_up]
            else:
                remaining = n_up - ind_draft_up.shape[0]
                ind_up = np.append(
                    ind_draft_up, ind_draft_dn[::-1][:remaining]
                )
            if n_dn <= ind_draft_dn.shape[0]:
                ind_dn = ind_draft_dn[:n_dn]
            else:
                remaining = n_dn - ind_draft_dn.shape[0]
                ind_dn = np.append(
                    ind_draft_dn, ind_draft_up[::-1][:remaining]
                )
                
        # Add electrons from outside to inside, however add additional electrons
        # from inside to outside.
        else:
            if n_up <= ind_draft_up.shape[0]:
                ind_up = ind_draft_up[::-1][:n_up]
            else:
                remaining = n_up - ind_draft_up.shape[0]
                ind_up = np.append(
                    ind_draft_up, ind_draft_dn[:remaining]
                )
            if n_dn <= ind_draft_dn.shape[0]:
                ind_dn = ind_draft_dn[::-1][:n_dn]
            else:
                remaining = n_dn - ind_draft_dn.shape[0]
                ind_dn = np.append(
                    ind_draft_dn, ind_draft_up[:remaining]
                )
        
        # Sort indices
        ind_up = np.sort(ind_up)
        ind_dn = np.sort(ind_dn)
        
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

    def _debug_plot_spins(self):
        s = self.a * 0.05
        
        fig = plt.figure(figsize=plt.figaspect(1.0))
        
        ax = []
        ax.append(fig.add_subplot(1, 1, 1))
        
        ax[-1].scatter(self.pos[self.ind_up,0], self.pos[self.ind_up,1], c='r')        
        ax[-1].scatter(
            self.pos[self.ind_dn,0], self.pos[self.ind_dn,1], marker='x', c='b'
        )
        
        for i in range(self.pos.shape[0]):
            string = '%d' %(i) 
            ax[-1].text(self.pos[i,0] + s, self.pos[i,1] + s, string)
        
        plt.show()
        sys.exit()


