"""
Input Generator class

Author: Gokhan Oztarhan
Created date: 18/06/2019
Last modified: 24/12/2023
"""

from os import urandom
import time
import logging

import numpy as np

from .auconverter import AUConverter
from .latgen import Lattice
from .inputwriter import inputwriter
from .inputplotter import inputplotter


logger = logging.getLogger(__name__)


class InputGenerator():
    _FLK_TYPE = {
        0: 'hexagonal_zigzag',
        1: 'hexagonal_armchair',
        2: 'triangular_zigzag',
        3: 'triangular_armchair',
        4: 'nanoribbon',
    }                 
    _FLK_TYPE_SHORT = {
        0: 'hexa_zz',
        1: 'hexa_armc',
        2: 'tri_zz',
        3: 'tri_armc',
        4: 'ribbon',
    }
    # np.polyval coefficients for estimating delta
    _DELTA_LOW = [4.14030228, 0.00847788]      
    _DELTA_MID = [
        -8.17462212e-04,  1.16552195e-02,
        -3.74729187e-02, -1.90940829e-01,
         1.60492769e+00, -4.27313546e+00,
         9.00066805e+00, -1.99794162e+00,
    ]
    _DELTA_HIGH = [2.89073666, 4.2907438]   

    def __init__(self):
        # [file]
        self.champ_input_file_name = 'input_file'
        
        # [title]
        self.info = None # optional additional info
        
        # [au]
        self.m_r = 0.067 # Electron mass ratio, m_eff / m_e
        self.kappa = 12.4  # Dielectric constant, varepsilon / varepsilon_0
        
        # [units]
        self.eunit = 'meV'
        self.lunit = 'nm'
        
        # [potential]
        self.gndot_v0 = -25.28
        self.gndot_rho = 20
        self.gndot_s = 1.40 # unitless
        self.gndot_k = 0 # unit: [ENERGY] / [LENGTH]^2
        
        # [basis]
        self.gauss_sigma = 10.97 # width guess for Gaussian basis
        
        # [lattice]
        self.a = 50
        self.n_side = 4
        self.width = 4 # for nanoribbon
        self.bc = 'xy' # for nanoribbon
        self.lat_type = 'honeycomb'
        self.flk_type = 1 # 0: hexagonal_zigzag, 
                          # 1: hexagonal_armchair, 
                          # 2: triangular_zigzag, 
                          # 3: triangular_armchair, 
                          # 4: nanoribbon
        
        # [electron]
        self.total_charge = None # set total number of electrons, nelec
                                 # None or 0 for charge neutral system
        self.Sz = None # total spin;
                       # to calculate the number of up and down electrons
                       # None to arrange nup and ndn according to Lieb's theorem
        self.spin_order = 'AFM' # electrons are located in a spin order
                                # AFM: antiferromagnetic, FM: ferromagnetic
        self.spin_order_direction = 1 # the direction in which electrons are 
                                      # located
                                      # 0: add electrons from outside to inside, 
                                      #    additional electrons from inside to 
                                      #    outside.
                                      # 1: add electrons from inside to outside, 
                                      #    additional electrons from outside to 
                                      #    inside.
        
        # [orbital]
        # orb_dot_coef, >0: overwrites spin_order
        # 0: gauss, 1: orbitals read from orb_dot_coef file 
        self.orb_dot_coef = 0
        
        # [random seed champ]
        # irn: a string of 16 digit integer
        #   'auto': automatically set random seed using os.urandom
        self.irn = 'auto'
        
        # [run]
        self.nstep = 100 # number of steps per block
        self.nblk = 15
        self.nblkeq = 5
        
        # [dmc]
        self.etrial = 0.1
        self.n_walkers = 25
        self.tau = 0.1
        
        # [jastrow]
        self.scalek = 0.2 # 0.2 is default
        self.nctype_of_edges = 1 # 1: same iwtype for all dots
                                 # 2: dots at the edges have different iwtype
                                 #    (corners not included)
        
        # [optional champ]
        self.ifixe = 0 #  0: do not write 2d density, 
                       # -1: write out 2d density
                       # -3: write out full pair density and 2d density
        self.xfix_pos = None # fix electron position, ex: [0.1, 0.4]
                             # None to find automatically
        self.xfix_angle = 60 # fix electron symmetry, default is 60
        
        # [opt]
        self.opt_mode = 0 # 0: both, 1: only width, 2: only jastrow
        self.opt_constraint = 1
        self.nopt_iter = 30
        self.nblk_max = 100 # maximum number of blocks 
                            # which can be increased during the optimization
        self.add_diag = 1e-4 # CHAMP uses abs(add_diag)
                             # negative sign for fixed add_diag
                             # positive sign for optimization of add_diag
        self.p_var = 0.1 # 0: energy, 1:variance
        self.tol_energy = 1e-12 # energy tolerance to finish optimization
        self.iopt = '00002' # last digit 2 is newton, 
                            # 01002 also a good choice, 
                            # 31101 is linear (bad choice)
        self.ipr_opt = -2 # -2: minimal output, 2: maximum output

    def config(self, config_dict):
        self.__dict__.update(config_dict) 
        
    def write(self):
        inputwriter(self)
        
    def plot(self, plot_fname='model.png', plot_dpi=300, plot_interactive=0):
        inputplotter(
            self, plot_fname=plot_fname, plot_dpi=plot_dpi, 
            plot_interactive=plot_interactive
        )

    def set(self):     
        self.set_units()
        self.set_delta()
        self.set_lattice()
        self.set_nctype()
        self.set_opt()
        self.set_irn()
        
        # Gaussian widths are implemented as 1/(width**2) in champ
        self.gauss_width_2 = 1.0 / self.gauss_sigma_au**2 
        
        # champ finds the nearest mesh point for a given coordinate
        # find the nearest lattice point to origin in the positive region
        # 60 degrees rotation for hexagon symmetry
        if self.xfix_pos is None:
            maskx = self.pos[:,0] > 0
            masky = self.pos[:,1] > 0
            mask = np.logical_and(maskx, masky)
            r = np.sqrt(self.pos[mask,0]**2 + self.pos[mask,1]**2)
            ind = r.argmin()
            x = self.pos[mask,0][ind]
            y = self.pos[mask,1][ind]
        else:
            x, y = self.xfix_pos[0], self.xfix_pos[1]    
            
        self.xfix = [x, y, self.xfix_angle]
        
        # Set title
        self.set_title()
        
    def set_units(self):
        eunit = self.eunit
        lunit = self.lunit
        auconverter = AUConverter(self.m_r, self.kappa)
                 
        self.gndot_v0_au = auconverter.energy_to_au(self.gndot_v0, eunit)
        self.gndot_rho_au = auconverter.length_to_au(self.gndot_rho, lunit)
        
        one_unit_e = auconverter.energy_to_au(1.0, eunit)
        one_unit_l = auconverter.length_to_au(1.0, lunit)
        self.gndot_k_au = self.gndot_k * one_unit_e / one_unit_l**2
        
        self.gauss_sigma_au = auconverter.length_to_au(self.gauss_sigma, lunit)
                                             
        self.a_au = auconverter.length_to_au(self.a, lunit)
        
    def set_delta(self):
        # Calculate polynomial fitted delta for around 0.5 acceptance
        if self.gauss_sigma_au <= 1:
            self.delta = np.polyval(self._DELTA_LOW, self.gauss_sigma_au) 
        elif self.gauss_sigma_au <= 5.5:
            self.delta = np.polyval(self._DELTA_MID, self.gauss_sigma_au)
        else:
            self.delta = np.polyval(self._DELTA_HIGH, self.gauss_sigma_au)

    def set_lattice(self):
        # Shift lattice center of mass to origin for CHAMP
        if 'triangular' in self._FLK_TYPE[self.flk_type]:
            com_to_origin = True
        else:
            com_to_origin = False

        lattice = Lattice(
            self.lat_type, self._FLK_TYPE[self.flk_type],
            self.a_au, self.n_side, width=self.width, bc=self.bc,
            com_to_origin=com_to_origin
        )
        
        lattice.set(
            total_charge=self.total_charge,
            Sz=self.Sz,
            spin_order=self.spin_order,
            spin_order_direction=self.spin_order_direction
        )
        
        self.ncent = lattice.n_site
        self.nelec = lattice.n_elec
        
        self.pos = lattice.pos
        self.ind_NN = lattice.ind_NN
        
        self.nup = lattice.n_up
        self.ndn = lattice.n_dn
        
        # Fortran indices start from 1
        if self.orb_dot_coef:
            self.ind_up = np.arange(1, self.nup + 1)
            self.ind_dn = np.arange(self.nup + 1, self.nelec + 1)
        else:
            self.ind_up = lattice.ind_up + 1
            self.ind_dn = lattice.ind_dn + 1
        
        self.nbasis = self.ncent
        
        if self.orb_dot_coef:
            self.norb = self.nelec
        else:
            self.norb = self.nbasis
        
        self.znuc = float(self.nelec) / self.ncent
        
        self.latmax = np.abs(lattice.pos).max().max()
        self.xmax = self.latmax + self.gauss_sigma_au + self.a_au
        
    def set_nctype(self):
        self.nctype = 1
        self.iwtype = np.full(self.ncent, 1)
        
        self.ind_edge = np.nan
        
        if self.nctype_of_edges == 2:
            # Find edge indices
            unique, counts = \
                np.unique(self.ind_NN.flatten(), return_counts=True)
            ind_sorted = np.argsort(unique)
            counts = counts[ind_sorted]
            self.ind_edge = np.where(counts < 3)[0]
            
            # Drop corner sites from edge indices for triangular zigzag flake
            include_corners = False
            
            if not include_corners:
                dist = np.sqrt((self.pos[self.ind_edge,:]**2).sum(axis=1))
                self.ind_edge = self.ind_edge[dist < dist.max() - self.a_au / 8]
            
            # Set nctype and iwtype of edges
            self.nctype = 2
            self.iwtype[self.ind_edge] = 2
        
    def set_opt(self):
        if self.opt_constraint:
            self.nparmo_3 = -1
        else:
            self.nparmo_3 = self.nbasis
        
        self.nparml = 0
        self.nparma = 4
        self.nparmb = 6
        self.nparmc = 15
        self.nparmf = 0
        self.nparmcsf = 0
        self.nparms = 0
        self.nparmg = 0
        self.nparmo_1 = 0
        self.nparmo_2 = 0
            
        if self.opt_mode == 1:
            self.nparma = 0
            self.nparmb = 0
            self.nparmc = 0
        elif self.opt_mode == 2:
            self.nparmo_3 = 0
        
        nparm = self.nparml + self.nparmb + self.nparmcsf \
            + self.nparms + self.nparmg + self.nparmo_1 + self.nparmo_2 \
            + np.abs(self.nparmo_3)
        
        for i in range(self.nctype):
            nparm += self.nparma
            nparm += self.nparmc
            nparm += self.nparmf
            
        self.nparm = np.int(nparm)
  
    def set_title(self):
        self.title = "'" \
            + 'TITLE:' \
            + 'ft=%s,' %self._FLK_TYPE_SHORT[self.flk_type] \
            + 'm_r=%.5f,' %self.m_r \
            + 'kappa=%.2f,' %self.kappa \
            + 'a=%.5f,' %self.a_au \
            + 'odc=%s,' %self.orb_dot_coef \
            + 'info=%s' %str(self.info) \
            + "'"
        
    def set_irn(self):
        if self.irn == 'auto':
            np.random.seed(
                int.from_bytes(urandom(4), 'big', signed=False)
            )
            ints = np.random.randint(1, 10, 1)
            ints = np.append(ints, np.random.randint(0, 10, 15))
            self.irn = ''.join([str(i) for i in ints])[:17]
        

