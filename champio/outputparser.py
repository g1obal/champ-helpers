"""
CHAMP Output Parser

Author: Gokhan Oztarhan
Created date: 27/01/2022
Last modified: 24/12/2023
"""

import os
from copy import deepcopy

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from .auconverter import AUConverter
from .densityparser import parse_den, parse_pairden


class OutputParser():
    FALLBACK = {
        bool: None,
        str: None,
        int: np.nan,
        float: np.nan,
    }

    def __init__(
        self, 
        fname='output_file', 
        path='.',
        m_r=None,
        kappa=None, 
        a=None,
        parse_density=True, 
        calculate_ss_corr=True,
        calculate_edge_pol=True,
        calculate_U_onsite=True,
    ):
        self.parent_path = None
        self.run_dir = None
        self.run_mode = None
        self.mpi_mode = None
        self.flk_type = None
        self.orb_dot_coef = np.nan
        self.info = None
        self.m_r = np.nan
        self.kappa = np.nan
        self.a = np.nan
        self.nelec = np.nan
        self.nup = np.nan
        self.ndn = np.nan
        self.norb = np.nan
        self.nbasis = np.nan
        self.gndot_v0 = np.nan
        self.gndot_rho = np.nan
        self.gndot_s = np.nan
        self.gndot_k = np.nan
        self.gauss_sigma = np.nan
        self.scalek = np.nan
        self.nopt_iter = np.nan
        self.add_diag = np.nan
        self.p_var = np.nan
        self.opt_type = None
        self.etrial = np.nan
        self.eshift = np.nan
        self.gauss_width_max = np.nan
        self.gauss_sigma_best = np.nan
        self.ss_corr_den = np.nan
        self.ss_corr_pairden = np.nan
        self.edge_pol_den = np.nan
        self.edge_pol_pairden = np.nan
        self.U_onsite_den = np.nan
        self.U_onsite_pairden = np.nan
        self.tot_E_1st = np.nan
        self.tot_E = np.nan 
        self.tot_E_err = np.nan
        self.pop_err = np.nan
        self.stdev_1st = np.nan
        self.stdev = np.nan
        self.stdev_err = np.nan
        self.ratio_int_kin = np.nan
        self.pot_E = np.nan
        self.pot_E_err = np.nan
        self.int_E = np.nan
        self.int_E_err = np.nan
        self.jf_kin_E = np.nan
        self.jf_kin_E_err = np.nan
        self.pb_kin_E = np.nan
        self.pb_kin_E_err = np.nan
        self.acceptance_1st = np.nan
        self.acceptance = np.nan 
        self.delta = np.nan
        self.Tcorr_1st = np.nan
        self.Tcorr = np.nan
        self.tau = np.nan
        self.overlap_dmc = np.nan
        self.nwalk_eff = np.nan
        self.nwalk_eff_f = np.nan
        self.nwalk_eff_fs = np.nan
        self.wts = np.nan
        self.wts_f = np.nan
        self.wts_fs = np.nan
        self.n_walkers_last = np.nan 
        self.n_walkers_global = np.nan
        self.n_walkers = np.nan
        self.n_cpu = np.nan
        self.n_block_eq = np.nan
        self.n_block = np.nan
        self.n_steps = np.nan
        self.n_block_max = np.nan
        self.tot_cpu_time = None

        # Variables declared until this line
        self.features = list(self.__dict__.keys())
        
        # Path
        self.path = path
        self.fname = fname
        
        # m_r, kappa and a
        if m_r is not None:
            self.m_r = m_r
        if kappa is not None:
            self.kappa = kappa
        if a is not None:
            self.a = a
        
        # Positions of the electrons or centers
        self.ncent = np.nan
        self.pos = np.nan
        
        # Best step of VMC optimization runs
        self.ind_best_run = np.nan
        self.ind_best_wf = np.nan
        
        # Indices of nearest neighbors, edges and bulk (interior points)
        self.ind_NN = np.nan
        self.ind_edge = np.nan
        self.ind_bulk = np.nan
        
        # Density variables
        self.den2d_parse_err = None
        self.pairden_parse_err = None
        self.den2d_t, self.den2d_s, self.den2d_nelec_calc, \
        self.pairden_dd, self.pairden_dt, self.pairden_du, \
        self.pairden_ud, self.pairden_ut, self.pairden_uu, \
        self.pairden_t, self.pairden_s, self.pairden_xfix = \
            [np.nan for i in range(12)]
        
        # Boolean for density parser
        self.parse_density = parse_density
        
        # Boolean for calculation of spin-spin correlations
        self.calculate_ss_corr = calculate_ss_corr
        
        # Boolean for calculation of edge polarization
        self.calculate_edge_pol = calculate_edge_pol
        
        # Boolean for calculation of U onsite (Hubbard U)
        self.calculate_U_onsite = calculate_U_onsite
        
    def dataframe(self):
        data = {}
        for feature in self.features:
            data[feature] = self.__dict__[feature]
        return pd.DataFrame(data, index=[0])
        
    def parse(self):
        path_split = os.path.split(self.path)
        self.parent_path = os.path.join(*path_split[:-1])
        self.run_dir = path_split[-1]
        
        try:
            with open(os.path.join(self.path, self.fname), 'r') as f:
                self.data = f.readlines()
                if self._data_corrupted():
                    self.data = None
                    print('[OutputParser] File corrupted: %s' \
                        %os.path.join(self.path, self.fname))
        except OSError:
            self.data = None
            print('[OutputParser] File read error: %s' \
                %os.path.join(self.path, self.fname))
        
        if self.data is not None:
            self._input_info()
            
            if self.run_mode == 'vmc':
                self._output_info_vmc()
                
            elif self.run_mode == 'dmc':
                self._split_data()
                self._output_info_dmc()
                
            if self.parse_density:
                self._den2d()
                self._pairden()
                self._indices_NN_edge_bulk()
                
            if self.calculate_ss_corr:
                self.ss_corr_den = self._ss_corr(self.den2d_t, self.den2d_s)
                self.ss_corr_pairden = \
                    self._ss_corr(self.pairden_t, self.pairden_s)
            
            if self.calculate_edge_pol:
                self.edge_pol_den = self._edge_pol(self.den2d_s)
                self.edge_pol_pairden = self._edge_pol(self.pairden_s)
            
            if self.calculate_U_onsite:
                self.U_onsite_den = self._U_onsite(self.den2d_t)
                self.U_onsite_pairden = self._U_onsite(self.pairden_t)
                
            self._cpu_time()
            
    def _data_corrupted(self):
        corrupted = False
        # Check the run start indicator line.
        string = '*********** START VMC CALCULATION  ***********'
        ind, line = get_lines(string, self.data, first_occurance=True)
        if ind is None:
            string = '*********** START DMC CALCULATION  ***********'
            ind, line = get_lines(string, self.data, first_occurance=True)
            if ind is None:
                corrupted = True
        return corrupted
        
    def _input_info(self):
        _feature = self._feature
        data = self.data
            
        self._run_mode()
        
        self._title_info()
        if self.flk_type is None and np.isnan(self.orb_dot_coef) \
            and np.isnan(self.m_r) and np.isnan(self.kappa) \
            and np.isnan(self.a) and self.info is None:
            self._title_info_backward_compatible()

        string = 'no. of electrons (all,up,dn) ='
        self.nelec = _feature(string, data, -3, int)
        self.nup = _feature(string, data, -2, int)
        self.ndn = _feature(string, data, -1, int)
        
        self.ncent = _feature('nctype,ncent =', data, -1, int)
        self._positions()

        self.norb = _feature('no. of orbitals =', data, -1, int)
        self.nbasis = _feature('no. of basis states =', data, -1, int)
        
        self.gndot_v0 = _feature('gndot_v0=', data, -1, float)
        self.gndot_rho = _feature('gndot_rho=', data, -1, float)
        self.gndot_s = _feature('gndot_s=', data, -1, float)
        try:
            ind, line = get_lines('GNDOT_K=', data, first_occurance=True)
            self.gndot_k = float(line.replace(',', '').split('=')[-1])
        except (IndexError, KeyError, AttributeError, TypeError, ValueError):
            self.gndot_k = np.nan
 
        gauss_sigma = _feature(
            'Floating gaussian widths:', data, 0, float, 
            relative_line_index=1
        )
        self.gauss_sigma = 1 / np.sqrt(gauss_sigma)
        
        self.scalek = _feature('scalek(1),a21=', data, -2, float)
        
        self.n_cpu = _feature('processors', data, -2, int)

        string = 'nopt_iter,nblk_max,add_diag(1),p_var,tol_energy='
        self.nopt_iter = _feature(string, data, 1, int)
        self.add_diag = _feature(string, data, 3, float, replace=['D','E'])
        self.p_var = _feature(string, data, 4, float, replace=['D','E'])
        if self.nopt_iter > 0:
            self.n_block_max = _feature(string, data, 2, int)
            self.gauss_width_max = _feature('GAUSS_WIDTH_MAX=', data, 1, float)
            self._opt_type()
        
        self.n_block_eq = _feature('no. of blocks before eq. =', data, -1, int)
        self.n_block = _feature('no. of blocks after eq.=', data, -1, int)
        self.n_steps = _feature('no. of steps/block =', data, -1, int)
        
        self.delta = _feature('step size =', data, -1, float)
        
        self.etrial = _feature('etrial', data, -1, float)
        self.tau = _feature('nfprod,tau', data, -1, float)
        
        self.n_walkers = _feature('target walker', data, -1, int)
        if self.mpi_mode is not None:
            if self.mpi_mode == 'mpi1':
                self.n_walkers_global = self.n_walkers
            elif self.mpi_mode == 'mpi2':
                self.n_walkers_global = self.n_walkers * self.n_cpu
        
        if self.gndot_k != 0 and not np.isnan(self.gndot_k):
            self.estimate_eshift()
            
        if np.isnan(self.a):
            try:
                distances = np.sqrt(
                    (self.pos[0,0] - self.pos[1:,0])**2 \
                    + (self.pos[0,1] - self.pos[1:,1])**2
                )
                self.a = distances.min()
            except (IndexError, AttributeError, TypeError, ValueError):
                self.a = np.nan
                
        if not np.isnan(self.orb_dot_coef):
            string = 'orbital coefficients - rows: norb, columns: nbasis'
            ind, line = get_lines(string, data, first_occurance=True)
            if self.orb_dot_coef == 1:
                if ind is None and line is None:
                    if isinstance(self.info, str):
                        self.info = '[WRONG orb_dot_coef] ' + self.info
                    else:
                        self.info = '[WRONG orb_dot_coef]'
            elif self.orb_dot_coef == 0:
                if ind is not None and line is not None:
                    if isinstance(self.info, str):
                        self.info = '[WRONG orb_dot_coef] ' + self.info
                    else:
                        self.info = '[WRONG orb_dot_coef]'
        
    def _output_info_vmc(self):
        """
        Parses VMC output.
        Reads all optimization steps.
        """
        _feature = self._feature
        _feature_all = self._feature_all
        _1st_best = self._1st_best
        data = self.data
        
        # tot_E (here to obtain best opt. step)
        self.tot_E_all = _feature_all('total E =', data, 3, float)
        self.tot_E_err_all = _feature_all('total E =', data, 5, float)
        
        # stdev (here to obtain best opt. step)
        self.stdev_all = _feature_all('total E =', data, -4, float)
        
        # optimization related variables
        if self.nopt_iter > 0:
            # best optimization step
            self._opt_step_best()

            # gauss_sigma_best
            self._gauss_sigma_best()
        
        # tot_E (remainig features)
        self.tot_E_1st, self.tot_E = _1st_best(self.tot_E_all, float)
        _, self.tot_E_err = _1st_best(self.tot_E_err_all, float)
        
        # stdev (remainig features)
        self.stdev_err_all = _feature_all('total E =', data, -2, float)
        self.stdev_1st, self.stdev = _1st_best(self.stdev_all, float)
        _, self.stdev_err = _1st_best(self.stdev_err_all, float)
        
        # Tcorr
        self.Tcorr_all = _feature_all('total E =', data, -1, float)
        self.Tcorr_1st, self.Tcorr = _1st_best(self.Tcorr_all, float)
        
        # pot_E
        self.pot_E_all = _feature_all('potential E', data, 3, float)
        self.pot_E_err_all = _feature_all('potential E', data, 5, float)
        _, self.pot_E = _1st_best(self.pot_E_all, float)
        _, self.pot_E_err = _1st_best(self.pot_E_err_all, float)
        
        # int_E
        self.int_E_all = _feature_all('interaction E', data, 3, float)
        self.int_E_err_all = _feature_all('interaction E', data, 5, float)
        _, self.int_E = _1st_best(self.int_E_all, float)
        _, self.int_E_err = _1st_best(self.int_E_err_all, float)
        
        # jf_kin_E
        self.jf_kin_E_all = _feature_all('jf kinetic E', data, 4, float)
        self.jf_kin_E_err_all = _feature_all('jf kinetic E', data, 6, float)
        _, self.jf_kin_E = _1st_best(self.jf_kin_E_all, float)
        _, self.jf_kin_E_err = _1st_best(self.jf_kin_E_err_all, float)
        
        # pb_kin_E
        self.pb_kin_E_all = _feature_all('pb kinetic E', data, 4, float)
        self.pb_kin_E_err_all = _feature_all('pb kinetic E', data, 6, float)
        _, self.pb_kin_E = _1st_best(self.pb_kin_E_all, float)
        _, self.pb_kin_E_err = _1st_best(self.pb_kin_E_err_all, float)
        
        # acceptance
        self.acceptance_all = _feature_all(
            'acceptance', data, -1, float, startswith=True
        )
        self.acceptance_1st, self.acceptance = _1st_best(
            self.acceptance_all, float
        )
        
        # ratio_int_kin
        self._ratio_int_kin_vmc()
            
    def _output_info_dmc(self):
        """
        Parses DMC output.
        """
        _feature = self._feature
        _feature_all = self._feature_all
        _1st_best = self._1st_best
        data_vmc = self.data_vmc
        data_dmc = self.data_dmc
        
        # nwalk_eff
        self.nwalk_eff = _feature('nwalk_eff/nwalk', data_dmc, -2, float)
        self.nwalk_eff_f = \
            _feature('nwalk_eff/nwalk with f', data_dmc, -2, float)
        self.nwalk_eff_fs = \
            _feature('nwalk_eff/nwalk with fs', data_dmc, -2, float)
        
        # weights
        self.wts = _feature('weights =', data_dmc, 2, float)
        self.wts_f = _feature('wts with f =', data_dmc, 4, float)
        self.wts_fs = _feature('wts with fs =', data_dmc, 4, float)
    
        # tot_E vmc
        self.tot_E_1st = _feature('total E =', data_vmc, 3, float)
        
        # tot_E_0 dmc
        tot_E_0 = _feature('total energy (   0)', data_dmc, -6, float)
        
        # tot_E dmc
        self.tot_E = _feature('total energy (  50)', data_dmc, -6, float)
        self.tot_E_err = _feature('total energy (  50)', data_dmc, -4, float)      
        
        # pop_err
        self.pop_err = np.abs(self.tot_E - tot_E_0)
        
        # stdev vmc
        self.stdev_1st = _feature('total E =', data_vmc, -4, float)
        
        # stdev dmc
        self.stdev = _feature('total energy (  50)', data_dmc, -2, float)  
        
        # Tcorr dmc
        self.Tcorr = _feature('total energy (  50)', data_dmc, -1, float)  
        
        # pot_E
        self.pot_E = _feature('potential energy', data_dmc, 3, float)
        self.pot_E_err = _feature('potential energy', data_dmc, 5, float)
        
        # int_E
        self.int_E = _feature('interaction energy', data_dmc, 3, float)
        self.int_E_err = _feature('interaction energy', data_dmc, 5, float)
        
        # jf_kin_E
        self.jf_kin_E = _feature('jf kinetic energy', data_dmc, 4, float)
        self.jf_kin_E_err = _feature('jf kinetic energy', data_dmc, 6, float)
        
        # pb_kin_E
        self.pb_kin_E = _feature('pb kinetic energy', data_dmc, 4, float)
        self.pb_kin_E_err = _feature('pb kinetic energy', data_dmc, 6, float)
        
        # acceptance vmc
        # acceptance_1st is the acceptance of the vmc run at the start of dmc
        # run if mcconfigs is not provided to dmc. That's why, value_best
        # is assigned to acceptance_1st, below.
        acceptance_vmc = _feature_all(
            'acceptance', data_vmc, -1, float, startswith=True
        )
        _, self.acceptance_1st = _1st_best(acceptance_vmc, float)
        
        # acceptance dmc
        self.acceptance = _feature(
            'No/frac. of node crossings,acceptance=', data_dmc, -1, float
        )
        
        # ratio_int_kin
        self._ratio_int_kin_dmc()
        
        # n_walkers_last
        self.n_walkers_last = _feature(
            'No. of walkers at end of run=', data_dmc, -1, float
        )
        
        # overlap_dmc
        string = 'approx. normalized overlap of FN and trial wave functions='
        self.overlap_dmc = _feature(string, data_dmc, -1, float)
        
    #------------------------#
    # Class helper functions #
    #------------------------#
    
    def _run_mode(self):
        run_mode = self._feature('Run in mode:', self.data, -1, _type=str)
        
        try:
            self.run_mode = run_mode.split('_')[0] 
        except (IndexError, AttributeError):
            self.run_mode = None
            
        try:
            mpi_mode = run_mode.split('_')[-1] 
            if 'mpi' in mpi_mode:
                self.mpi_mode = mpi_mode
            else:
                self.mpi_mode = None
        except (IndexError, AttributeError):
            self.mpi_mode = None

    def _title_info(self):
        _feature = self._feature
        data = self.data

        ind, line = get_lines('TITLE:', data, first_occurance=True)
        
        if line is not None:
            lines = line.split(':')[-1].replace('\n', '').split(',')
            
            for string in lines:        
                if string.startswith('ft='):
                    self.flk_type = string.split('=')[-1]
                    
                elif string.startswith('m_r='):
                    if np.isnan(self.m_r):
                        self.m_r = float(string.split('=')[-1])
                    
                elif string.startswith('kappa='):
                    if np.isnan(self.kappa):
                        self.kappa = float(string.split('=')[-1])
            
                elif string.startswith('a='):
                    if np.isnan(self.a):
                        self.a = float(string.split('=')[-1])
                     
                elif string.startswith('odc='):
                    self.orb_dot_coef = int(string.split('=')[-1])
                    
                elif string.startswith('info='):
                    self.info = string.split('=')[-1]
                    if self.info == 'None':
                        self.info = None 
                        
    def _title_info_backward_compatible(self):
        _feature = self._feature
        data = self.data
        
        self.flk_type = _feature('rtitle', data, 2, str)
        
        orb_dot_coef = _feature('rtitle', data, 3, str)
        if orb_dot_coef is not None:
            if orb_dot_coef == 'gauss':
                self.orb_dot_coef = 0
            elif orb_dot_coef == 'tb' or orb_dot_coef == 'mfh':
                self.orb_dot_coef = 1
        
        ind, line = get_lines('rtitle', data, first_occurance=True)
        if line is not None:
            lines = line.split()
            
            for string in lines:
                if string.startswith('m_r='):
                    if np.isnan(self.m_r):
                        self.m_r = float(string.split('=')[-1])
                    
                elif string.startswith('kappa='):
                    if np.isnan(self.kappa):
                        self.kappa = float(string.split('=')[-1])
            
                elif string.startswith('a='):
                    if np.isnan(self.a):
                        self.a = float(string.split('=')[-1])
                    
                elif string.startswith('addinfo='):
                    self.info = string.split('=')[-1]
                    if self.info == 'None':
                        self.info = None 
            
    def _opt_type(self):
        _feature = self._feature
        _feature_all = self._feature_all
        data = self.data
        
        jastrow = False
        gwidth = False
        
        string = 'nparm,nparml,nparmj,nparmcsf,nparms,nparmg,nparme='
        nparmj = _feature(string, data, -5, int)
        if nparmj > 0:
            jastrow = True
            
        string = 'orbital parameters varied='
        nparmo = _feature_all(string, data, -1, int)[2]
        if nparmo > 0:
            gwidth = True
            
        if jastrow and gwidth:
            self.opt_type = 'both'
        elif jastrow:
            self.opt_type = 'jastrow'
        elif gwidth:
            self.opt_type = 'gwidth'
            
    def _positions(self):
        data = self.data
        self.pos = np.zeros((self.ncent,2))
        
        ind, line = get_lines('center positions', data, first_occurance=True)
        i = 0
        for ind_data in range(ind + 1, ind + 1 + self.ncent):
            line = data[ind_data].split()
            self.pos[i,0], self.pos[i,1] = line[-2], line[-1]
            i += 1
            
    def _ratio_int_kin_vmc(self):
        try:
            self.ratio_int_kin_all = self.int_E_all / self.jf_kin_E_all
            self.ratio_int_kin = self.ratio_int_kin_all[-1]
        except (IndexError, AttributeError, TypeError, ValueError):
            try:
                self.ratio_int_kin_all = self.int_E_all / self.pb_kin_E_all
                self.ratio_int_kin = self.ratio_int_kin_all[-1]
            except (IndexError, AttributeError, TypeError, ValueError):
                self.ratio_int_kin_all = np.nan
                self.ratio_int_kin = np.nan
                
    def _ratio_int_kin_dmc(self):
        try:
            self.ratio_int_kin = self.int_E / self.jf_kin_E
        except (IndexError, AttributeError, TypeError, ValueError):
            try:
                self.ratio_int_kin = self.int_E / self.pb_kin_E
            except (IndexError, AttributeError, TypeError, ValueError):
                self.ratio_int_kin = np.nan
    
    def _opt_step_best(self):
        try:
            # CHAMP determines the best wave function according to 
            # energy_plus_err (ending of opt_wf subroutine in opt_wf.f90).
            # If the run is terminated before best wave function output,
            # use the parameters of the optimization step in which 
            # energy+3*error+p_var*sigma is minimum. The index has minus 1 
            # since energy is calculated after the parameters are set.
            energy_err_sigma = self.tot_E_all + 3 * self.tot_E_err_all \
                + self.p_var * self.stdev_all
            # nanargmin ignores NaN values
            self.ind_best_run = np.nanargmin(energy_err_sigma)
            self.ind_best_wf = self.ind_best_run - 1
        except (IndexError, AttributeError, TypeError, ValueError):
            self.ind_best_run = np.nan
            self.ind_best_wf = np.nan
    
    def _gauss_sigma_best(self):
        _feature = self._feature
        _feature_all = self._feature_all
        data = self.data
        
        string = '(floating_gauss_x_width_best(it,i),i=1,nbasis)'
        gauss_sigma_best = _feature(string, data[-50:], 0, float)
        
        string = '(floating_gauss_x_width_new(it,i),i=1,nbasis)'
        gauss_sigma_all = _feature_all(string, data, 0, float)
        if self.ind_best_wf >= 0:
            gauss_sigma_best_alt = gauss_sigma_all[self.ind_best_wf]
        else:
            gauss_sigma_best_alt = np.nan
        
        if np.isnan(gauss_sigma_best):
            gauss_sigma_best = gauss_sigma_best_alt
            
        self.gauss_sigma_all = 1 / np.sqrt(gauss_sigma_all)
        self.gauss_sigma_best = 1 / np.sqrt(gauss_sigma_best)
        
    def _split_data(self):
        string = '*********** START DMC CALCULATION  ***********'
        ind, line = get_lines(string, self.data, first_occurance=True)
        self.data_vmc = self.data[:ind]
        self.data_dmc = self.data[ind:]
        
    def _cpu_time(self):
        _feature = self._feature
            
        cpu_time = _feature(
            'The program ended normally.', self.data[-50:], -1, str, 
            relative_line_index=-1
        )
        
        try:
            hours = int(cpu_time.split()[-1][:4])
            mins = int(cpu_time.split()[-1][5:7])
            secs = int(cpu_time.split()[-1][8:10])
            self.tot_cpu_time = '%dh %dm %ds' %(hours,mins,secs)   
        except (IndexError, AttributeError, TypeError, ValueError):
            self.tot_cpu_time = None
    
    def to_SI(self, eunit='meV', lunit='nm'):
        auconverter = AUConverter(self.m_r, self.kappa)
        
        length = [
            'delta', 'a', 'gndot_rho', 'gauss_sigma', 'gauss_sigma_best'
        ]
        energy = [
            'gndot_v0', 'etrial', 'eshift', 'pop_err',
            'U_onsite_den', 'U_onsite_pairden',
            'tot_E_all', 'tot_E_1st', 'tot_E', 'tot_E_err',
            'pot_E_all', 'pot_E', 'pot_E_err', 
            'int_E_all', 'int_E', 'int_E_err', 
            'jf_kin_E_all', 'jf_kin_E', 'jf_kin_E_err', 
            'pb_kin_E_all', 'pb_kin_E', 'pb_kin_E_err',
        ]  
        
        err = (IndexError, KeyError, AttributeError, TypeError, ValueError)
             
        for feature in length:
            try:
                self.__dict__[feature] = auconverter.length_to_SI(
                    self.__dict__[feature], lunit)
            except err:
                self.__dict__[feature] = np.nan

        for feature in energy:
            try:
                self.__dict__[feature] = auconverter.energy_to_SI(
                    self.__dict__[feature], eunit)
            except err:
                self.__dict__[feature] = np.nan
                
        # Convert gndot_k from 1/au to [ENERGY]/[LENGTH]^2
        try:
            one_unit_SI_e = auconverter.energy_to_SI(1.0, eunit)    
            one_unit_SI_l = auconverter.length_to_SI(1.0, lunit) 
            self.gndot_k *= one_unit_SI_e / one_unit_SI_l**2    
        except err:
            self.gndot_k = np.nan        
                
    def opt_gauss_width(self):
        """
        Get optimized parameters as string
        """
        if self.data is not None:
            string = '(floating_gauss_x_width_best(it,i),i=1,nbasis)'
            ind, line = get_lines(string, self.data[-50:])
            best_wf_found = True
            
            if not line:
                string = '(floating_gauss_x_width_new(it,i),i=1,nbasis)'
                ind, line = get_lines(string, self.data)
                best_wf_found = False

                if self.ind_best_wf >= 0:
                    line = [line[self.ind_best_wf]]
                else:
                    line = []

            if line:
                return line[-1], best_wf_found
            else:
                return None
        else:
            return None
       
    def opt_jastrow(self):
        """
        Get optimized parameters as string
        """
        if self.data is not None:
            string_a = '(a_best(iparmj),iparmj=1,nparma)'
            string_b = '(b_best(iparmj),iparmj=1,nparmb)'
            string_c = '(c_best(iparmj),iparmj=1,nparmc)'
            ind, line_a = get_lines(string_a, self.data[-50:])
            ind, line_b = get_lines(string_b, self.data[-50:])
            ind, line_c = get_lines(string_c, self.data[-50:])
            best_wf_found = True
            
            if not (line_a and line_b and line_c):
                string_a = '(a_new(iparmj),iparmj=1,nparma)'
                string_b = '(b_new(iparmj),iparmj=1,nparmb)'
                string_c = '(c_new(iparmj),iparmj=1,nparmc)'
                ind, line_a = get_lines(string_a, self.data)
                ind, line_b = get_lines(string_b, self.data)
                ind, line_c = get_lines(string_c, self.data)
                best_wf_found = False
                
                if self.ind_best_wf >= 0:
                    line_a = [line_a[self.ind_best_wf]]
                    line_b = [line_b[self.ind_best_wf]]
                    line_c = [line_c[self.ind_best_wf]]
                else:
                    line_a, line_b, line_c = [], [], []
            
            if line_a and line_b and line_c:
                return line_a[-1], line_b[-1], line_c[-1], best_wf_found
            else:
                return None     
        else:
            return None
            
    def estimate_eshift(self):      
        self.eshift = (
            self.gndot_k * (self.pos[:,0]**2 + self.pos[:,1]**2)
        ).sum()
        
    def _indices_NN_edge_bulk(self):
        try:
            # Nearest neighbor indices of the lattice sites.
            # Algorithm also gives the neighbor with 0 distance, 
            # thus the point itself is also included in the NN array.  
            # In example, if there are 3 nearest neighbors, NN array for
            # that point consist of 4 elements.
            eps = self.a * 0.001
            neigh = NearestNeighbors(n_neighbors=4).fit(self.pos)
            self.ind_NN = neigh.radius_neighbors(
                self.pos, radius=self.a + eps, return_distance=False
            )
        except (IndexError, AttributeError, TypeError, ValueError):
            self.ind_NN = np.nan
            
        try:
            # Nearest neighbor count of each lattice site
            n_NN = np.array([array.shape[0] for array in self.ind_NN])
            
            # Indices of edge sites
            self.ind_edge = np.where(n_NN < 4)[0]

            # Indices of bulk sites
            self.ind_bulk = np.where(n_NN >= 4)[0] 
        except (IndexError, AttributeError, TypeError, ValueError):
            self.ind_edge = np.nan
            self.ind_bulk = np.nan
        
    def _den2d(self):
        try:
            self.den2d_t, self.den2d_s, self.den2d_nelec_calc = \
                parse_den(self.path, self.run_mode)
        except Exception as err:
            self.den2d_parse_err = err

    def _pairden(self):
        try:
            self.pairden_dd, self.pairden_dt, self.pairden_du, \
            self.pairden_ud, self.pairden_ut, self.pairden_uu, \
            self.pairden_t, self.pairden_s, self.pairden_xfix = \
            parse_pairden(
                self.path, self.run_mode, self.nelec, self.nup, self.ndn
            )
        except Exception as err:
            self.pairden_parse_err = err

    def _ss_corr(self, den_t, den_s, radius=None):
        """
        Normalized real space spin-spin correlation function
        g = <s_i * s_j> / <n_i * n_j> 
        where i and j are nearest neighbor sites in lattice.
        
        g = -1 : antiferromagnetic
        g =  0 : metallic
        g =  1 : ferromagnetic (or all spins align in the same direction)
        
        Normalization works since if the system seems antiferromagnetic
        but the numerical values of spin density are much lower than the
        total electron density, it is actually metallic. So, comparison of
        total densities and spin values is a check method for this purpose.
        """
        try:
            # Radius of spin values considered around potential centers 
            if radius is None:
                radius = self.a / 2
            
            # Electron densities and spin densities around potential centers
            dens = np.zeros(self.pos.shape[0])
            spins = np.zeros(self.pos.shape[0])
            for i in range(self.pos.shape[0]):
                dist = np.sqrt(
                    (self.pos[i,0] - den_s[:,:,0])**2 \
                    + (self.pos[i,1] - den_s[:,:,1])**2
                )
                indices = dist < radius
                dens[i] = den_t[indices, 2].mean()
                spins[i] = den_s[indices, 2].mean()
            
            # Loop over nearest neighbor indices to calculate ss_corr
            ss_sum = 0
            dd_sum = 0
            for i in range(self.ind_NN.shape[0]):
                for j in self.ind_NN[i]:
                    if i != j:
                        ss_sum += spins[i] * spins[j]
                        dd_sum += dens[i] * dens[j]

            # There is no need to calculate the mean of both ss_sum and dd_sum
            # since their denominator is the same.
            g_ss_corr = ss_sum / dd_sum
            
        except (IndexError, AttributeError, TypeError, ValueError):
            g_ss_corr = np.nan
        
        return g_ss_corr
        
    def _edge_pol(self, den_s, radius=None, include_corners=False):
        """
        Edge Polarization
        p = (<|s_edge|> - <|s_bulk|>) / (<|s_edge|> + <|s_bulk|>)
        where bulk is the interior lattice points.
        
        p =  1 : all total spin are polarized at the edges
        p =  0 : no polarization at the edges
        p = -1 : all total spin are polarized at the bulk
        """
        try:
            # Radius of spin values considered around potential centers 
            if radius is None:
                radius = self.a / 2
            
            # Spin densities around potential centers
            # Notice that np.abs() is used here. Thus, resulting array is |s|.
            spins = np.zeros(self.pos.shape[0])
            for i in range(self.pos.shape[0]):
                dist = np.sqrt(
                    (self.pos[i,0] - den_s[:,:,0])**2 \
                    + (self.pos[i,1] - den_s[:,:,1])**2
                )
                indices = dist < radius
                spins[i] = np.abs(den_s[indices, 2]).mean()
            
            # Drop corner sites from edge indices for triangular zigzag flake
            if not include_corners:
                dist = np.sqrt((self.pos[self.ind_edge,:]**2).sum(axis=1))
                ind_edge = self.ind_edge[dist < dist.max() - self.a / 8]
            else:
                ind_edge = self.ind_edge

            # <|s_edge|>
            s_edge = spins[ind_edge].mean()

            # <|s_bulk|>
            s_bulk = spins[self.ind_bulk].mean()
            
            # Edge polarization
            p_edge_pol = (s_edge - s_bulk) / (s_edge + s_bulk)
        
        except (IndexError, AttributeError, TypeError, ValueError):
            p_edge_pol = np.nan
        
        return p_edge_pol
        
    def _U_onsite(self, den_t, radius=None):
        """
        Onsite Coulomb Interaction (Hubbard U)
        Estimation from electron density
        
        U = \int n(x,y) * V(x,y) * dx * dy
        V(x',y') = \int n(x,y) / sqrt((x-x')^2 + (y - y')^2) * dx * dy
        where n is the electron density.
        """
        try:
            # Radius of spin values considered around potential centers 
            if radius is None:
                radius = self.a / 2
            
            # dx and dy for integration
            dx = den_t[1,0,0] - den_t[0,0,0]
            dy = den_t[0,1,1] - den_t[0,0,1]
            
            # Loop over potential centers
            U_onsite = 0
            for i in range(self.pos.shape[0]):
                dist = np.sqrt(
                    (self.pos[i,0] - den_t[:,:,0])**2 \
                    + (self.pos[i,1] - den_t[:,:,1])**2
                )
                indices = dist < radius
                x = den_t[indices,0]
                y = den_t[indices,1]
                n = den_t[indices,2]
                
                # Meshgrid for all combinations of points
                ind = np.arange(x.shape[0])
                ind1, ind2 = np.meshgrid(ind, ind)
                
                # Select non-diagonal elements
                mask = ~np.eye(ind1.shape[0], dtype=bool)
                ind1 = ind1[mask]
                ind2 = ind2[mask]
                
                # U_onsite for a single site omitting dx, dy
                r = np.sqrt((x[ind1] - x[ind2])**2 + (y[ind1] - y[ind2])**2)
                U_onsite += (n[ind1] * n[ind2] / r).sum()

            # Average value of U_onsite over all sites
            U_onsite *= dx * dy * dx * dy / self.pos.shape[0]
            
        except (IndexError, AttributeError, TypeError, ValueError):
            U_onsite = np.nan
        
        return U_onsite
    
    #-----------------------------------#
    # Class utilities:                  #
    # _1st_best, _feature, _feature_all #        
    #-----------------------------------#  
    
    def _1st_best(self, value_all, _type):
        try:
            if self.nopt_iter == 0:
                value_1st = np.nan
                value_best = value_all[-1]
            else:
                value_1st = value_all[0]
                value_best = value_all[self.ind_best_run]
        except (IndexError, AttributeError, TypeError, ValueError):
            value_1st = self.FALLBACK[_type]
            value_best = self.FALLBACK[_type]
        
        return value_1st, value_best
   
    def _feature(
        self, string, data, _index, _type, 
        relative_line_index=0, replace=None
    ):
        try:
            ind, line = get_lines(string, data, first_occurance=True)
            line = data[ind + relative_line_index]
            
            feature = line.split()[_index]
            
            if replace is not None:
                feature = feature.replace(*replace)
                
            feature = _type(feature)
        except (IndexError, AttributeError, TypeError, ValueError):
            feature = self.FALLBACK[_type]
        
        return feature
                      
    def _feature_all(
        self, string, data, _index, _type, ndarray=True, startswith=False
    ):
        try:
            ind, line = get_lines(string, data)
            
            if startswith:
                line = [s for s in line if s.startswith(string)]
            
            feature = []
            for _line in line:
                try:
                    value = _line.split()[_index]
                    feature.append(_type(value))
                except (IndexError, AttributeError, ValueError):
                    feature.append(self.FALLBACK[_type])
                    
            if ndarray == True:
                feature = np.array(feature)
        except (IndexError, AttributeError, TypeError, ValueError):
            feature = self.FALLBACK[_type]
        
        return feature

     
def get_lines(string, data, first_occurance=False):
    iterator = iter([i, s] for [i, s] in enumerate(data) if string in s)

    if first_occurance:
        return next(iterator, [None, None]) # index, line
    else:
        ind_line = [item for item in iterator]
        return [
            [item[0] for item in ind_line], # index
            [item[1] for item in ind_line], # line
        ]


