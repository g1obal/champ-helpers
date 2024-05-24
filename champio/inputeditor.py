"""
CHAMP Input Editor

Edits currently present inputs.
Updates optimized parameters using output file.

Author: Gokhan Oztarhan
Created date: 30/01/2022
Last modified: 24/05/2024
"""

import os

import numpy as np

from .outputparser import get_lines
from .inputwriter import generate_line


# Polynomial coefficients for fixing acceptance
POLY_COEF = np.array([
    -8.96914127e+05, 6.62592648e+06, -2.21200940e+07, 4.41408228e+07,
    -5.86875575e+07, 5.48595716e+07, -3.70830989e+07, 1.83644274e+07,
    -6.67791146e+06, 1.77055326e+06, -3.36530353e+05, 4.44860783e+04,
    -3.87221313e+03, 1.96201108e+02, -2.13410784e+00, 2.31992931e-01,
])


class InputEditor():
    def __init__(self, fname='input_file', path='.'):
        self.path = path
        
        try:
            with open(os.path.join(self.path, fname), 'r') as f:
                self.data = f.readlines()
        except OSError:
            self.data = None
            print(
                '[InputEditor] File read error: %s' \
                %os.path.join(self.path, fname)
            )
            
    def edit_irn(self, irn):
        if irn is not None:
            line = self.data[1].split()
            
            if not isinstance(irn, str):
                line[0] = '%d' %(irn)
            elif irn == 'auto':
                np.random.seed(
                    int.from_bytes(os.urandom(4), 'big', signed=False)  
                )
                ints = np.random.randint(1, 10, 1)
                ints = np.append(ints, np.random.randint(0, 10, 15))
                new_irn = ''.join([str(i) for i in ints])[:17]
                line[0] = '%s' %(new_irn)
                
            settings_str = ' '.join(line[:1])
            info_str = ' '.join(line[1:])
            self.data[1] = generate_line(settings_str, info_str)
               
    def edit_etrial(self, run_mode, tot_E, tot_E_err, pop_err, etrial):
        if etrial is not None:
            string = 'hb,etrial,eunit' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            energy_trial = np.nan
            
            if not isinstance(etrial, str):
                energy_trial = etrial
                
            elif etrial == 'auto' and run_mode is not None:
                if run_mode == 'vmc':
                    energy_trial = tot_E - tot_E_err
                    energy_trial = energy_trial - abs(energy_trial) * 0.005
                elif run_mode == 'dmc':
                    energy_trial = tot_E + max([tot_E_err, pop_err]) / 2
                
            if not np.isnan(energy_trial):    
                line[1] = '%.5f' %(energy_trial)            
                
                settings_str = ' '.join(line[:3]) + '  ' + line[3]
                info_str = ' '.join(line[4:])
                self.data[ind] = generate_line(settings_str, info_str)
                
                print(
                    '[InputEditor] edit_etrial: ' \
                    'tot_E = %.5f, etrial = %.5f' %(tot_E, float(line[1]))
                )
            
    def edit_nstep(self, nstep):
        if nstep is not None:
            string = 'nstep,nblk,nblkeq,nconf,nconf_new' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[0] = '%d' %(nstep)
            
            settings_str = ' '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str)
            
    def edit_nblk(self, nblk):
        if nblk is not None:
            string = 'nstep,nblk,nblkeq,nconf,nconf_new' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[1] = '%d' %(nblk)
            
            settings_str = ' '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str)       
            
    def edit_nblkeq(self, nblkeq):
        if nblkeq is not None:
            string = 'nstep,nblk,nblkeq,nconf,nconf_new' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[2] = '%d' %(nblkeq)
            
            settings_str = ' '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str) 
            
    def edit_nconf(self, nconf):
        if nconf is not None:
            string = 'nstep,nblk,nblkeq,nconf,nconf_new' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[3] = '%d' %(nconf)
            
            settings_str = ' '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str)  
            
    def edit_nconf_new(self, nconf_new):
        if nconf_new is not None:
            string = 'nstep,nblk,nblkeq,nconf,nconf_new' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[4] = '%d' %(nconf_new)
            
            settings_str = ' '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str)  
            
    def edit_isite(self, isite):
        if isite is not None:
            string = 'idump,irstar,isite,ipr' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[2] = '%d' %(isite)
            
            settings_str = ' '.join(line[:4])
            info_str = ' '.join(line[4:])
            self.data[ind] = generate_line(settings_str, info_str)
            
    def edit_delta(self, acceptance, delta):
        if delta is not None:
            string = 'imetro,delta,deltar,deltat,fbias' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
        
            if not isinstance(delta, str):
                line[1] = '%.2f' %(delta)
            elif delta == 'auto' and not np.isnan(acceptance):
                if acceptance < 0.05:
                    line[1] = '%.2f' %(float(line[1]) * 0.3) 
                elif acceptance > 0.95:
                    line[1] = '%.2f' %(float(line[1]) * 2.5) 
                elif acceptance < 0.49 or acceptance > 0.51:
                    ratio = np.polyval(POLY_COEF, acceptance)
                    line[1] = '%.2f' %(float(line[1]) * ratio)  
        
            settings_str = '  '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str)  
            
    def edit_tau(self, tau):
        if tau is not None:
            string = 'nfprod,tau' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[1] = '%.3f' %(tau)
            
            settings_str = '  '.join(line[:2])
            info_str = ' '.join(line[2:])
            self.data[ind] = generate_line(settings_str, info_str)  
            
    def edit_ifixe(self, ifixe):
        if ifixe is not None:
            string = 'ifixe=' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split('=')
    
            line[1] = '%d\n' %(ifixe)
            
            self.data[ind] = '='.join(line)
            
    def edit_xmax(self, xmax):
        if xmax is not None:
            string = 'xmax=' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split('=')
    
            line[1] = '%.2fd0\n' %(xmax)
            
            self.data[ind] = '='.join(line)
        
    def edit_nopt_iter(self, nopt_iter):
        if nopt_iter is not None:
            string = 'nopt_iter,nblk_max,add_diag(1),p_var,tol_energy' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[0] = '%3d' %(nopt_iter)
            
            settings_str = ' '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str)
            
    def edit_add_diag(self, add_diag):
        if add_diag is not None:
            string = 'nopt_iter,nblk_max,add_diag(1),p_var,tol_energy' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[2] = ('%.1e' %(add_diag)).replace('e', 'd')
            
            settings_str = ' '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str)
            
    def edit_p_var(self, p_var):
        if p_var is not None:
            string = 'nopt_iter,nblk_max,add_diag(1),p_var,tol_energy' 
            ind, line = get_lines(string, self.data, first_occurance=True)
            line = line.split()
            
            line[3] = '%.2f' %(p_var)
            
            settings_str = ' '.join(line[:5])
            info_str = ' '.join(line[5:])
            self.data[ind] = generate_line(settings_str, info_str)
            
    def edit_opt_mode(self, nbasis, constraint, opt_mode):
        if opt_mode is not None:
            if opt_mode == 0:
                self._nparm_jastrow(True)
                self._nparm_gauss_width(nbasis, constraint, True)
                self._nparm_pos(nbasis, True)
            elif opt_mode == 1:
                self._nparm_jastrow(False)
                self._nparm_gauss_width(nbasis, constraint, True)
                self._nparm_pos(nbasis, True)
            elif opt_mode == 2:
                self._nparm_jastrow(True)
                self._nparm_gauss_width(nbasis, constraint, True)
                self._nparm_pos(nbasis, False)
            elif opt_mode == 3:
                self._nparm_jastrow(True)
                self._nparm_gauss_width(nbasis, constraint, False)
                self._nparm_pos(nbasis, True)
            elif opt_mode == 4:
                self._nparm_jastrow(False)
                self._nparm_gauss_width(nbasis, constraint, False)
                self._nparm_pos(nbasis, True)
            elif opt_mode == 5:
                self._nparm_jastrow(False)
                self._nparm_gauss_width(nbasis, constraint, True)
                self._nparm_pos(nbasis, False)
            elif opt_mode == 6:
                self._nparm_jastrow(True)
                self._nparm_gauss_width(nbasis, constraint, False)
                self._nparm_pos(nbasis, False)
            
            self._nparm()
                
    def _nparm_jastrow(self, optimize):
        string = 'nparml,nparma,nparmb,nparmc,nparmf,nparmcsf,nparms,nparmg'
        ind, line = get_lines(string, self.data, first_occurance=True)
        line = line.split()
    
        if optimize:
            line[1] = '%i' %(4)
            line[2] = '%i' %(6)
            line[3] = '%i' %(15)
        else:
            line[1] = '%i' %(0)
            line[2] = '%i' %(0)
            line[3] = '%i' %(0)
        
        settings_str = ' '.join(line[:8]) + '   ' +  ' '.join(line[8:11])
        info_str = ' '.join(line[11:])
        self.data[ind] = generate_line(settings_str, info_str)  
        
    def _nparm_gauss_width(self, nbasis, constraint, optimize):
        string = 'nparml,nparma,nparmb,nparmc,nparmf,nparmcsf,nparms,nparmg'
        ind, line = get_lines(string, self.data, first_occurance=True)
        line = line.split()
    
        if optimize:
            if constraint:
                line[10] = '%i' %(-1)
            else:
                line[10] = '%i' %nbasis
        else:
            line[10] = '%i' %(0)
        
        settings_str = ' '.join(line[:8]) + '   ' +  ' '.join(line[8:11])
        info_str = ' '.join(line[11:])
        self.data[ind] = generate_line(settings_str, info_str)
        
    def _nparm_pos(self, nbasis, optimize):
        string = 'nparml,nparma,nparmb,nparmc,nparmf,nparmcsf,nparms,nparmg'
        ind, line = get_lines(string, self.data, first_occurance=True)
        line = line.split()
    
        if optimize:
            line[8] = '%i' %nbasis
            line[9] = '%i' %nbasis
        else:
            line[8] = '%i' %(0)
            line[9] = '%i' %(0)
        
        settings_str = ' '.join(line[:8]) + '   ' +  ' '.join(line[8:11])
        info_str = ' '.join(line[11:])
        self.data[ind] = generate_line(settings_str, info_str) 
        
    def _nparm(self):
        string = 'nparml,nparma,nparmb,nparmc,nparmf,nparmcsf,nparms,nparmg'
        ind, line = get_lines(string, self.data, first_occurance=True)
        line = line.split()
        
        nparm = sum([abs(int(item)) for item in line[:11]])
        
        string = 'NDATA,NPARM,icusp,icusp2,NSIG,NCALLS,iopt,ipr_opt'
        ind, line = get_lines(string, self.data, first_occurance=True)
        line = line.split()
        
        line[1] = '%i' %(nparm)
        
        settings_str = ' '.join(line[:8])
        info_str = ' '.join(line[8:])
        self.data[ind] = generate_line(settings_str, info_str)  

    def update_pos(self, new_params):
        if new_params is not None:
            string = '* Determinantal section'
            ind, line = get_lines(string, self.data, first_occurance=True)
            
            self.data[ind + 3] = new_params[0]
            self.data[ind + 4] = new_params[1]
            
            print('[InputEditor] update_pos: best_wf_found = %s' \
                %new_params[-1])
        else:
            print('[InputEditor] update_pos: new_params = None')

    def update_gauss_width(self, new_params):
        if new_params is not None:
            string = '* Determinantal section'
            ind, line = get_lines(string, self.data, first_occurance=True)
            
            self.data[ind + 5] = new_params[0]
            
            print('[InputEditor] update_gauss_width: best_wf_found = %s' \
                %new_params[-1])
        else:
            print('[InputEditor] update_gauss_width: new_params = None')
        
    def update_jastrow(self, new_params):
        if new_params is not None:
            string = '* Jastrow section'
            ind, line = get_lines(string, self.data, first_occurance=True)
            
            self.data[ind + 5] = new_params[0]
            self.data[ind + 6] = new_params[1]
            self.data[ind + 7] = new_params[2]
            
            print('[InputEditor] update_jastrow: best_wf_found = %s' \
                %new_params[-1])
        else:
            print('[InputEditor] update_jastrow: new_params = None')
            
    def generate_mc_configs(self, root, files, n_mc_configs, pick_random):
        mc_configs = []
        i = 0
        for f in files:
            if 'mc_configs_new' in f:
                with open(os.path.join(root, f), 'r') as f_temp:
                    mc_configs += f_temp.readlines()
                    i += 1
        
        if not mc_configs:
            f_callback = os.path.join(root, 'mc_configs_start')
            if os.path.exists(f_callback):
                with open(f_callback, 'r') as f_temp:
                    mc_configs += f_temp.readlines()
                    i += 1
        
        print(
            '[InputEditor] generate_mc_configs: ' \
            + 'mc_configs files found = %i\n' %i \
            + '[InputEditor] generate_mc_configs: ' \
            + 'total # of lines found = %i' %len(mc_configs)
        )
        
        if mc_configs:
            if pick_random:
                ind = np.random.permutation(np.arange(len(mc_configs)))
                mc_configs_shuffle = [mc_configs[j] for j in ind]
            else:
                mc_configs_shuffle = []
                iskip = len(mc_configs) // i
                for j in range(iskip):
                    mc_configs_shuffle += mc_configs[j::iskip]

            if n_mc_configs <= len(mc_configs_shuffle):
                mc_configs_shuffle = mc_configs_shuffle[:n_mc_configs]
                
            with open(os.path.join(root, 'mc_configs'), 'w') as f:
                f.writelines(mc_configs_shuffle)
            
            print('[InputEditor] generate_mc_configs: ' \
                + '# of lines in mc_configs = %i' %len(mc_configs_shuffle))

    def commit(self, fname):
        with open(os.path.join(self.path, fname), 'w') as f:
            f.writelines(self.data)
        
    
