"""
Writer module for InputGenerator

Author: Gokhan Oztarhan
Created date: 11/06/2019
Last modified: 14/01/2023
"""

import time
import logging

import numpy as np


logger = logging.getLogger(__name__)


def generate_line(settings_str, info_str):
    skiplen = 40 - len(settings_str)
    string = '%-s' %settings_str \
             + ''.join([' ' for i in range(max(1, skiplen))]) \
             + '%-s' %info_str \
             + '\n'
    return string


def inputwriter(inputgenerator):
    tic = time.time()
    
    # Initialize input string
    input_str = ''
    
    # General Settings 
    settings_str = inputgenerator.title
    info_str = 'title'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = inputgenerator.irn[:16]
    info_str = 'irn'
    input_str += generate_line(settings_str, info_str)
    
    s1 = '{:d}'.format(0) # iperiodic
    s2 = '{:d}'.format(4) # ibasis
    s3 = "'gaussian'"
    settings_str = ' '.join([s1,s2,s3])
    info_str = 'iperiodic,ibasis'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = "0.5 {:.5f} '  Hartrees'".format(inputgenerator.etrial)
    info_str = 'hb,etrial,eunit'
    input_str += generate_line(settings_str, info_str) 
     
    s1 = '{:d}'.format(inputgenerator.nstep) # nstep
    s2 = '{:d}'.format(inputgenerator.nblk) # nblk
    s3 = '{:d}'.format(inputgenerator.nblkeq) # nblkeq
    s4 = '{:d}'.format(inputgenerator.n_walkers) # nconf
    s5 = '{:d}'.format(0) # nconf_new
    settings_str = ' '.join([s1,s2,s3,s4,s5])
    info_str = 'nstep,nblk,nblkeq,nconf,nconf_new'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '0  0  1  -4'
    info_str = 'idump,irstar,isite,ipr'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '1  {:.2f}  4.  1.  1.'.format(inputgenerator.delta)
    info_str = 'imetro,delta,deltar,deltat,fbias'
    input_str += generate_line(settings_str, info_str)
     
    settings_str = '2 1 1 1 1 0 0 0 0'
    info_str = 'idmc,ipq,itau_eff,iacc_rej,icross,icuspg,idiv_v,icut_br,icut_e'
    input_str += generate_line(settings_str, info_str) 
    
    settings_str = '50  {:.3f}'.format(inputgenerator.tau)
    info_str = 'nfprod,tau'
    input_str += generate_line(settings_str, info_str)
    
    s1 = '{:d}'.format(-7) # nloc
    s2 = '{:d}'.format(0) # numr
    s3 = '{:d}'.format(1) # nforce
    s4 = '{:d}'.format(0) # nefp
    settings_str = '  '.join([s1,s2,s3,s4])  
    info_str = 'nloc,numr,nforce,nefp'
    input_str += generate_line(settings_str, info_str)
     
    s1 = '{:.5f}'.format(inputgenerator.gndot_v0_au) # gndot_v0
    s2 = '{:.5f}'.format(inputgenerator.gndot_rho_au) # gndot_rho
    s3 = '{:.5f}'.format(inputgenerator.gndot_s) # gndot_s
    settings_str = '  '.join([s1,s2,s3]) 
    info_str = 'gndot_v0,gndot_rho,gndot_s'
    input_str += generate_line(settings_str, info_str)
    
    s1 = '{:d}'.format(inputgenerator.nelec) # nelec
    s2 = '{:d}'.format(inputgenerator.nup) # nup
    settings_str = ' '.join([s1,s2]) 
    info_str = 'nelec,nup'
    input_str += generate_line(settings_str, info_str)
    
    input_str += '\n'
    
    # Geometry Section 
    input_str += "'* Geometry section'\n"
    
    s1 = '{:d}'.format(2) # ndim
    settings_str = s1
    info_str = 'ndim'
    input_str += generate_line(settings_str, info_str)
    
    s1 = '{:d}'.format(1) # nctype
    s2 = '{:d}'.format(inputgenerator.ncent) # ncent
    settings_str = ' '.join([s1,s2]) 
    info_str = 'nctype,ncent'
    input_str += generate_line(settings_str, info_str)
    
    s1 = ''
    for i in range(inputgenerator.ncent):
        s1 += '{:d} '.format(1)
    settings_str = s1
    info_str = '(iwtype(i),i=1,ncent)'
    input_str += generate_line(settings_str, info_str)
    
    s1 = '{:.8f} '.format(inputgenerator.znuc) 
    settings_str = s1
    info_str = '(znuc(i),i=1,nctype)'
    input_str += generate_line(settings_str, info_str)
    
    info_str = '(cent(k,ic),k=1,ndim)'
    for i in range(inputgenerator.ncent):
        settings_str = '{0: .12f} {1: .12f}'.format(
            inputgenerator.pos[i,0], inputgenerator.pos[i,1]
        )
        input_str += generate_line(settings_str, info_str)
    
    input_str += '\n'

    # Determinantal section
    input_str += "'* Determinantal section'\n"  
    
    s1 = '{:d}'.format(0) # inum_orb
    settings_str = s1
    info_str = 'inum_orb'
    input_str += generate_line(settings_str, info_str)
    
    s1 = '{:d}'.format(1) # ndet
    s2 = '{:d}'.format(inputgenerator.nbasis) # nbasis
    s3 = '{:d}'.format(inputgenerator.norb) # norb
    settings_str = ' '.join([s1,s2,s3]) 
    info_str = 'ndet,nbasis,norb'
    input_str += generate_line(settings_str, info_str)    
    
    s1 = ''
    for i in range(inputgenerator.nbasis):
        s1 += '{: .12f} '.format(inputgenerator.pos[i,0])
    settings_str = s1
    info_str = '(floating_gauss_x_pos(it,i),i=1,nbasis)'
    input_str += generate_line(settings_str, info_str)
    
    s1 = ''
    for i in range(inputgenerator.nbasis):
        s1 += '{: .12f} '.format(inputgenerator.pos[i,1])
    settings_str = s1
    info_str = '(floating_gauss_y_pos(it,i),i=1,nbasis)'
    input_str += generate_line(settings_str, info_str)
    
    s1 = ''
    for i in range(inputgenerator.nbasis):
        s1 += '{: .12f} '.format(inputgenerator.gauss_width_2)
    settings_str = s1
    info_str = '(floating_gauss_width(it,i),i=1,nbasis)' \
               + ' - these widths are "1/(width**2)"'
    input_str += generate_line(settings_str, info_str)
    
    s1 = ''
    for i in inputgenerator.ind_up:
        s1 += '{:d} '.format(i)
    s2 = ''
    for i in inputgenerator.ind_dn:
        s2 += '{:d} '.format(i)
    settings_str = ' '.join([s1,s2])
    info_str = '(iworbd(j,idet),j=1,nelec)'
    input_str += generate_line(settings_str, info_str)
    
    input_str += '\n'
    settings_str = '1'
    info_str = 'ncsf'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '1.0'
    info_str = '(csf_coef(icsf),icsf=1,ncsf)'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '1'
    info_str = '(ndet_in_csf(icsf),icsf=1,ncsf)'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '1'
    info_str = '(iwcsf(idet_in_csf,1),idet_in_csf=1,ndet_in_csf(1))'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '1.0'
    info_str = '(cdet_csf(idet_in_csf,1),idet_in_csf=1,ndet_in_csf(1))'
    input_str += generate_line(settings_str, info_str)
    
    input_str += '\n'
    
    # Jastrow Section
    input_str += "'* Jastrow section'\n"
    
    s1 = '{:d}'.format(1) # ianalyt_lap
    settings_str = s1
    info_str = 'ianalyt_lap'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '4 2 1 1 5 0'
    info_str = 'ijas,isc,nspin1,nspin2,nord,ifock'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '6 6 5'
    info_str = 'norda,nordb,nordc'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '{:.2f} 0.0'.format(inputgenerator.scalek)
    info_str = 'scalek,a21'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '0. 0. 0. 0. 0. 0. 0. '
    info_str = '(a(iparmj),iparmj=1,nparma)'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '1. 0.2 0. 0. 0. 0. 0. '
    info_str = '(b(iparmj),iparmj=1,nparmb)'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '0. 0. 0. 0. 0. 0. 0. 0. 0.  0. 0. 0. 0. 0. 0. 0.' \
                   + '  0. 0. 0. 0.  0. 0. 0. '
    info_str = '(c(iparmj),iparmj=1,nparmc)'
    input_str += generate_line(settings_str, info_str)
    
    input_str += '\n'
    
    # Optional Features
    input_str += "'* Optional features'\n"
    
    s1 = '&opt_list ifixe={:d}\n'.format(inputgenerator.ifixe) \
         + '          xfix=%.5f %.5f %.2f\n' \
            %(inputgenerator.xfix[0], 
              inputgenerator.xfix[1], 
              inputgenerator.xfix[2]) \
         + '          iperturb=0\n' \
         + '          icoosys=1\n' \
         + '          xmax={:.1f}d0\n'.format(inputgenerator.xmax) \
         + '          iantiferromagnetic=0\n' \
         + '          izigzag=0\n' \
         + '          gndot_k={:.5e}\n'.format(
            inputgenerator.gndot_k_au).replace('e', 'd') \
         + '/\n'
    input_str += s1

    input_str += '\n'

    # Optimization Section
    input_str += "'* Optimization section'\n"

    s1 = '{:3d}'.format(inputgenerator.nopt_iter) # nopt_iter
    s2 = '{:d}'.format(inputgenerator.nblk_max) # nblk_max
    s3 = '{:.0e}'.format(-5.e-1).replace('e', 'd') # add_diag(1)
    s4 = '{:.2f}'.format(inputgenerator.p_var) # p_var
    s5 = '{:.0e}'.format(1.e-8).replace('e', 'd') # tol_energy
    settings_str = ' '.join([s1,s2,s3,s4,s5])
    info_str = 'nopt_iter,nblk_max,add_diag(1),p_var,tol_energy'
    input_str += generate_line(settings_str, info_str)

    s1 = '{:d}'.format(2000) # NDATA
    s2 = '{:d}'.format(inputgenerator.nparm) # NPARM
    s3 = '{:d}'.format(-1) # icusp
    s4 = '{:d}'.format(1) # icusp2
    s5 = '{:d}'.format(2) # NSIG
    s6 = '{:d}'.format(1000) # NCALLS
    s7 = '{:s}'.format(inputgenerator.iopt) # iopt
    s8 = '{:d}'.format(inputgenerator.ipr_opt) # ipr
    settings_str = ' '.join([s1,s2,s3,s4,s5,s6,s7,s8])
    info_str = 'NDATA,NPARM,icusp,icusp2,NSIG,NCALLS,iopt,ipr_opt'
    input_str += generate_line(settings_str, info_str)

    settings_str = '0 0 0 0'
    info_str = 'i3body,irewgt,iaver,istrech'
    input_str += generate_line(settings_str, info_str)

    settings_str = '0 0 0 0 0 0 0 0 0 0'
    info_str = 'ipos,idcds,idcdr,idcdt,id2cds,id2cdr,id2cdt,idbds,idbdr,idbdt'
    input_str += generate_line(settings_str, info_str)
    
    s1 = ''
    for i in range(inputgenerator.norb):
        s1 += '{:d} '.format(-1) 
    settings_str = s1
    info_str = '(lo(iorb),iorb=1,norb)'
    input_str += generate_line(settings_str, info_str)
    
    s1 = '{:d}'.format(inputgenerator.nparml)
    s2 = '{:d}'.format(inputgenerator.nparma)
    s3 = '{:d}'.format(inputgenerator.nparmb)
    s4 = '{:d}'.format(inputgenerator.nparmc)
    s5 = '{:d}'.format(inputgenerator.nparmf)
    s6 = '{:d}'.format(inputgenerator.nparmcsf)
    s7 = '{:d}'.format(inputgenerator.nparms)
    s8 = '{:d}'.format(inputgenerator.nparmg)
    s9 = '{:d}'.format(inputgenerator.nparmo_1)
    s10 = '{:d}'.format(inputgenerator.nparmo_2)
    s11 = '{:d}'.format(inputgenerator.nparmo_3)
    s1 = ' '.join([s1,s2,s3,s4,s5,s6,s7,s8])
    s2 = ' '.join([s9,s10,s11])
    settings_str = '   '.join([s1,s2])
    info_str = 'nparml,nparma,nparmb,nparmc,nparmf,nparmcsf,nparms,nparmg,' \
               + '(nparmo(i),i=1,notype) - notype is related to ibasis'
    input_str += generate_line(settings_str, info_str)    
    
    settings_str = ' '
    info_str = '(iwo(iparm),iparm=1,nparmo(1))'
    input_str += generate_line(settings_str, info_str)     
    
    settings_str = ' '
    info_str = '(iwo(iparm),iparm=1,nparmo(2))'
    input_str += generate_line(settings_str, info_str)  

    s1 = ''
    for i in range(1,np.abs(inputgenerator.nparmo_3)+1):
        s1 += '{:d} '.format(i)
    settings_str = s1
    info_str = '(iwo(iparm),iparm=1,nparmo(3))'
    input_str += generate_line(settings_str, info_str)  
    
    s1 = ''
    for i in range(inputgenerator.nparml):
        s1 += '{:d} '.format(0)
    settings_str = s1
    info_str = '(iworb(iparm),iwbasi(iparm),iparm=1,nparml)'
    input_str += generate_line(settings_str, info_str)    
    
    s1 = ''
    for i in range(inputgenerator.nparm-inputgenerator.nparml):
        s1 += '{:d} '.format(0) 
    settings_str = s1
    info_str = '(iwbase(iparm),iparm=1,nparm-nparml)'
    input_str += generate_line(settings_str, info_str)    

    s1 = ''
    for i in range(inputgenerator.nparmcsf):
        s1 += '{:d} '.format(0)
    settings_str = s1
    info_str = '(iwcsf(iparm),iparm=1,nparmcsf)'
    input_str += generate_line(settings_str, info_str)   

    settings_str = '    3 4 5 6 7'
    info_str = '(iwjasa(iparm),iparm=1,nparma)'
    input_str += generate_line(settings_str, info_str)   
    
    settings_str = '  2 3 4 5 6 7'
    info_str = '(iwjasb(iparm),iparm=1,nparmb)'
    input_str += generate_line(settings_str, info_str)   
    
    settings_str = '    3   5   7 8 9    11    13 14 15 16 17 18' \
                   + '    20 21    23    '
    info_str = '(iwjasc(iparm),iparm=1,nparmc)'
    input_str += generate_line(settings_str, info_str)   

    settings_str = '0 0 0'
    info_str = 'necn,nebase - moved from fit'
    input_str += generate_line(settings_str, info_str) 
    
    settings_str = ' '
    info_str = '((ieorb(iorb,ibasis),iebasi(iorb,ibasis),iorb=1,2),' \
               + 'ibasis=1,necn) - moved from fit'
    input_str += generate_line(settings_str, info_str) 
    
    settings_str = ' '
    info_str = '((iebase(iorb,ibasis),iorb=1,2),ibasis=1,nebase)' \
               + ' - moved from fit'
    input_str += generate_line(settings_str, info_str) 
    
    s1 = ''
    for i in range(inputgenerator.norb):
        s1 += '{:d} '.format(0) 
    settings_str = s1
    info_str = '(ipivot(j),j=1,norb) (not needed if necn < 0) - moved from fit'
    input_str += generate_line(settings_str, info_str)
    
    settings_str = '0'
    info_str = 'eave (guess for energy, used in fit) - moved from fit'
    input_str += generate_line(settings_str, info_str) 
    
    settings_str = '1.d-6 5. 1 50 1'
    info_str = 'pmarquardt,tau,noutput,nstep,ibold (used in fit)' \
               + ' - moved from fit'
    input_str += generate_line(settings_str, info_str) 

    settings_str = 'T F'
    info_str = 'analytic,cholesky (used in fit) - moved from fit'
    input_str += generate_line(settings_str, info_str) 
    
    settings_str = '0 0 0'
    if inputgenerator.opt_constraint:
        settings_str = '0 0 {:d}'.format(inputgenerator.nbasis - 1)
    info_str = '(norb_constraints(i),i=1,notype)' \
               + ' - 0 means no constraint to orbitals'
    input_str += generate_line(settings_str, info_str) 
    
    settings_str = ' '
    info_str = '((orb_constraints(1,i,j),j=1,2),i=1,norb_constrants(1))'
    input_str += generate_line(settings_str, info_str)  

    settings_str = ' '
    info_str = '((orb_constraints(1,i,j),j=1,2),i=1,norb_constrants(2))'
    input_str += generate_line(settings_str, info_str)  

    settings_str = ''
    if inputgenerator.opt_constraint:
        for i in range(2,inputgenerator.nbasis + 1):
            settings_str += '{:d} 1 '.format(i) 
    else:
        settings_str += ' '
    info_str = '((orb_constraints(1,i,j),j=1,2),i=1,norb_constrants(3))'
    input_str += generate_line(settings_str, info_str)      

    input_str += 'end\n'
    input_str += 'exit\n'

    # Write to input_file and close
    with open(inputgenerator.champ_input_file_name, 'w') as input_file:
        input_file.write(input_str)

    toc = time.time()     
    string = '\nwrite_input done. (%.3f s)\n' %(toc-tic)
    logger.info(string)
    
