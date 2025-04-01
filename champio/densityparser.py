"""
CHAMP Density Parser

Author: Gokhan Oztarhan
Created date: 20/06/2022
Last modified: 26/03/2025
"""

import os

import numpy as np


def parse_den(path, run_mode):
    # Load density files
    fname_den2d_t = os.path.join(path, 'den2d_t_' + run_mode)
    fname_den2d_u = os.path.join(path, 'den2d_u_' + run_mode)
    fname_den2d_d = os.path.join(path, 'den2d_d_' + run_mode)
    
    # Total, up and down densities
    den2d_t = np.loadtxt(fname_den2d_t, dtype=float)
    den2d_u = np.loadtxt(fname_den2d_u, dtype=float)
    den2d_d = np.loadtxt(fname_den2d_d, dtype=float)
    
    # Parse and preprocess densities
    xdim = np.unique(den2d_t[:,0]).shape[0]
    ydim = np.unique(den2d_t[:,1]).shape[0]

    xx = np.reshape(den2d_t[:,0], (xdim,ydim))
    yy = np.reshape(den2d_t[:,1], (xdim,ydim))
    
    zz = np.reshape(den2d_t[:,2], (xdim,ydim)) 
    den2d_t = np.stack([xx, yy, zz], axis=2)
    
    up = np.reshape(den2d_u[:,2], (xdim,ydim))
    dn = np.reshape(den2d_d[:,2], (xdim,ydim))
    spin = up - dn
    den2d_s = np.stack([xx, yy, spin], axis=2)
    
    # Integral of electron density
    dx = xx[1,0] - xx[0,0]
    dy = yy[0,1] - yy[0,0]
    nelec_calc = np.trapz(
        np.trapz(den2d_t[:,:,2], dx=dx, axis=0), dx=dy, axis=0
    )
    
    return den2d_t, den2d_s, nelec_calc


def parse_pairden(path, run_mode, nelec, nup, ndn):
    # Load pair density files
    fname_pairden_dd = os.path.join(path, 'pairden_dd_' + run_mode)
    fname_pairden_dt = os.path.join(path, 'pairden_dt_' + run_mode)
    fname_pairden_du = os.path.join(path, 'pairden_du_' + run_mode)
    fname_pairden_ud = os.path.join(path, 'pairden_ud_' + run_mode)
    fname_pairden_ut = os.path.join(path, 'pairden_ut_' + run_mode)
    fname_pairden_uu = os.path.join(path, 'pairden_uu_' + run_mode)                    
    
    # Pair densities
    pairden_ud = np.loadtxt(fname_pairden_ud, skiprows=1, dtype=float)
    pairden_ut = np.loadtxt(fname_pairden_ut, skiprows=1, dtype=float)
    pairden_uu = np.loadtxt(fname_pairden_uu, skiprows=1, dtype=float)
    try:
        pairden_dd = np.loadtxt(fname_pairden_dd, skiprows=1, dtype=float)
        pairden_dt = np.loadtxt(fname_pairden_dt, skiprows=1, dtype=float)
        pairden_du = np.loadtxt(fname_pairden_du, skiprows=1, dtype=float)
    except (OSError, IndexError, AttributeError, TypeError, ValueError):
        pairden_dd = np.zeros(pairden_ut.shape)
        pairden_dt = np.zeros(pairden_ut.shape)
        pairden_du = np.zeros(pairden_ut.shape)         
    
    with open(fname_pairden_ut, 'r') as f:
        line = f.readline()
        xfix = [
            float(line.split()[-3]), 
            float(line.split()[-2]), 
            float(line.split()[-1])
        ]
        
    # Parse and preprocess pair densities
    xdim = np.unique(pairden_ut[:,0]).shape[0]
    ydim = np.unique(pairden_ut[:,1]).shape[0]

    xx = np.reshape(pairden_ut[:,0], (xdim,ydim))
    yy = np.reshape(pairden_ut[:,1], (xdim,ydim))
    
    zz = np.reshape(pairden_dd[:,2], (xdim,ydim)) 
    pairden_dd = np.stack([xx, yy, zz], axis=2)
    
    zz = np.reshape(pairden_dt[:,2], (xdim,ydim)) 
    pairden_dt = np.stack([xx, yy, zz], axis=2)
    
    zz = np.reshape(pairden_du[:,2], (xdim,ydim)) 
    pairden_du = np.stack([xx, yy, zz], axis=2)
    
    zz = np.reshape(pairden_ud[:,2], (xdim,ydim)) 
    pairden_ud = np.stack([xx, yy, zz], axis=2)
    
    zz = np.reshape(pairden_ut[:,2], (xdim,ydim)) 
    pairden_ut = np.stack([xx, yy, zz], axis=2)
    
    zz = np.reshape(pairden_uu[:,2], (xdim,ydim)) 
    pairden_uu = np.stack([xx, yy, zz], axis=2)

    # Integration variables
    dx = xx[1,0] - xx[0,0]
    dy = yy[0,1] - yy[0,0]

    # Normalization of dd, du, ud and uu should be done according to integral
    # of their sum, since (dt + ut) / 2 = (ud + du + uu + dd) / 2.
    # These are not assumptions, they are tested.
    # Their sum should give the total electron density (multiplied by 2).
    # !!! Their individual normalizations are not important since CHAMP does 
    # not keep track of their individual normalizations !!!
    # If up electrons are at the fixed positions more frequently, 
    # ud and uu should have larger weights than dd and du.
    # This normalization keeps track of that, as well.
    # For example; if only up electrons occupies the fixed point, 
    # all points in dt, du and dd will be close to zero (or equals to zero).
    # If the normalization of dt, du and dd are restored independently, 
    # the noise in those pair density calculations will be boosted.
    # Thus, the end results will be completely dominated by noise.
    # That's why, we need to normalize (ud + du + uu + dd) and (dt + ut).
    if nup == ndn:
        zz = pairden_ud[:,:,2] + pairden_du[:,:,2] \
            + pairden_uu[:,:,2] + pairden_dd[:,:,2]
        integral = np.trapz(np.trapz(zz, dx=dx, axis=0), dx=dy, axis=0)
        multiplier = 2 * (nelec - 1) / integral
        pairden_dd[:,:,2] = restore_normalization(pairden_dd[:,:,2], multiplier)
        pairden_du[:,:,2] = restore_normalization(pairden_du[:,:,2], multiplier)
        pairden_ud[:,:,2] = restore_normalization(pairden_ud[:,:,2], multiplier)
        pairden_uu[:,:,2] = restore_normalization(pairden_uu[:,:,2], multiplier)
        
        # Normalization of dt and ut
        zz = pairden_dt[:,:,2] + pairden_ut[:,:,2]
        integral = np.trapz(np.trapz(zz, dx=dx, axis=0), dx=dy, axis=0)
        multiplier = 2 * (nelec - 1) / integral
        pairden_dt[:,:,2] = restore_normalization(pairden_dt[:,:,2], multiplier)
        pairden_ut[:,:,2] = restore_normalization(pairden_ut[:,:,2], multiplier)
        
        # Total Electron Density
        zz = (pairden_dt[:,:,2] + pairden_ut[:,:,2]) / 2
        pairden_t = np.stack([xx, yy, zz], axis=2)
        
        # Spin density
        zz = pairden_uu[:,:,2] + pairden_dd[:,:,2] \
             - (pairden_ud[:,:,2] + pairden_du[:,:,2])
        zz /= 2
        pairden_s = np.stack([xx, yy, zz], axis=2)
    else:
        zz = pairden_ud[:,:,2] + pairden_uu[:,:,2]
        integral = np.trapz(np.trapz(zz, dx=dx, axis=0), dx=dy, axis=0)
        multiplier = (nelec - 1) / integral
        pairden_ud[:,:,2] = restore_normalization(pairden_ud[:,:,2], multiplier)
        pairden_uu[:,:,2] = restore_normalization(pairden_uu[:,:,2], multiplier)

        # Normalization of ut
        zz = pairden_ut[:,:,2]
        integral = np.trapz(np.trapz(zz, dx=dx, axis=0), dx=dy, axis=0)
        multiplier = (nelec - 1) / integral
        pairden_dt[:,:,2] = restore_normalization(pairden_dt[:,:,2], multiplier)
        pairden_ut[:,:,2] = restore_normalization(pairden_ut[:,:,2], multiplier)
        
        # Total Electron Density
        zz = pairden_ut[:,:,2]
        pairden_t = np.stack([xx, yy, zz], axis=2)
        
        # Spin density
        zz = pairden_uu[:,:,2] - pairden_ud[:,:,2]
        pairden_s = np.stack([xx, yy, zz], axis=2)

    return [
        pairden_dd, pairden_dt, pairden_du,
        pairden_ud, pairden_ut, pairden_uu,
        pairden_t, pairden_s, xfix
    ]
    
    
def restore_normalization(density, multiplier):
    with np.errstate(divide='ignore', invalid='ignore'):
        density *= multiplier
        density = \
            np.nan_to_num(density, copy=True, nan=0.0, posinf=0.0, neginf=0.0)

    return density

