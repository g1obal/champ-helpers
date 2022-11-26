"""
Hartree Atomic Units

Converts quantities;
SI units to atomic units
Atomic units to SI units

Author: Gokhan Oztarhan
Created date: 17/09/2020
Last modified: 02/08/2022
"""

import numpy as np


def convert_length(length, unit, m_r, kappa):
    # Convert length to meters
    if unit == 'm':
        pass
    elif unit == 'nm':
        length = length * 1e-9
    elif unit == 'A':
        length = length * 1e-10
       
    # a_0 = ( 4 * pi * varepsilon_0 * hbar^2 ) / ( m_e * e^2  ) 
    #     = 5.291772109 x 10^−11 m
    a_0 = 5.291772109e-11
    
    # m_r = m_eff / m_e, kappa = varepsilon / varepsilon_0
    a_0_eff = a_0 * (kappa / m_r)
    
    length_au = length / a_0_eff
    
    return length_au


def convert_energy(energy, unit, m_r, kappa):
    # E_h = m * ( ( e^2 ) / ( 4 * pi * varepsilon_0 * hbar ) )^2 
    #     = 27.211386245988 eV 
    #     = 4.3597447222071 x 10^-18 J
    if 'eV' in unit:
        E_h = 27.211386245988
    elif unit == 'J':
        E_h = 4.3597447222071e-18
        
    # Convert unit to eV
    if unit == 'eV':
        pass
    elif unit == 'meV':
        energy = energy * 1e-3
    
    # m_r = m_eff / m_e, kappa = varepsilon / varepsilon_0
    E_h_eff = E_h * ( m_r / kappa**2 )
    
    energy_au = energy / E_h_eff
    
    return energy_au

    
def convert_length_inverse(length_au, unit, m_r, kappa):       
    # a_0 = ( 4 * pi * varepsilon_0 * hbar^2 ) / ( m_e * e^2  ) 
    #     = 5.291772109 x 10^−11 m
    a_0 = 5.291772109e-11
    
    # m_r = m_eff / m_e, kappa = varepsilon / varepsilon_0
    a_0_eff = a_0 * (kappa / m_r)
    
    length = length_au * a_0_eff

    # Convert length to meters
    if unit == 'm':
        pass
    elif unit == 'nm':
        length = length / 1e-9
    elif unit == 'A':
        length = length / 1e-10
    
    return length
    

def convert_energy_inverse(energy_au, unit, m_r, kappa):
    # E_h = m * ( ( e^2 ) / ( 4 * pi * varepsilon_0 * hbar ) )^2 
    #     = 27.211386245988 eV 
    #     = 4.3597447222071 x 10^-18 J
    if 'eV' in unit:
        E_h = 27.211386245988
    elif ( unit == 'J' ):
        E_h = 4.3597447222071e-18
    
    # m_r = m_eff / m_e, kappa = varepsilon / varepsilon_0
    E_h_eff = E_h * ( m_r / kappa**2 )
    
    energy = energy_au * E_h_eff 
    
    # Convert unit to eV
    if unit == 'eV':
        pass
    elif unit == 'meV':
        energy = energy / 1e-3
    
    return energy

