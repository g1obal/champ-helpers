"""
Atomic Units Converter

Converts quantities;
SI units to atomic units
Atomic units to SI units

Bohr radius:
    a_0 = (4 * pi * varepsilon_0 * hbar^2) / (m_e * e^2) 
        = 5.291772109 x 10^−11 m

Hartree energy:
    E_h = m_e * (e^2 / (4 * pi * varepsilon_0 * hbar))^2 
        = 27.211386245988 eV
        = 4.3597447222071 x 10^-18 J

Effective electron mass ratio:
    m_r = m_eff / m_e

Dielectric constant:
    kappa = varepsilon / varepsilon_0

Effective Bohr radius:
    a_0_eff = a_0 * (kappa / m_r)
    
Effective Hartree energy:
    E_h_eff = E_h * (m_r / kappa**2)

** The electronvolt (eV) is not an SI unit however it is more practical to use 
   eV for our purposes. For this reason and for readability, we also convert 
   Joule to eV before atomic unit conversions. 1 J = 6.241509074461 x 10^18 eV.    

Author: Gokhan Oztarhan
Created date: 17/09/2020
Last modified: 18/02/2023
"""

LUNITS = {
    'm': 1,
    'cm': 1e-2,
    'mm': 1e-3,
    'um': 1e-6,
    'nm': 1e-9,
    'A': 1e-10,
    'pm': 1e-12,
    'fm': 1e-15,
}

EUNITS = {
    'J': 6.241509074461e18,
    'GeV': 1e9,
    'MeV': 1e6,
    'keV': 1e3,
    'eV': 1,
    'meV': 1e-3,
}

# Bohr radius
a_0 = 5.291772109e-11

# Hartree energy
E_h = 27.211386245988


class AUConverter():
    """Atomic Units Converter"""
    def __init__(self, m_r=1, kappa=1):
        # Effective Bohr radius
        self.a_0_eff = a_0 * (kappa / m_r)

        # Effective Hartree energy
        self.E_h_eff = E_h * (m_r / kappa**2)
        
    def length_to_au(self, length, unit):
        # Convert to meters, and convert to au
        return length * LUNITS[unit] / self.a_0_eff
        
    def energy_to_au(self, energy, unit):
        # Convert to eV, and convert to au
        return energy * EUNITS[unit] / self.E_h_eff
        
    def length_to_SI(self, length_au, unit):
        # Convert to meters, and convert to given units
        return length_au * self.a_0_eff / LUNITS[unit]
        
    def energy_to_SI(self, energy_au, unit):
        # Convert to eV, and convert to given units
        return energy_au * self.E_h_eff / EUNITS[unit]


