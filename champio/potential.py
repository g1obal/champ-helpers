"""
Potential functions for generate_input

Author: Gokhan Oztarhan
Created date: 10/06/2019
Last modified: 30/01/2022
"""

import numpy as np


def cyldot(x, v, s, rho):
    return (v/2.0) * (np.tanh(s * (x + rho)) - np.tanh(s * (x - rho)))


def radius_cyldot(k, s, rho):
    return np.arccosh((k - 1) * np.cosh(2 * s * rho) + k) / (2 * s)
    

def gndot(x, v0, rho, s):
    return v0 * np.exp(-((x * x) / (rho * rho))**s)
