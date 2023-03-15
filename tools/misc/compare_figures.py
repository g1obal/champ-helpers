"""
Author: Gokhan Oztarhan
Created date: 15/03/2023
Last modified: 15/03/2023
"""

import os

import numpy as np
from matplotlib.image import imread


root_dir = 'density_plots'
root_dir_old = 'density_plots_old'

diff = 0
for root, dirs, files in sorted(os.walk(root_dir)):
    root_old = root.replace(root_dir, root_dir_old)
    
    for fname in files:
        img = imread(os.path.join(root, fname))
        img_old = imread(os.path.join(root_old, fname))
        
        diff += np.sqrt((img - img_old)**2).sum()
        
print(diff)
