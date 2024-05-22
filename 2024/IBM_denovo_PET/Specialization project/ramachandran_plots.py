# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 14:12:13 2023

@author: Ingrid Moa Bolstad

The function number_of_res and codelines 45-50 was written using ChatGPT 3.5
"""

#import libraries
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import warnings
warnings.filterwarnings('ignore')# suppress MDAnalysis warnings about PSF files
import os

#------Load files----------
abspath = os.path.abspath(__file__)
dirn = os.path.dirname(abspath)
os.chdir(dirn)
#pdb files to plot Ramachandran distributions for:
pathB18 = r"./Data/RFdiffusion_output/B18.result/outputs/B18/best.pdb"
pathB22 = r"./Data/RFdiffusion_output/B22.result/outputs/B22/best.pdb"
paths=[pathB18, pathB22]
names=['B18', 'B22']


#------------------------------Functions---------------------

def number_of_res(pdb):
    with open(r'%s' %pdb, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                columns = line.split()
                final_res = columns[5][-1] #finds the number of residues
                return final_res
            
def ramachandran(design_dir):
    u = mda.Universe(r"%s" %design_dir) 
    r = u.select_atoms("protein")
    R = Ramachandran(r).run()
    return R

fig, axes = plt.subplots(ncols=len(paths), figsize=(4*len(paths), 4))
for i, path in enumerate(paths):
    R = ramachandran(path)
    ax = axes[i] if len(paths) > 1 else axes
    R.plot(ax=ax, color='000000', marker='o', ref='True')
    ax.set_title(f"{names[i]}")
    
plt.tight_layout()
plt.savefig('./Figures/ramachandran_B18_B22.png',  dpi=300, bbox_inches='tight')
plt.show()






 