# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 23:31:00 2023

@author: Ingrid Moa Bolstad
The function AA_sequence was co-written by ChatGPT 3.5
"""

import os
from Bio import PDB

#----------------------------Load all files----------------------------------------
abspath = os.path.abspath(__file__)
dirn = os.path.dirname(abspath)
os.chdir(dirn)
# Get the list of all files and directories
path = r"./Data/RFdiffusion_output" #all 
dir_list = os.listdir(path)
LCC_path= r"./Data/6tht_mhet.pdb"

dirs=[] #paths to pdb files will be loaded into this list
namelist=[] #names for all designs
for folder in dir_list:
    newpath=path+'/'+folder+'/outputs/'
    dir_list=os.listdir(newpath)
    name=dir_list[0]
    namelist.append(name)
    newpath=newpath+name+'/'+'best.pdb'
    dirs.append(newpath)
    
#------Find AA sequences----------

def AA_sequence(pdb_path):
    structure = PDB.PDBParser().get_structure('protein', pdb_path)
    amino_acid_sequence = ""
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    amino_acid_sequence += PDB.Polypeptide.three_to_one(residue.get_resname())
    return amino_acid_sequence


sequences=[] #index corresponds with namelist
for d in dirs: 
    sequences.append(AA_sequence(d))

#----------Find average sequence lengths----------
protein_lengthA=[]
protein_lengthB=[]
protein_lengthC=[]
for i in range(len(sequences)):
    if 'A' in namelist[i]:
        protein_lengthA.append(len(sequences[i]))
    if 'B' in namelist[i]:
        protein_lengthB.append(len(sequences[i]))
    if 'C' in namelist[i]:
        protein_lengthC.append(len(sequences[i]))
    
#find average lengths for A, B and C:
avg_A=round(sum(protein_lengthA)/len(protein_lengthA))
avg_B=round(sum(protein_lengthB)/len(protein_lengthB))
avg_C=round(sum(protein_lengthC)/len(protein_lengthC))

print(avg_A, avg_B, avg_C)
print(len(protein_lengthB))

#-----------Find longest and shortest sequences--------
lengths=[]
for seq in sequences:
    lengths.append(len(seq))
combined_data = list(zip(lengths, namelist))
sorted_lengths = sorted(combined_data, key=lambda x: x[0], reverse=True)

print(sorted_lengths)