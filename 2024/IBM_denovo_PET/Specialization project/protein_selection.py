# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 16:14:23 2023

@author: Ingrid Moa Bolstad
"""
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import math

#---------Open file------------------------
abspath = os.path.abspath(__file__)
dirn = os.path.dirname(abspath)
os.chdir(dirn)
file=open(r"Analysis_results.txt", 'r')
lines=file.readlines()
file.close()

#-----Extract MolProbity scores and Ramachandran outliers---------
num_of_designs=len(lines)/8
names=[]
molprobity_scores=[]
molprob_names_and_scores=[]
for i in range(0,len(lines), 8):
    name=(lines[i].split(' ')[0])
    names.append(name)
    score=('%s.%s'%(re.findall(r'\d+', lines[i+7])[0],re.findall(r'\d+', lines[i+7])[1]))
    molprobity_scores.append(float(score))
    d={'%s'%name: score}
    molprob_names_and_scores.append(d)
    
rama_outliers=[]
rama_names_and_scores=[]
for i in range(0,len(lines), 8):
    name=(lines[i].split(' ')[0])
    names.append(name)
    score=float(lines[i][-10:-5].strip())
    rama_outliers.append(float(score))
    d={'%s'%name: score}
    rama_names_and_scores.append(d)
    
#clean the names in the namelist
names = [name.split('_')[0] for name in names]

#sort proteins by lowest percent outliers:
sorted_rama = sorted(rama_names_and_scores, key=lambda x: float(list(x.values())[0]))

#sort proteins by lowest molpriboty score:
sorted_molprob = sorted(molprob_names_and_scores, key=lambda x: float(list(x.values())[0]))


#-------------Write to file-----------------
#Ramachandran outliers
file = open("./Data/Ramachandran_outliers_sorted.txt", 'w')
for i in range(len(sorted_rama)):
    file.write(str(sorted_rama[i])+'\n')
file.close()

#Molprobity scores
file = open("./Data/Molprobity_scores_sorted.txt", 'w')
for i in range(len(sorted_molprob)):
    file.write(str(sorted_molprob[i])+'\n')
file.close()

#------------Plot values----------------------
ref_score=1.04 #reference LCC-ICCG MolProbity Score

def barplot(sorted_data, molprob=False, ramachandran=False): #set either molprob OR ramachandran to be True
    keys_list = [list(entry.keys())[0] for entry in sorted_data]
    values_list = [list(entry.values())[0] for entry in sorted_data]
    keys_B=[]
    keys_C=[]
    keys_A=[]
    values_B=[]
    values_C=[]
    values_A=[]
    for i in range(len(keys_list)):
        if 'B' in keys_list[i]:
            keys_B.append(keys_list[i])
            values_B.append(float(values_list[i]))
        if 'C' in keys_list[i]:
            keys_C.append(keys_list[i])
            values_C.append(float(values_list[i]))
        if 'A' in keys_list[i]:
            keys_A.append(keys_list[i])
            values_A.append(float(values_list[i]))
    if molprob:        
        plt.yticks(np.linspace(0, math.ceil(float(values_list[-1])), num=6))
        plt.ylim(0, float(values_list[-1])+1.7) # y-axis limits 
        plt.bar('ref', ref_score, width=0.7, color='#925FC9', label='LCC-ICCG', align='center')
        plt.bar(keys_A, values_A, width=0.7, color='#E5767C', label='A', align='center')
        plt.bar(keys_B, values_B, width=0.7, color='#6FCC71', label='B', align='center')
        plt.bar(keys_C, values_C, width=0.7, color='#1FC4C2', label='C', align='center')
        plt.xlabel('Protein design')
        plt.ylabel('MolProbity Score')
        plt.legend(loc='upper left')
        plt.xticks([])
        plt.savefig('./Figures/MolProbity_scores_plot.png',  dpi=300, bbox_inches='tight')
        plt.show()
    if ramachandran:
        plt.yticks(np.linspace(0, math.ceil(float(values_list[-1])+1), num=6))
        plt.ylim(0, float(values_list[-1])+0.2) # y-axis limits 
        plt.bar(keys_A, values_A, width=0.7, color='#E5767C', label='A', align='center')
        plt.bar(keys_B, values_B, width=0.7, color='#6FCC71', label='B', align='center')
        plt.bar(keys_C, values_C, width=0.7, color='#1FC4C2', label='C', align='center')
        plt.xlabel('Protein design')
        plt.ylabel('Ramachandran outliers [%]')
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xticks([])
        plt.savefig('./Figures/Ramachandran_outliers_plot.png',  dpi=300, bbox_inches='tight')
        plt.show()
        
barplot(sorted_molprob, molprob=True)
barplot(sorted_rama, ramachandran=True)

