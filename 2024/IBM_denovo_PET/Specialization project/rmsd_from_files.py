# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:11:40 2023

@author: Ingrid Moa Bolstad
"""
#import libraries
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import os
import pickle
import warnings
warnings.filterwarnings('ignore')#suppress some MDAnalysis warnings about PSF files

#--------Functions-------------
def find_activesite(c): #c is in the format c=str('74-74/A90-100/94-94/A160-170/46-46/A185-195/16-16/A205-215/16-16/A237-247/25-25')
    cs = c.split('/')
    inserts = []
    sites = []
    for i in cs:
        css = i.split('-')
        if css[0][0] == 'A':
            sites.append([int(css[0][1:]), int(css[1])])
        else:
            inserts.append(int(css[0]))
    #as1 is active site 1, same for as2 and as3
    as1=1+inserts[0] + (sites[0][1]-sites[0][0]+1) + inserts[1] + 165-sites[1][0]
    
    as2=1+inserts[0] + (sites[0][1]-sites[0][0]+1) + inserts[1] + (sites[1][1]-sites[1][0]+1) + inserts[2] + \
    (sites[2][1]-sites[2][0]+1) + inserts[3] + 210-sites[3][0]
    
    as3=1+inserts[0] + (sites[0][1]-sites[0][0]+1) + inserts[1] + (sites[1][1]-sites[1][0]+1) + inserts[2] + \
    (sites[2][1]-sites[2][0]+1) + inserts[3] + (sites[3][1]-sites[3][0]+1) + inserts[4] + 242-sites[4][0]
    cat_triad=[as1, as2, as3]
    return(cat_triad)

def find_bindingsite(c): #c is in the format c=str('74-74/A90-100/94-94/A160-170/46-46/A185-195/16-16/A205-215/16-16/A237-247/25-25')
    cs = c.split('/')
    inserts = []
    sites = []
    for i in cs:
        css = i.split('-')
        if css[0][0] == 'A':
            sites.append([int(css[0][1:]), int(css[1])])
        else:
            inserts.append(int(css[0]))
   #bs1 is bindingsite 1, same for bs2 and bs3
    bs1=1+inserts[0] + 95-sites[0][0]
    bs2=1+inserts[0] + (sites[0][1]-sites[0][0]+1) + inserts[1] + 166-sites[1][0] 
    bs3=1+inserts[0] + (sites[0][1]-sites[0][0]+1) + inserts[1] + (sites[1][1]-sites[1][0]+1) + inserts[2] + \
    + 190-sites[2][0]
    bindingsite=[bs1, bs2, bs3]
    return(bindingsite)

    
def rmsd(design_dir, activeres): #takes in the directory of the pdb file and a list of the catalytic triad and returns rmsd
    s=activeres[0]
    d=activeres[1]
    h=activeres[2]
    u = mda.Universe(r"%s" %design_dir)  
    ref = mda.Universe(r"C:\Users\Bruker\OneDrive - NTNU\femtis\Prosjektoppgave\6tht_mhet (3).pdb")
    return rms.rmsd(u.select_atoms('resid %s or resid %s or resid %s' %(s,d,h)).positions,  # coordinates to align
             ref.select_atoms('resid 165 or resid 210 or resid 242').positions,  # reference coordinates
             center=True,  # subtract the center of geometry
             superposition=True)  # superimpose coordinates

def rmsd_bindingsite(design_dir, bindingres): #takes in the directory of the pdb file and a list of three binding residues and returns rmsd
    s=bindingres[0]
    d=bindingres[1]
    h=bindingres[2]
    u = mda.Universe(r"%s" %design_dir)  
    ref = mda.Universe(r"C:\Users\Bruker\OneDrive - NTNU\femtis\Prosjektoppgave\6tht_mhet (3).pdb") 
    return rms.rmsd(u.select_atoms('resid %s or resid %s or resid %s' %(s,d,h)).positions,  # coordinates to align
             ref.select_atoms('resid 95 or resid 166 or resid 190').positions,  # reference coordinates
             center=True,  # subtract the center of geometry
             superposition=True)  # superimpose coordinates


#----------------------------Load all files----------------------------------------
abspath = os.path.abspath(__file__)
dirn = os.path.dirname(abspath)
os.chdir(dirn)
# Get the list of all files and directories
path = r"./Data/RFdiffusion_output" #all 
dir_list = os.listdir(path)

dir1=[] #paths to trb files will be loaded into this list
dir2=[] #paths to pdb files will be loaded into this list
namelist=[] #names for all designs
for folder in dir_list:
    newpath=path+'/'+folder+'/outputs/'
    dir_list2=os.listdir(newpath)
    name=dir_list2[0]
    namelist.append(name)
    newpath2=newpath+name+'/'+'best.pdb'
    dir2.append(newpath2)
    a=os.listdir('%s' % newpath)
    for file in a:
        if '.trb' in file:
            dir1.append('%s'% newpath + '%s' %file)

#clean the names in the namelist
namelist = [name.split('_')[0] for name in namelist]

#open string of residues from .trb files and find active sites:
con=[]
activesites=[]
for d in dir1:
    with open('%s' %d, 'rb') as f:
        data=pickle.load(f)  
        c=str(data['config']['contigmap']['contigs']).replace("['","").replace("']","")
        con.append(c) 
    activesites.append(find_activesite(c)) 

#same for binding sites:
bindingsites=[]
for d in dir1:
    with open('%s' %d, 'rb') as f:
        data=pickle.load(f)  
        c=str(data['config']['contigmap']['contigs']).replace("['","").replace("']","") 
    bindingsites.append(find_bindingsite(c))

    
#Find rmsd of all pdb files and sort them from low to high:
def find_rmsd(dirlist, activeres, names, catalytic=False, binding=False): #Specify either catalytic OR binding to be True to find catalytic/binding rmsd
    results=[]
    for i in range(len(dirlist)):
        if catalytic:
            r=rmsd(dirlist[i], activeres[i])
        if binding:
            r=rmsd_bindingsite(dirlist[i], activeres[i])
        result={'%s'%names[i]: r}
        results.append(result)
        sorted_results = sorted(results, key=lambda x: float(list(x.values())[0]))
    return sorted_results #sorted_results is a list of dictionaries with name as key and rmsd as value

def all_activesites(dirlist, activeres, names):
    nameandsites=[]
    for i in range(len(dirlist)):
        n={'%s'%names[i]: activeres[i]}
        nameandsites.append(n)

    
#Write sorted rmsd to file:
cat_rmsd=find_rmsd(dir2, activesites, namelist, catalytic=True) #catalytic triad
binding_rmsd=find_rmsd(dir2, bindingsites, namelist, binding=True) #binding sites
file = open("./Data/rmsd_catalytic.txt", 'w')
file.write(str(cat_rmsd))
file = open("./Data/rmsd_binding.txt", 'w')
file.write(str(binding_rmsd))
file.close()

#Find rmsd for AlphaFold prediction of reference:
LCC_AF_dir=r"./Data/LCC-ICCG_AFprediction.pdb"
LCC_catres=(130, 175, 207) #The indices for the functional residues are skewed -35 residues as the reference (6THT) starts at residue 36
LCC_bindres=(60, 131, 155) #The indices for the functional residues are skewed -35 residues as the reference (6THT) starts at residue 36
cat_LCC_AF_rmsd=rmsd(LCC_AF_dir, LCC_catres)
binding_LCC_AF_rmsd=rmsd_bindingsite(LCC_AF_dir, LCC_bindres)
print(cat_LCC_AF_rmsd, binding_LCC_AF_rmsd)


#---------------------Plot results------------------------------------------------------
def rmsd_barplot(rmsd_list):
    keys_list = [list(entry.keys())[0] for entry in rmsd_list]
    values_list = [list(entry.values())[0] for entry in rmsd_list]
    keys_B=[]
    keys_C=[]
    keys_A=[]
    values_B=[]
    values_C=[]
    values_A=[]
    for i in range(len(keys_list)):
        if 'B' in keys_list[i]:
            keys_B.append(keys_list[i])
            values_B.append(values_list[i])
        if 'C' in keys_list[i]:
            keys_C.append(keys_list[i])
            values_C.append(values_list[i])
        if 'A' in keys_list[i]:
            keys_A.append(keys_list[i])
            values_A.append(values_list[i])
    plt.bar(keys_A, values_A, width=0.7, color='#E5767C', label='A', align='center')
    plt.bar(keys_B, values_B, width=0.7, color='#6FCC71', label='B', align='center')
    plt.bar(keys_C, values_C, width=0.7, color='#1FC4C2', label='C', align='center')
    plt.xlabel('Protein design')
    #plt.title('Catalytic triad deviation relative to LCC-ICCG')
    plt.legend()
    plt.xticks([])
    if rmsd_list==cat_rmsd:
        plt.ylabel('Catalytic triad RMSD [Å]')
        plt.savefig('./Figures/rmsd_catalytic.png',  dpi=300, bbox_inches='tight')
    if rmsd_list==binding_rmsd:
        plt.ylabel('Binding residues RMSD [Å]')
        plt.savefig('./Figures/rmsd_binding.png',  dpi=300, bbox_inches='tight')
    plt.show()

rmsd_barplot(cat_rmsd)
rmsd_barplot(binding_rmsd)





    
