# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 19:13:52 2023

@author: Ingrid Moa Bolstad
"""

import mmtbx.model
from mmtbx.validation import molprobity
import iotbx.pdb
import os
from io import StringIO
import sys

#---------------------------------------------------------------
#Uncomment to run ONE pdb at a time:

# pdb_inp = iotbx.pdb.input(file_name="./6tht.pdb") 
# def molprobity_analysis(pdb):
#     pdb_inp = iotbx.pdb.input(file_name='%s' %pdb)

# model = mmtbx.model.manager(model_input = pdb_inp)
# a=molprobity.molprobity(model)
# molprobity_score=a.molprobity_score()


# buffer = StringIO()
# sys.stderr = buffer
# a.show_summary(out=sys.stderr)
# print_output = buffer.getvalue()
# sys.stderr = sys.__stderr__


# file = open("6THT_results.txt", 'a')
# file.write(print_output)
# file.close()

#-------------------------------------------
#To run ALL pdb files:
    
path = r"./RFdiffusion_output" #path to folder with all the output from RFdiffusion
dir_list = os.listdir(path) #all folders

pdb_paths=[] #paths to pdb files will be loaded into this list
namelist=[] #names for all designs
for folder in dir_list:
    newpath=path+'/'+folder+'/outputs/'
    dir_list=os.listdir(newpath)
    name=dir_list[0]
    namelist.append(name)
    newpath=newpath+name+'/'+'best.pdb'
    pdb_paths.append(newpath)

pdb_input=[]
for i in range(len(pdb_paths)):
    pdb_input.append(iotbx.pdb.input(file_name='%s'%pdb_paths[i])) #file is the path of the pdb file


def molprobity_analysis(pdb_path):
    pdb=iotbx.pdb.input(file_name='%s'%pdb_path)
    model = mmtbx.model.manager(model_input = pdb)
    a=molprobity.molprobity(model)
    
    buffer = StringIO()
    sys.stderr = buffer
    a.show_summary(out=sys.stderr)
    print_output = buffer.getvalue()
    sys.stderr = sys.__stderr__
    
    return print_output

file = open("Analysis_results.txt", 'a')
for i in range(len(namelist)):
    file.write(namelist[i])
    analysis=molprobity_analysis(pdb_paths[i])
    file.write(analysis)
file.close()











