# -*- coding: utf-8 -*-
"""
Created on Wed May 16 19:52:46 2018
Dynamic Flux Analysis & Constrained Flux Regulation

@author: scott
"""
from __future__ import print_function
import sys, os
import cobra, cobra.test, json
from escher import Builder
from gurobipy import *
from os.path import join
from cobra.io.mat import model_to_pymatbridge
from pymatbridge import Matlab

# if you're going to be working with both Matlab and Python:
mlab = Matlab()
mlab.start()
results = mlab.run_code('a=1;')
model.solver = 'gurobi'

# redirect towards the datafiles
data_dir = r'C:\\Users\\scott\\Google Drive\\Work\\UM\\Research\\Sriram\\Projects\\Lyssiotis_Proteomics'

# we need to use a recon1 model with an objective function:
recon1_model = cobra.io.load_matlab_model(join(data_dir, "model_human_duarte.mat"))
model_to_pymatbridge(recon1_model, variable_name="metabolicmodel")

recon1_model.metabolicmodel.S
#Palsson_core = cobra.io.load_matlab_model(join(data_dir, "Core_Model_Palsson.mat"), variable_name="core_genecomb")

def flux_activity_coeff(model, timecourse_metabolomics, sheetname='Sheet1', kappa=None, kappa2=None, genedelflag=None, rxndelflag=None, normquantile=None):
    '''docstring
    flux_activity_coeff 
    
    '''
    # default params
    if kappa is None:
        kappa = 1
    if kappa2 is None:
        kappa2 = 1E-3
    if sheetname is None:
        sheetname = 'Sheet1'
    if genedelflag is None:
        genedelflag = 0
    if rxndelflag is None:
        rxndelflag = 0
    if normquantile is None:
        normquantile = 0
    
    filetype = os.path.splitext(model)[-1]
    # MATLAB functionality
    if filetype is '.mat':
        # test model
        model_test = cobra.test.create_test_model('model)
        
def contrain_flux_regulation(model1, onreactions, offreactions, kappa, rho, epsilon, mode, epsilon2, minfluxflag):
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
